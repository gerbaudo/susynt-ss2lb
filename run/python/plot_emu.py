#!/bin/env python

# make basic plots for the emu selection
#
# davide.gerbaudo@gmail.com
# Sep 2014

import collections
import datetime
import math
import optparse
import os
import pprint

import dataset
from rootUtils import (drawAtlasLabel
                       ,dummyHisto
                       ,getBinContents
                       ,getMinMax
                       ,graphWithPoissonError
                       ,increaseAxisFont
                       ,topRightLegend
                       ,importRoot
                       ,integralAndError
                       ,setAtlasStyle
                       ,topRightLabel
                       ,writeObjectsToFile
                       )
r = importRoot()
from utils import (first
                   ,getCommandOutput
                   ,mkdirIfNeeded
                   ,filterWithRegexp
                   ,remove_duplicates
                   ,sortedAs
                   )

import utils
from kin import addTlv, computeCollinearMassLepTau

susyntutils = utils.import_susyntutils()
r = susyntutils.import_root()
susyntutils.load_packages()

from CutflowTable import CutflowTable

import systUtils

usage="""
This code is used either (1) to fill the histos, or (2) to make plots
and tables. The output of (1) is used as input of (2).

Required inputs: trees produced with `submitJobs.py --seltuple` and
with `--matrix-prediction`.

Example usage ('fill' mode):
%prog \\
 --input-other  out/selection_tuple/Jul_26 \\
 --input-fake out/matrix_prediction/Jul_26/ \\
 --output-dir out/plot_emu/Jul_26/histos \\
 --verbose \\
 2>&1 | tee log/plot_emu/Jul_26/fill.log

Example usage ('plot' mode):
%prog \\
 --input-dir out/plot_emu/Jul_26/histos \\
 --output-dir out/plot_emu/Jul_26/ \\
 --verbose \\
 2>&1 | tee log/plot_emu/Jul_26/plot.log"""

def main() :
    parser = optparse.OptionParser(usage=usage)
    parser.add_option('-b', '--batch',  action='store_true', help='submit to batch (used in fill mode)')
    parser.add_option('-g', '--group', help='group to be processed (used only in fill mode)')
    parser.add_option('-f', '--input-fake', help='location of fake trees')
    parser.add_option('-O', '--input-other', help='location other trees')
    parser.add_option('-i', '--input-dir')
    parser.add_option('-o', '--output-dir')
    parser.add_option('--samples-dir', default='samples/',
                      help='directory with the list of samples; default ./samples/')
    parser.add_option('-v', '--verbose', action='store_true', default=False)
    parser.add_option('--debug', action='store_true', default=False)

    (opts, args) = parser.parse_args()
    inOtherSpecified, inDirSpecified = opts.input_other!=None, opts.input_dir!=None
    eitherMode = inOtherSpecified != inDirSpecified
    if not eitherMode : parser.error("Run either in 'fill' or 'plot' mode")
    mode = 'fill' if inOtherSpecified else 'plot' if inDirSpecified else None
    requiredOptions = (['input_fake', 'input_other', 'output_dir'] if mode=='fill'
                       else ['input_dir', 'output_dir'])
    allOptions = [x.dest for x in parser._get_all_options()[1:]]
    def optIsNotSpecified(o) : return not hasattr(opts, o) or getattr(opts,o) is None
    if any(optIsNotSpecified(o) for o in requiredOptions):
        parser.error('Missing required option\n'
                     +'\n'.join(["%s : %s"%(o, getattr(opts, o)) for o in requiredOptions]))
    if opts.verbose:
        print ('\nUsing the following options:\n'
               +'\n'.join("%s : %s"%(o, str(getattr(opts, o))) for o in allOptions))

    if opts.debug : dataset.Dataset.verbose_parsing = True
    groups = dataset.DatasetGroup.build_groups_from_files_in_dir(opts.samples_dir)
    groups.append(first([g for g in groups if g.is_data]).clone_data_as_fake())
    if opts.group : groups = [g for g in groups if g.name==opts.group]
    print '\n'.join("group {0} : {1} samples".format(g.name, len(g.datasets)) for g in groups)

    if   mode=='fill' : runFill(opts, groups)
    elif mode=='plot' : runPlot(opts, groups)

def runFill(opts, groups) :
    batchMode    = opts.batch
    inputFakeDir = opts.input_fake
    inputGenDir  = opts.input_other
    outputDir    = opts.output_dir
    verbose      = opts.verbose

    if verbose : print "filling histos"
    mkdirIfNeeded(outputDir)
    systematics = ['NOM']
    if verbose : print "about to loop over these systematics:\n %s"%str(systematics)
    if batchMode:
        for syst in systematics: # todo : probably not needed, leave it in for now
            for group in groups :
                submit_batch_fill_job_per_group(group, opts)
    else:
        selections = regions_to_plot()
        variables = variables_to_plot()
        group_names = [g.name for g in groups]
        counters_all_groups = bookCounters(group_names, selections)
        histos_all_groups = bookHistos(variables, group_names, selections)
        for group in groups:
            tree_name = 'hlfv_tuple'
            chain = r.TChain(tree_name)
            input_dir = opts.input_fake if group.name=='fake' else opts.input_other
            for dataset in group.datasets:
                chain.Add(os.path.join(input_dir, dataset.name+'.root'))
            if opts.verbose:
                print "{0} : {1} entries from {2} samples".format(group.name,
                                                                  chain.GetEntries(),
                                                                  len(group.datasets))
            histos = histos_all_groups[group.name]
            counters = counters_all_groups[group.name]
            for iEntry, event in enumerate(chain):
                if iEntry>10000:break
                run_num = event.pars.runNumber
                evt_num = event.pars.eventNumber
                weight =  event.pars.weight
                l0 = addTlv(event.l0)
                l1 = addTlv(event.l1)
                met = addTlv(event.met)
                l0_is_el, l0_is_mu = l0.isEl, l0.isMu
                l1_is_el, l1_is_mu = l1.isEl, l1.isMu
                is_emu = int(l0_is_el and l1_is_mu)
                is_mue = int(l0_is_mu and l1_is_el)
                is_same_sign = int((l0.charge * l1.charge)>0)
                is_opp_sign  = not is_same_sign
                l0_pt, l1_pt = l0.p4.Pt(), l1.p4.Pt()
                dphi_l0_met = abs(l0.p4.DeltaPhi(met.p4))
                dphi_l1_met = abs(l1.p4.DeltaPhi(met.p4))
                dphi_l0_l1 = abs(l0.p4.DeltaPhi(l1.p4))
                dpt_l0_l1 = l0.p4.Pt()-l1.p4.Pt()
                m_coll = computeCollinearMassLepTau(l0.p4, l1.p4, met.p4)
                for sel in selections:
                    pass_sel = eval(selection_formulas()[sel])
                    if not pass_sel : continue

                    histos[sel]['onebin'].Fill(1.0, weight)
                    histos[sel]['pt0'].Fill(l0.p4.Pt(), weight)
                    histos[sel]['pt1'].Fill(l1.p4.Pt(), weight)
                    histos[sel]['mcoll'].Fill(m_coll, weight)
                    counters[sel] += (weight) # if passSels[sel] else 0.0)

            for v in ['onebin', 'pt0', 'pt1']:
                for sel in selections:
                    h = histos[sel][v]
                    print "{0}: integral {1}, entries {2}".format(h.GetName(), h.Integral(), h.GetEntries())
        plotting_groups = dict([(g.name, systUtils.Group(g.name)) for g in groups])
        saveHistos(plotting_groups, histos_all_groups, outputDir, opts.verbose)
        # print counters

def runPlot(opts, groups) :
    inputDir     = opts.input_dir
    outputDir    = opts.output_dir
    sysOption    = 'NOM' #opts.syst
    excludedSyst = None #opts.exclude
    verbose      = opts.verbose
    mkdirIfNeeded(outputDir)
    buildTotBkg = systUtils.buildTotBackgroundHisto
    buildStat = systUtils.buildStatisticalErrorBand
    buildSyst = systUtils.buildSystematicErrorBand

    selections = regions_to_plot()
    variables = variables_to_plot()
    plot_groups = [systUtils.Group(g.name) for g in groups]
    for group in plot_groups :
        group.setHistosDir(inputDir)
        group.exploreAvailableSystematics(verbose)
        group.filterAndDropSystematics(sysOption, excludedSyst, verbose)

    mkdirIfNeeded(outputDir)
    systematics = ['NOM']
    anySys = sysOption==None
    if sysOption=='fake'   or anySys : systematics += systUtils.fakeSystVariations()
    if sysOption=='object' or anySys : systematics += systUtils.mcObjectVariations()
    if sysOption=='weight' or anySys : systematics += systUtils.mcWeightVariations()
    if sysOption and sysOption.count(','):
        systematics = [s for s in systUtils.getAllVariations() if s in sysOption.split(',')]
    elif sysOption in systUtils.getAllVariations() : systematics = [sysOption]
    if not anySys and len(systematics)==1 and sysOption!='NOM' : raise ValueError("Invalid syst %s"%str(sysOption))
    if excludedSyst : systematics = [s for s in systematics if s not in filterWithRegexp(systematics, excludedSyst)]
    if verbose : print "using the following systematics : %s"%str(systematics)

    fakeSystematics = [s for s in systematics if s in systUtils.fakeSystVariations()]
    mcSystematics = [s for s in systematics if s in systUtils.mcObjectVariations() + systUtils.mcWeightVariations()]

    simBkgs = [g for g in plot_groups if g.isMcBkg]
    data = findByName(plot_groups, 'data')
    fake = findByName(plot_groups, 'fake')
    signal = findByName(plot_groups, 'signal')
    print 'names_stacked_groups to be improved'
    names_stacked_groups = [g.name for g in groups if g.name not in ['data', 'signal']]
    for sel in selections :
        if verbose : print '-- plotting ',sel
        for var in variables :
            if verbose : print '---- plotting ',var
            for g in plot_groups : g.setSystNominal()
            nominalHistoData    = data.getHistogram(variable=var, selection=sel, cacheIt=True)
            nominalHistoSign    = signal.getHistogram(variable=var, selection=sel, cacheIt=True)
            nominalHistoFakeBkg = fake.getHistogram(variable=var, selection=sel, cacheIt=True)
            nominalHistosSimBkg = dict([(g.name, g.getHistogram(variable=var, selection=sel, cacheIt=True))
                                        for g in simBkgs])
            nominalHistosBkg    = dict([('fake', nominalHistoFakeBkg)] +
                                       [(g, h) for g, h in nominalHistosSimBkg.iteritems()])
            nominalHistoTotBkg  = buildTotBkg(histoFakeBkg=nominalHistoFakeBkg,
                                              histosSimBkgs=nominalHistosSimBkg)
            statErrBand = buildStat(nominalHistoTotBkg)
            systErrBand = buildSyst(fake=fake, simBkgs=simBkgs, variable=var, selection=sel,
                                    fakeVariations=fakeSystematics, mcVariations=mcSystematics,
                                    verbose=verbose)

            plotHistos(histoData=nominalHistoData,# histoSignal=nominalHistoSign,
                       histoTotBkg=nominalHistoTotBkg, histosBkg=nominalHistosBkg,
                       statErrBand=statErrBand, systErrBand=systErrBand,
                       stack_order=names_stacked_groups,
                       topLabel=sel,
                       canvasName=(sel+'_'+var), outdir=outputDir, verbose=verbose)
    for group in plot_groups :
        group.printVariationsSummary()

def submit_batch_fill_job_per_group(group, opts):
    verbose = opts.verbose

    group_name = group.name if hasattr(group, 'name') else group
    newOptions  = " --input-other %s" % opts.input_other
    newOptions += " --input-fake %s" % opts.input_fake
    newOptions += " --output-dir %s" % opts.output_dir
    newOptions += " --group %s" % group_name
    newOptions += (" --verbose " if opts.verbose else '')
    template = 'batch/templates/plot_emu.sh'
    log_dir = mkdirIfNeeded('log/plot_emu')
    script_dir = mkdirIfNeeded('batch/plot_emu')
    script_name = os.path.join(script_dir, group_name+'.sh')
    script_file = open(script_name, 'w')
    script_file.write(open(template).read()
                      .replace('%(opt)s', newOptions)
                      .replace('%(logfile)s', log_dir+'/'+group_name+'.log')
                      .replace('%(jobname)s', group_name))
    script_file.close()
    cmd = "sbatch %s"%script_name
    if verbose : print cmd
    out = getCommandOutput(cmd)
    if verbose : print out['stdout']
    if out['stderr'] : print  out['stderr']

def countAndFillHistos(samplesPerGroup={}, syst='', verbose=False, outdir='./') :

    selections = allRegions()
    variables = variablesToPlot()

    mcGroups, fakeGroups = mcDatasetids().keys(), ['fake']
    objVariations, weightVariations, fakeVariations = systUtils.mcObjectVariations(), systUtils.mcWeightVariations(), systUtils.fakeSystVariations()
    def groupIsRelevantForSys(g, s) :
        isRelevant = (s=='NOM' or (g in mcGroups and s in objVariations+weightVariations) or (g in fakeGroups and s in fakeVariations))
        if verbose and not isRelevant : print "skipping %s for %s"%(g, s)
        return isRelevant
    def dropIrrelevantGroupsForThisSys(groups, sys) : return dict((g, samples) for g, samples in groups.iteritems() if groupIsRelevantForSys(g, syst))
    def dropSamplesWithoutTree(samples) : return [s for s in samples if s.hasInputHftTree(msg='Warning! ')]
    def dropGroupsWithoutSamples(groups) : return dict((g, samples) for g, samples in groups.iteritems() if len(samples))
    samplesPerGroup = dropIrrelevantGroupsForThisSys(samplesPerGroup, syst)
    samplesPerGroup = dict((g, dropSamplesWithoutTree(samples)) for g, samples in samplesPerGroup.iteritems())
    samplesPerGroup = dropGroupsWithoutSamples(samplesPerGroup)

    groups = samplesPerGroup.keys()
    counters = bookCounters(groups, selections)
    histos = bookHistos(variables, groups, selections)
    for group, samplesGroup in samplesPerGroup.iteritems() :
        logLine = "---->"
        if verbose : print 1*' ',group
        histosGroup = histos  [group]
        countsGroup = counters[group]
        for sample in samplesGroup :
            if verbose : logLine +=" %s"%sample.name
            fillAndCount(histosGroup, countsGroup, sample, blind=False)
        if verbose : print logLine
    if verbose : print 'done'
    return counters, histos

def printCounters(counters):
    countTotalBkg(counters)
    blindGroups   = [g for g in counters.keys() if g!='data']
    unblindGroups = [g for g in counters.keys()]
    tableSr  = CutflowTable(samples=blindGroups,   selections=signalRegions(), countsSampleSel=counters)
    tablePre = CutflowTable(samples=blindGroups,   selections=controlRegions(), countsSampleSel=counters)
    tableBld = CutflowTable(samples=unblindGroups, selections=blindRegions(), countsSampleSel=counters)
    for table in [tableSr, tablePre, tableBld] : table.nDecimal = 6
    print 4*'-',' sig regions ',4*'-'
    print tableSr.csv()
    print 4*'-',' pre regions ',4*'-'
    print tablePre.csv()
    print 4*'-',' blind regions ',4*'-'
    print tableBld.csv()
#___________________________________________________________
def allGroups(noData=False, noSignal=True) :
    return ([k for k in mcDatasetids().keys() if k!='signal' or not noSignal]
            + ([] if noData else ['data'])
            + ['fake']
            )

def selection_formulas(sel=None):
    pt_req = 'l0_pt>45.0 and l1_pt>12.0'
    common_req = (pt_req+' and '+
                  'dphi_l1_met<0.7 and dphi_l0_l1>2.3 and '+
                  'dpt_l0_l1>7.0 and dphi_l0_met>2.5')
    formulas = {
        'pre_emu' : 'is_emu and '+pt_req,
        'pre_mue' : 'is_mue and '+pt_req,
        'pre_emu_mue' : '(is_emu or is_mue) and '+pt_req,
        'sr_emu' : 'l0_is_el and l1_is_mu and '+common_req,
        'sr_mue' : 'l0_is_mu and l1_is_el and '+common_req,
        'sr_emu_mue' : '(is_emu or is_mue) and '+common_req,
        }
    formulas = dict([(k+'_'+ssos, v+' and '+ssos_expr)
                     for k, v in formulas.iteritems()
                     for ssos, ssos_expr in [('ss', 'is_same_sign'), ('os', 'is_opp_sign')]])
    # symmetric selection
    pt_sym_req = 'l0_pt>20.0 and l1_pt>20.0'
    for lf, lf_expr in [('emu', 'is_emu'), ('mue', 'is_mue'), ('emu_mue', '(is_emu or is_mue)')]:
        for ssos, ssos_expr in [('ss', 'is_same_sign'), ('os', 'is_opp_sign')]:
            formulas['sym_'+lf+'_'+ssos] = pt_sym_req+' and '+lf_expr+' and '+ssos_expr
    # validation region used by Matt in the 2L paper, see sec6.4 ATL-COM-PHYS-2012-1808
    formulas_vrss_btag = 'num_b_jets==1 and et_miss_rel>50.0 and abs(m_ll-91.2)>10.0 if is_ee else True) and ((m_ll<90.0 or m_ll>120) if is_mumu else True)'
    # formulas['vrss_btag'] = formulas_vrss_btag
    return formulas[sel] if sel else formulas


def fillAndCount(histos, counters, sample, blind=True) :
    group    = sample.group
    filename = sample.filenameHftTree
    treename = sample.hftTreename
    file = r.TFile.Open(filename)
    tree = file.Get(treename)
    selections = allRegions()
    selWeights = dict((s, r.TTreeFormula(s, selectionFormulas(s), tree)) for s in selections)
    weightFormula = r.TTreeFormula('weightFormula', sample.weightLeafname, tree)
    l1 = r.TLorentzVector()
    l2 = r.TLorentzVector()
    met = r.TLorentzVector()
    for iEvent, event in enumerate(tree) :
        weight = weightFormula.EvalInstance()
        passSels = dict((s, selWeights[s].EvalInstance()) for s in selections)
        for sel in selections : counters[sel] += (weight if passSels[sel] else 0.0)
        for sel in selections :
            fillHisto = passSels[sel]
            if blind and sample.isData :
                if sel in signalRegions() : fillHisto = False
                else : fillHisto = passSels[blindRegionFromAnyRegion(sel)] and not passSels[signalRegionFromAnyRegion(sel)]
            oneJet = event.L2nCentralLightJets==1
            mev2gev = 1.0e-3
            mljj = mev2gev*(event.mlj if oneJet else event.mljj)
            l1.SetPtEtaPhiM(event.lept1Pt*mev2gev, event.lept1Eta, event.lept1Phi, 0.0) # massless here is good enough
            l2.SetPtEtaPhiM(event.lept2Pt*mev2gev, event.lept2Eta, event.lept2Phi, 0.0)
            met.SetPtEtaPhiM(event.met*mev2gev,               0.0, event.metPhi,   0.0)
            ll = l1+l2
            ptll = ll.Pt()
            mll = ll.M()
            l1IsMu = event.lept1Flav==1
            l2IsMu = event.lept2Flav==1
            dphil0met = abs(l1.DeltaPhi(met)) if l1.Pt()>l2.Pt() else abs(l2.DeltaPhi(met))
            if fillHisto :
                histos[sel]['mll'   ].Fill(mll, weight)
                histos[sel]['mljj'  ].Fill(mljj, weight)
                histos[sel]['ptll'  ].Fill(ptll, weight)
                histos[sel]['onebin'].Fill(1.0,  weight)
                histos[sel]['dphil0met'].Fill(dphil0met, weight)
                if l1IsMu or l2IsMu:
                    dphimumet = abs(l1.DeltaPhi(met)) if l1IsMu else abs(l2.DeltaPhi(met))
                    histos[sel]['dphimumet'].Fill(dphimumet, weight)
            # checks
            if (True and fillHisto
                and sel in signalRegions()
                and (sample.isData or sample.isFake)) :
                channel = 'ee' if event.isEE else 'mm' if event.isMUMU else 'em'
                dataOrFake = 'data' if sample.isData else 'fake' if sample.isFake else 'other'
                print "ev %d run %d channel %s sel %s sample %s weight %f"%(event.runNumber, event.eventNumber, channel, sel, dataOrFake, weight)
    file.Close()

def dataSampleNames() :
    return ["period%(period)s.physics_%(stream)s"%{'period':p, 'stream':s}
            for p in ['A','B','C','D','E','G','H','I','J','L']
            for s in ['Egamma','Muons']]
def mcDatasetids() :
    print 'todo: now it is an attribute Dataset.dsid parsed from the dataset name'

def allSamplesAllGroups() :
    Sample = systUtils.Sample
    asg = dict( [(group, [Sample(groupname=group, name=dsid) for dsid in dsids]) for group, dsids in mcDatasetids().iteritems()]
               +[('data', [Sample(groupname='data', name=s) for s in dataSampleNames()])]
               +[('fake', [Sample(groupname='fake', name=s) for s in dataSampleNames()])])
    return asg
def allGroups() :
    return [systUtils.Group(g) for g in mcDatasetids().keys()+['data']+['fake']]

def stackedGroups(groups) :
    return [g for g in allSamplesAllGroups().keys() if g not in ['data', 'signal']]

def variablesToPlot() :
    return ['onebin','mljj', 'ptll', 'mll', 'dphil0met', 'dphimumet']
    return ['pt0','pt1','mll','mtmin','mtmax','mtllmet','ht','metrel','dphill','detall',
            'mt2j','mljj','dphijj','detajj']
def histoName(sample, selection, variable) : return "h_%s_%s_%s"%(variable, sample, selection)
def bookHistos(variables, samples, selections) :
    "book a dict of histograms with keys [sample][selection][var]"
    def histo(variable, sam, sel) :
        twopi = +2.0*math.pi
        mljjLab = 'm_{lj}' if '1j' in sel else 'm_{ljj}'
        h = None
        if   v=='onebin'  : h = r.TH1F(histoName(sam, sel, 'onebin' ), ';; entries',                             1, 0.5,   1.5)
        elif v=='pt0'     : h = r.TH1F(histoName(sam, sel, 'pt0'    ), ';p_{T,l0} [GeV]; entries/bin',          12, 0.0, 240.0)
        elif v=='pt1'     : h = r.TH1F(histoName(sam, sel, 'pt1'    ), ';p_{T,l1} [GeV]; entries/bin',          12, 0.0, 240.0)
        elif v=='mcoll'   : h = r.TH1F(histoName(sam, sel, 'mcoll'  ), ';m_{coll,l0,l1} [GeV]; entries/bin',    12, 0.0, 240.0)
        elif v=='mll'     : h = r.TH1F(histoName(sam, sel, 'mll'    ), ';m_{l0,l1} [GeV]; entries/bin',         12, 0.0, 240.0)
        elif v=='ptll'    : h = r.TH1F(histoName(sam, sel, 'ptll'   ), ';p_{T,l0+l1} [GeV]; entries/bin',       12, 0.0, 240.0)
        elif v=='dphil0met': h= r.TH1F(histoName(sam, sel, 'dphil0met'),';#Delta#phi(l0, met) [rad]; entries/bin',  10, 0.0, twopi)
        elif v=='dphil1met': h= r.TH1F(histoName(sam, sel, 'dphil1met'),';#Delta#phi(l1, met) [rad]; entries/bin',  10, 0.0, twopi)
        else : print "unknown variable %s"%v
        h.Sumw2()
        h.SetDirectory(0)
        return h
    return dict([(sam, dict([(sel, dict([(v, histo(v, sam, sel)) for v in variables]))
                         for sel in selections]))
                 for sam in samples])
def bookCounters(samples, selections) :
    "book a dict of counters with keys [sample][selection]"
    return dict((s, dict((sel, 0.0) for sel in selections)) for s in samples)
def countTotalBkg(counters={'sample' : {'sel':0.0}}) :
    backgrounds = [g for g in counters.keys() if g!='signal' and g!='data']
    selections = first(counters).keys()
    counters['totBkg'] = dict((s, sum(counters[b][s] for b in backgrounds)) for s in selections)
def getGroupColor(g) :
    oldColors = [('data', r.kBlack), ('diboson',r.kSpring+2), ('higgs',r.kAzure-4),
                 ('signal',r.kMagenta), ('top', r.kRed+1), ('zjets', r.kOrange-2),
                 ('fake',r.kGray)]
    newColors = [] #[('signal',r.kMagenta), ('WW',r.kAzure-9), ('Higgs',r.kYellow-9)]
    colors = dict((g,c) for g,c in  oldColors + newColors)
    return colors[g]

def regions_to_plot():
    return selection_formulas().keys()
    # return ['emu_ss', 'emu_os']
def variables_to_plot():
    return ['onebin', 'pt0', 'pt1', 'mcoll']

def plotHistos(histoData=None, histoSignal=None, histoTotBkg=None, histosBkg={},
               statErrBand=None, systErrBand=None, # these are TGraphAsymmErrors
               canvasName='canvas', outdir='./', verbose=False,
               stack_order=[],
               topLabel='',
               drawStatErr=False, drawSystErr=False,
               drawYieldAndError=False) :
    "Note: blinding can be required for only a subrange of the histo, so it is taken care of when filling"
    setAtlasStyle()
    padMaster = histoData
    if verbose : print "plotting ",padMaster.GetName()
    can = r.TCanvas(canvasName, padMaster.GetTitle(), 800, 600)
    can.cd()
    can._hists = [padMaster]
    padMaster.Draw('axis')
    can.Update() # necessary to fool root's dumb object ownership of the stack
    stack = r.THStack('stack_'+padMaster.GetName(), '')
    r.SetOwnership(stack, False)
    can._hists.append(stack)
    leg = topRightLegend(can, 0.225, 0.325)
    can._leg = leg
    leg.SetBorderSize(0)
    leg._reversedEntries = []
    def integralWou(h):
        "Integral with underflow and overflow"
        return h.Integral(0, h.GetNbinsX()+1)
    for group, histo in sortedAs(histosBkg, stack_order) :
        histo.SetFillColor(getGroupColor(group))
        histo.SetLineWidth(2)
        histo.SetLineColor(r.kBlack)
        stack.Add(histo)
        can._hists.append(histo)
        leg._reversedEntries.append((histo, "{0}: {1:.2f}".format(group, integralWou(histo)), 'F'))
    leg._reversedEntries.append((dummyHisto(), "{0}, {1:.2f}".format('bkg', sum([integralWou(h) for h in stack.GetHists()])), 'l'))
    leg._reversedEntries.append((histoData, "{0}, {1:.2f}".format('data', integralWou(histoData)), 'p'))
    for h, g, o in leg._reversedEntries[::-1] : leg.AddEntry(h, g, o) # stack goes b-t, legend goes t-b
    stack.Draw('hist same')
    histoData.SetMarkerStyle(r.kFullCircle)
    histoData.SetLineWidth(2)
    dataGraph = graphWithPoissonError(histoData)
    dataGraph.Draw('same p')
    if histoSignal :
        histoSignal.SetLineColor(getGroupColor('signal'))
        histoSignal.SetLineWidth(2)
        histoSignal.Draw('histo same')
        leg.AddEntry(histoSignal, '(m_{C1},m_{N1})=(130, 0)GeV', 'l')
    if statErrBand and drawStatErr :
        statErrBand.SetFillStyle(3006)
        statErrBand.Draw('E2 same')
        leg.AddEntry(statErrBand, 'stat', 'f')
    if systErrBand and drawSystErr :
        systErrBand.SetFillStyle(3007)
        systErrBand.Draw('E2 same')
        leg.AddEntry(systErrBand, 'syst', 'f')
    totErrBand = systUtils.addErrorBandsInQuadrature(statErrBand, systErrBand)
    if totErrBand :
        totErrBand.Draw('E2 same')
        totErrBand.SetFillStyle(3005)
        leg.AddEntry(totErrBand, 'stat+syst', 'f')
    leg.Draw('same')
    can.Update()
    tex = r.TLatex()
    tex.SetTextSize(0.5 * tex.GetTextSize())
    tex.SetNDC(True)
    label  = "%s tot bkg : "%(can.GetName())
    label += "%.3f #pm %.3f (stat)"%(integralAndError(histoTotBkg))
    if systErrBand :
        sysUp, sysDo = systUtils.totalUpDownVariation(systErrBand)
        label += "#pm #splitline{%.3f}{%.3f} (syst)"%(sysUp, sysDo)
    if drawYieldAndError :
        tex.DrawLatex(0.10, 0.95, label)
        can.SetTopMargin(2.0*can.GetTopMargin())
    drawAtlasLabel(can, xpos=0.125, align=13)
    if topLabel : topRightLabel(can, topLabel, ypos=1.0)
    yMin, yMax = getMinMax([histoData, dataGraph, histoTotBkg, histoSignal, totErrBand])
    padMaster.SetMinimum(0.0)
    padMaster.SetMaximum(1.1 * yMax)
    increaseAxisFont(padMaster.GetXaxis())
    increaseAxisFont(padMaster.GetYaxis())
    can.RedrawAxis()
    can.Update() # force stack to create padMaster
    for ext in ['png','eps'] : can.SaveAs(outdir+'/'+can.GetName()+'.'+ext)

def saveHistos(samplesPerGroup={}, histosPerGroup={}, outdir='./', verbose=False) :
    for groupname, histosThisGroup in histosPerGroup.iteritems() :
        group = samplesPerGroup[groupname].setHistosDir(outdir)
        outFilename = group.filenameHisto
        writeObjectsToFile(outFilename, histosThisGroup, verbose)

if __name__=='__main__' :
    main()
