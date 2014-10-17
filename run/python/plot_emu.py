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
    parser.add_option('-s', '--syst', help="variations to process (default all)."
                      " Give a comma-sep list or say 'weight', 'object', or 'fake'")
    parser.add_option('-l', '--list-systematics', action='store_true', default=False,
                      help='list what is already in output_dir')
    parser.add_option('-L', '--list-all-systematics', action='store_true', default=False,
                      help='list all possible systematics')
    parser.add_option('-e', '--exclude', help="skip some systematics, example 'EL_FR_.*'")
    parser.add_option('-v', '--verbose', action='store_true', default=False)
    parser.add_option('--debug', action='store_true', default=False)

    (opts, args) = parser.parse_args()
    if opts.list_all_systematics :
        print "All systematics:\n\t%s"%'\n\t'.join(systUtils.getAllVariations())
        return
    if opts.list_systematics :
        print listExistingSyst(opts.input_dir)
        return

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

    if   mode=='fill' : runFill(opts)
    elif mode=='plot' : runPlot(opts)

def runFill(opts) :
    batchMode    = opts.batch
    inputFakeDir = opts.input_fake
    inputGenDir  = opts.input_other
    outputDir    = opts.output_dir
    verbose      = opts.verbose

    if opts.debug : dataset.Dataset.verbose_parsing = True
    groups = dataset.DatasetGroup.build_groups_from_files_in_dir(opts.samples_dir)
    groups.append(first([g for g in groups if g.is_data]).clone_data_as_fake())
    if opts.group : groups = [g for g in groups if g.name==opts.group]
    print '\n'.join("group {0} : {1} samples".format(g.name, len(g.datasets)) for g in groups)

    if verbose : print "filling histos"
    mkdirIfNeeded(outputDir)
    systematics = get_list_of_syst_to_fill(opts)
    if verbose : print "about to loop over these systematics:\n %s"%str(systematics)
    if batchMode:
        for group in groups:
            for systematic in systematics:
                if systUtils.Group(group.name).isNeededForSys(systematic):
                    opts.syst = systematic
                    submit_batch_fill_job_per_group(group, opts)
    else:
        for group in groups:
            systematics = [s for s in systematics if systUtils.Group(group.name).isNeededForSys(s)]
            if not systematics : print "warning, empty syst list. You should have at least the nominal"
            for systematic in systematics:
                # note to self: here you will want to use a modified Sample.setHftInputDir
                # for now we just have the fake syst that are in the nominal tree
                tree_name = 'hlfv_tuple'
                chain = r.TChain(tree_name)
                input_dir = opts.input_fake if group.name=='fake' else opts.input_other
                for ds in group.datasets:
                    chain.Add(os.path.join(input_dir, ds.name+'.root'))
                if opts.verbose:
                    print "{0} : {1} entries from {2} samples".format(group.name,
                                                                      chain.GetEntries(),
                                                                      len(group.datasets))
                counters, histos = count_and_fill(chain=chain, sample=group.name,
                                                  syst=systematic, verbose=verbose)
                out_filename = systUtils.Group(group.name).setSyst(systematic).setHistosDir(outputDir).filenameHisto
                print 'out_filename: ',out_filename
                writeObjectsToFile(out_filename, histos, verbose)
        # print counters

def runPlot(opts) :
    inputDir     = opts.input_dir
    outputDir    = opts.output_dir
    verbose      = opts.verbose
    mkdirIfNeeded(outputDir)
    buildTotBkg = systUtils.buildTotBackgroundHisto
    buildStat = systUtils.buildStatisticalErrorBand
    buildSyst = systUtils.buildSystematicErrorBand
    selections = regions_to_plot()
    variables = variables_to_plot()

    groups = dataset.DatasetGroup.build_groups_from_files_in_dir(opts.samples_dir)
    groups.append(first([g for g in groups if g.is_data]).clone_data_as_fake())
    plot_groups = [systUtils.Group(g.name) for g in groups]
    for group in plot_groups :
        group.setHistosDir(inputDir)
        group.exploreAvailableSystematics(verbose)
        group.filterAndDropSystematics(opts.syst, opts.exclude, verbose)
    available_systematics = sorted(list(set([s for g in plot_groups for s in g.systematics])))
    systematics_to_use = get_list_of_syst_to_fill(opts)
    systematics = [s for s in systematics_to_use if s in available_systematics]
    if verbose :
        print "using the following systematics : {0}".format(systematics)
        print "missing the following systematics : {0}".format([s for s in systematics_to_use if s not in available_systematics])
    fakeSystematics = [s for s in systematics if s in systUtils.fakeSystVariations()]
    mcSystematics = [s for s in systematics if s in systUtils.mcObjectVariations() + systUtils.mcWeightVariations()]

    mkdirIfNeeded(outputDir)
    findByName = systUtils.findByName
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
    systematic = opts.syst if opts.syst else None
    newOptions  = " --input-other %s" % opts.input_other
    newOptions += " --input-fake %s" % opts.input_fake
    newOptions += " --output-dir %s" % opts.output_dir
    newOptions += " --group %s" % group_name
    newOptions += (" --verbose " if opts.verbose else '')
    newOptions += (" --syst {0}".format(opts.syst) if opts.syst else '')
    newOptions += (" --exclude {0}".format(opts.exclude) if opts.exclude else '')
    template = 'batch/templates/plot_emu.sh'
    log_dir = mkdirIfNeeded('log/plot_emu')
    script_dir = mkdirIfNeeded('batch/plot_emu')
    script_name = os.path.join(script_dir, group_name+("_{0}".format(systematic) if systematic else '')+'.sh')
    log_name = log_dir+'/'+group_name+("_{0}".format(systematic) if systematic else '')+'.log'
    script_file = open(script_name, 'w')
    script_file.write(open(template).read()
                      .replace('%(opt)s', newOptions)
                      .replace('%(logfile)s', log_name)
                      .replace('%(jobname)s', group_name))
    script_file.close()
    cmd = "sbatch %s"%script_name
    if verbose : print cmd
    out = getCommandOutput(cmd)
    if verbose : print out['stdout']
    if out['stderr'] : print  out['stderr']

#-------------------
def count_and_fill(chain, sample='', syst='', verbose=False):
    """
    count and fill for one sample (or group), one syst.
    """
    sysGroup = systUtils.Group(sample).setSyst(syst)
    selections = regions_to_plot()
    counters = book_counters(selections)
    histos = book_histograms(sample_name=sample, variables=variables_to_plot(),
                             systematics=[syst], selections=selections
                             )[syst]
    weight_expr = 'event.pars.weight'
    weight_expr = sysGroup.weightLeafname
    print 'weight_expr: ',weight_expr
    for iEntry, event in enumerate(chain):
        weight = eval(weight_expr)
        run_num = event.pars.runNumber
        evt_num = event.pars.eventNumber
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
    if verbose:
        for v in ['onebin', 'pt0', 'pt1']:
            for sel in selections:
                h = histos[sel][v]
                print "{0}: integral {1}, entries {2}".format(h.GetName(), h.Integral(), h.GetEntries())
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
#___________________________________________________________
def variablesToPlot() :
    return ['onebin','mljj', 'ptll', 'mll', 'dphil0met', 'dphimumet']
def book_histograms(sample_name='', variables=[], systematics=[], selections=[]) :
    "book a dict of histograms with keys [systematics][selection][var]"
    histoName = systUtils.BaseSampleGroup.histoname
    def histo(variable, sam, sys, sel) :
        twopi = +2.0*math.pi
        mljjLab = 'm_{lj}' if '1j' in sel else 'm_{ljj}'
        h = None
        if   v=='onebin'  : h = r.TH1F(histoName(sam, sys, sel, 'onebin' ), ';; entries',                             1, 0.5,   1.5)
        elif v=='pt0'     : h = r.TH1F(histoName(sam, sys, sel, 'pt0'    ), ';p_{T,l0} [GeV]; entries/bin',          12, 0.0, 240.0)
        elif v=='pt1'     : h = r.TH1F(histoName(sam, sys, sel, 'pt1'    ), ';p_{T,l1} [GeV]; entries/bin',          12, 0.0, 240.0)
        elif v=='mcoll'   : h = r.TH1F(histoName(sam, sys, sel, 'mcoll'  ), ';m_{coll,l0,l1} [GeV]; entries/bin',    12, 0.0, 240.0)
        elif v=='mll'     : h = r.TH1F(histoName(sam, sys, sel, 'mll'    ), ';m_{l0,l1} [GeV]; entries/bin',         12, 0.0, 240.0)
        elif v=='ptll'    : h = r.TH1F(histoName(sam, sys, sel, 'ptll'   ), ';p_{T,l0+l1} [GeV]; entries/bin',       12, 0.0, 240.0)
        elif v=='dphil0met': h= r.TH1F(histoName(sam, sys, sel, 'dphil0met'),';#Delta#phi(l0, met) [rad]; entries/bin',  10, 0.0, twopi)
        elif v=='dphil1met': h= r.TH1F(histoName(sam, sys, sel, 'dphil1met'),';#Delta#phi(l1, met) [rad]; entries/bin',  10, 0.0, twopi)
        else : print "unknown variable %s"%v
        h.Sumw2()
        h.SetDirectory(0)
        return h
    return dict([(sys,
                  dict([(sel,
                         dict([(v, histo(v, sample_name, sys, sel))
                               for v in variables]))
                        for sel in selections]))
                 for sys in systematics])
def book_counters(selections) :
    "book a dict of counters with keys [systematic][selection]"
    return dict((sel, 0.0) for sel in selections)
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

def listExistingSyst(dirname) :
    print "listing systematics from ",dirname
    print "...not implemented yet..."

def get_list_of_syst_to_fill(opts):
    systematics = ['NOM']
    sysOption    = opts.syst
    excludedSyst = opts.exclude
    anySys       = sysOption==None
    if sysOption=='fake'   or anySys : systematics += systUtils.fakeSystVariations()
    if sysOption=='object' or anySys : systematics += systUtils.mcObjectVariations()
    if sysOption=='weight' or anySys : systematics += systUtils.mcWeightVariations()
    if sysOption and sysOption.count(','):
        systematics = [s for s in systUtils.getAllVariations() if s in sysOption.split(',')]
    elif sysOption in systUtils.getAllVariations(): systematics = [sysOption]
    elif not anySys and len(systematics)==1 and sysOption!='NOM':
        raise ValueError("Invalid syst %s"%str(sysOption))
    if excludedSyst:
        systematics = [s for s in systematics if s not in filterWithRegexp(systematics, excludedSyst)]
    return systematics


if __name__=='__main__' :
    main()
