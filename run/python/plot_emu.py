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
import time
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
                       ,buildRatioHistogram
                       ,buildBotTopPads
                       ,getXrange
                       ,referenceLine
                       )
r = importRoot()
from utils import (first
                   ,getCommandOutput
                   ,mkdirIfNeeded
                   ,filterWithRegexp
                   ,remove_duplicates
                   ,sortedAs
                   )
import fakeUtils as fakeu
import utils
from kin import addTlv, computeCollinearMassLepTau, computeRazor, computeMt

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

You can also add '--batch --log-dir log/plot_emu/Jul_26'.

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
    parser.add_option('--log-dir', default='log/plot_emu', help='directory where the batch logs will be')
    parser.add_option('-e', '--exclude', help="skip some systematics, example 'EL_FR_.*'")
    parser.add_option('-T', '--tight-def', help='on-the-fly tight def, one of defs in fakeUtils.py: fakeu.lepIsTight_std, etc.')
    parser.add_option('-v', '--verbose', action='store_true', default=False)
    parser.add_option('--debug', action='store_true', default=False)
    parser.add_option('--unblind', action='store_true', default=False)
    parser.add_option('--require-tight-tight', action='store_true', default=False, help='fill histos only when both leps are tight')

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
    def optIsNotSpecified(o) : return not hasattr(opts, o) or getattr(opts,o) is None
    if any(optIsNotSpecified(o) for o in requiredOptions):
        parser.error('Missing required option\n'
                     +'\n'.join(["%s : %s"%(o, getattr(opts, o)) for o in requiredOptions]))
    if opts.verbose : utils.print_running_conditions(parser, opts)

    if   mode=='fill' : runFill(opts)
    elif mode=='plot' : runPlot(opts)

def runFill(opts) :
    batchMode    = opts.batch
    inputFakeDir = opts.input_fake
    inputGenDir  = opts.input_other
    outputDir    = opts.output_dir
    verbose      = opts.verbose
    debug        = opts.debug
    blinded      = not opts.unblind
    tightight    = opts.require_tight_tight

    if debug : dataset.Dataset.verbose_parsing = True
    groups = dataset.DatasetGroup.build_groups_from_files_in_dir(opts.samples_dir)
    # groups.append(dataset.DatasetGroup.build_qflip_from_simulated_samples(groups))
    groups.append(first([g for g in groups if g.is_data]).clone_data_as_fake())
    if opts.group : groups = [g for g in groups if g.name==opts.group]
    if verbose : print '\n'.join("group {0} : {1} samples".format(g.name, len(g.datasets)) for g in groups)
    if debug :
        print '\n'.join("group {0} : {1} samples: {2}".format(g.name,
                                                              len(g.datasets),
                                                              '\n\t'+'\n\t'.join(d.name for d in g.datasets))
                        for g in groups)
    if verbose : print "filling histos"
    # eval will take care of aborting on typos
    onthefly_tight_def = eval(opts.tight_def) if opts.tight_def else None
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
                                                  syst=systematic, verbose=verbose,
                                                  debug=debug, blinded=blinded,
                                                  onthefly_tight_def=onthefly_tight_def,
                                                  tightight=tightight)
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
    # groups.append(dataset.DatasetGroup.build_qflip_from_simulated_samples(groups)) # not ready yet
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
    newOptions += (" --unblind " if opts.unblind else '')
    newOptions += (" --verbose " if opts.verbose else '')
    newOptions += (" --debug " if opts.debug else '')
    newOptions += (" --syst {0}".format(opts.syst) if opts.syst else '')
    newOptions += (" --exclude {0}".format(opts.exclude) if opts.exclude else '')
    newOptions += (" --require-tight-tight " if opts.require_tight_tight else '')
    newOptions += (" --tight-def {0}".format(opts.tight_def) if opts.tight_def else '')
    newOptions += (" --samples-dir {0}".format(opts.samples_dir) if opts.samples_dir else '')

    print 'todo: re-implement submit_batch_fill_job_per_group (just mod group opt)'
    template = 'batch/templates/plot_emu.sh'
    log_dir = mkdirIfNeeded(opts.log_dir)
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
def count_and_fill(chain, sample='', syst='', verbose=False, debug=False, blinded=True,
                   onthefly_tight_def=None, tightight=False):
    """
    count and fill for one sample (or group), one syst.
    """
    sysGroup = systUtils.Group(sample).setSyst(syst)
    is_mc = systUtils.Group(sample).isMc
    is_data = systUtils.Group(sample).isData
    is_qflip_sample = dataset.DatasetGroup(sample).is_qflip
    selections = regions_to_plot()
    counters = book_counters(selections)
    histos = book_histograms(sample_name=sample, variables=variables_to_fill(),
                             systematics=[syst], selections=selections
                             )[syst]
    if is_qflip_sample : # for qflip, only fill ss histos
        selections = [s for s in selections if s.endswith('_ss')]
    weight_expr = 'event.pars.weight'
    weight_expr = sysGroup.weightLeafname
    qflip_expr = 'event.pars.qflipWeight'
    print 'weight_expr: ',weight_expr
    start_time = time.clock()
    num_processed_entries = 0
    fields_to_print = ['l0_pt', 'l1_pt', 'l0_eta', 'l1_eta',
                       'met_pt',
                       'm_ll', 'pt_ll', 'dpt_l0_l1',
                       'dphi_l0_met', 'dphi_l1_met', 'dphi_l0_l1',
                       'mt0', 'mt1',
                       'n_soft_jets',
                       'eta_csj0', 'phi_csj0', 'eta_csj1', 'phi_csj1']
    if debug : print ",".join(fields_to_print)
    for iEntry, event in enumerate(chain):
        run_num = event.pars.runNumber
        evt_num = event.pars.eventNumber
        l0 = addTlv(event.l0)
        l1 = addTlv(event.l1)
        met = addTlv(event.met)
        l0_is_el, l0_is_mu = l0.isEl, l0.isMu
        l1_is_el, l1_is_mu = l1.isEl, l1.isMu
        l0_is_t = onthefly_tight_def(l0) if onthefly_tight_def else l0.isTight
        l1_is_t = onthefly_tight_def(l1) if onthefly_tight_def else l1.isTight
        is_emu = int(l0_is_el and l1_is_mu)
        is_mue = int(l0_is_mu and l1_is_el)
        is_mumu = int(l0_is_mu and l1_is_mu)
        is_ee = int(l0_is_el and l1_is_el)
        is_same_sign = int((l0.charge * l1.charge)>0)
        is_opp_sign  = not is_same_sign
        is_qflippable = is_opp_sign and (l0_is_el or l1_is_el) and is_mc
        weight = eval(weight_expr)
        qflip_prob = eval(qflip_expr)
        # print "event : same sign {0}, opp_sign {1}, qflippable {2}, qflip_prob {3}".format(is_same_sign, is_opp_sign, is_qflippable, eval(qflip_expr))
        l0_pt, l1_pt = l0.p4.Pt(), l1.p4.Pt()
        d_pt0_pt1 = l0_pt - l1_pt
        l0_eta, l1_eta = l0.p4.Eta(), l1.p4.Eta()
        l0_phi, l1_phi = l0.p4.Phi(), l1.p4.Phi()
        met_pt = met.p4.Pt()
        m_ll = (l0.p4 + l1.p4).M()
        pt_ll = (l0.p4 + l1.p4).Pt()
        dphi_l0_met = abs(l0.p4.DeltaPhi(met.p4))
        dphi_l1_met = abs(l1.p4.DeltaPhi(met.p4))
        dphi_l0_l1 = abs(l0.p4.DeltaPhi(l1.p4))
        dpt_l0_l1 = l0.p4.Pt()-l1.p4.Pt()
        m_coll = computeCollinearMassLepTau(l0.p4, l1.p4, met.p4)
        mt0, mt1 = computeMt(l0.p4, met.p4), computeMt(l1.p4, met.p4)
        dphillbeta, mdr = computeRazor(l0.p4, l1.p4, met.p4)
        def jet_pt2(j) : return j.px*j.px+j.py*j.py
        n_cl_jets = sum(1 for j in event.jets if jet_pt2(j)>30.*30.)
        n_jets = n_cl_jets + event.pars.numFjets + event.pars.numBjets
        # n_jets = event.pars.numFjets + event.pars.numBjets
        soft_jets = [addTlv(j) for j in event.jets if jet_pt2(j)<30.**2]
        n_soft_jets = len(soft_jets)
        csj0 = first(sorted(soft_jets, key=lambda j : j.p4.DeltaR(l0.p4)))
        csj1 = first(sorted(soft_jets, key=lambda j : j.p4.DeltaR(l1.p4)))
        eta_csj0 = csj0.p4.Eta() if csj0 else -5.0
        phi_csj0 = csj0.p4.Phi() if csj0 else -5.0
        eta_csj1 = csj1.p4.Eta() if csj1 else -5.0
        phi_csj1 = csj1.p4.Phi() if csj1 else -5.0
        drl0csj  = csj0.p4.DeltaR(l0.p4) if csj0 else None
        drl1csj  = csj1.p4.DeltaR(l1.p4) if csj1 else None
        pass_sels = {}
        if tightight and not (l0_is_t and l1_is_t) : continue
        for sel in selections:
            pass_sel = eval(selection_formulas()[sel])
            pass_sels[sel] = pass_sel
            is_ss_sel = sel.endswith('_ss')
            as_qflip = is_qflippable and (is_opp_sign and is_ss_sel)
            if is_qflip_sample and not as_qflip : pass_sel = False
            if not is_qflip_sample and as_qflip : pass_sel = False
            if not pass_sel : continue
            # <isElectron 1> <isElectron 2> <isTight 1> <isTight 2> <pt 1> <pt 2> <eta 1> <eta 2>
            lltype = "{0}{1}".format('e' if l0_is_el else 'mu', 'e' if l1_is_el else 'mu')
            qqtype = "{0}{1}".format('T' if l0_is_t else 'L', 'T' if l1_is_t else 'L')
            if debug : print ','.join([str(eval(_)) for _ in fields_to_print])
            def fmt(b) : return '1' if b else '0'
            # --- begin dbg
            # print "event: {0:12s} {1} {2} {3} {4} {5} {6} {7} {8}".format(lltype+' '+qqtype, #+' '+sel,
            #                                                               fmt(l0_is_el), fmt(l1_is_el),
            #                                                               fmt(l0_is_t), fmt(l1_is_t),
            #                                                               l0_pt, l1_pt,
            #                                                               l0.p4.Eta(), l1.p4.Eta())
            # print "event: {0:12s} {1} {2} {3:.2f} {4:.2f}".format(lltype+' '+qqtype+' '+sel,
            #                                                       run_num, evt_num,
            #                                                       l0_pt, l1_pt)
            # --- end dbg
            fill_weight = (weight * qflip_prob) if as_qflip else weight
            h = histos[sel]
            h['onebin'   ].Fill(1.0, fill_weight)
            h['njets'    ].Fill(n_jets, fill_weight)
            h['pt0'      ].Fill(l0_pt, fill_weight)
            h['pt1'      ].Fill(l1_pt, fill_weight)
            h['d_pt0_pt1'].Fill(d_pt0_pt1, fill_weight)
            h['eta0'     ].Fill(l0_eta, fill_weight)
            h['eta1'     ].Fill(l1_eta, fill_weight)
            h['phi0'     ].Fill(l0_phi, fill_weight)
            h['phi1'     ].Fill(l1_phi, fill_weight)
            h['mll'      ].Fill(m_ll, fill_weight)
            h['ptll'     ].Fill(pt_ll, fill_weight)
            h['met'      ].Fill(met_pt, fill_weight)
            h['dphil0met'].Fill(dphi_l0_met, fill_weight)
            h['dphil1met'].Fill(dphi_l1_met, fill_weight)

            h['pt0_vs_pt1'      ].Fill(l1_pt, l0_pt, fill_weight)
            h['met_vs_pt1'      ].Fill(l1_pt, met.p4.Pt(), fill_weight)
            h['dphil0met_vs_pt1'].Fill(l1_pt, dphi_l1_met, fill_weight)
            h['dphil0met_vs_pt1'].Fill(l1_pt, dphi_l1_met, fill_weight)
            h['nsj'             ].Fill(n_soft_jets, fill_weight)
            if n_soft_jets:
                h['drl0csj'].Fill(drl0csj, fill_weight)
                h['drl1csj'].Fill(drl1csj, fill_weight)
            if is_data and (blinded and 100.0<m_coll and m_coll<150.0) : pass
            else :
                h['mcoll'].Fill(m_coll, fill_weight)
                h['mcoll_vs_pt1'].Fill(l1_pt, m_coll, fill_weight)
            counters[sel] += (fill_weight)
        # print ('e' if l0_is_el else 'm'),('e' if l1_is_el else 'm'),' : ',
        # print ' is_opp_sign: ',is_opp_sign,
        # print ' is_qflippable: ',is_qflippable,
        # print pass_sels
        num_processed_entries += 1
    end_time = time.clock()
    delta_time = end_time - start_time
    if verbose:
        print ("processed {0:d} entries ".format(num_processed_entries)
               +"in "+("{0:d} min ".format(int(delta_time/60)) if delta_time>60 else
                       "{0:.1f} s ".format(delta_time))
               +"({0:.1f} kHz)".format((num_processed_entries/delta_time) if delta_time else 1.0e9)
               )
    if verbose:
        for v in ['onebin']: #, 'pt0', 'pt1']:
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
def selection_formulas():
    pt_req = 'l0_pt>45.0 and l1_pt>12.0 '
    pt_req += ' and abs(l0_eta)<2.4 and abs(l1_eta)<2.4' # tmp to be dropped, fixed upstream in ee8b767
    common_req_sr = (pt_req+
                     ' and n_jets==0'+
                     ' and dphi_l1_met<0.7'+
                     ' and dphi_l0_l1>2.3 '+
                     ' and dphi_l0_met>2.5'+
                     ' and (l0_pt-l1_pt)>7.0')
    common_req_vr = (pt_req+
                     ' and n_jets==0'+
                     ' and (dphi_l1_met>0.7 or dphi_l0_l1<2.3 or dphi_l0_met<2.5 or (l0_pt-l1_pt)<7.0)')
    formulas = {
        'sr_emu' : 'l0_is_el and l1_is_mu and '+common_req_sr,
        'sr_mue' : 'l0_is_mu and l1_is_el and '+common_req_sr,
        'vr_emu' : 'l0_is_el and l1_is_mu and '+common_req_vr,
        'vr_mue' : 'l0_is_mu and l1_is_el and '+common_req_vr,
        'sr_emu_mue' : '(is_emu or is_mue) and '+common_req_sr,
        'sr_ee' : 'l0_is_el and l1_is_el and '+common_req_sr,
        'sr_mumu' : 'l0_is_mu and l1_is_mu and '+common_req_sr,
        'vr_emu_mue' : '(is_emu or is_mue) and '+common_req_vr,
        'vr_ee' : 'l0_is_el and l1_is_el and '+common_req_vr,
        'vr_mumu' : 'l0_is_mu and l1_is_mu and '+common_req_vr,
        # 'pre_emu' : 'is_emu and '+pt_req,
        # 'pre_mue' : 'is_mue and '+pt_req,
        # 'pre_emu_mue' : '(is_emu or is_mue) and '+pt_req,
        }
    os_expr = '(is_opp_sign)'
    ss_expr = '(is_same_sign or is_qflippable)'
    formulas = dict([(k+'_'+ssos, v+' and '+ssos_expr)
                     for k, v in formulas.iteritems()
                     for ssos, ssos_expr in [('ss', ss_expr), ('os', os_expr)]])
    # # symmetric selection
    # pt_sym_req = 'l0_pt>20.0 and l1_pt>20.0'
    # for lf, lf_expr in [('emu', 'is_emu'), ('mue', 'is_mue'), ('emu_mue', '(is_emu or is_mue)')]:
    #     for ssos, ssos_expr in [('ss', ss_expr), ('os', os_expr)]:
    #         formulas['sym_'+lf+'_'+ssos] = pt_sym_req+' and '+lf_expr+' and '+ssos_expr
    # validation region used by Matt in the 2L paper, see sec6.4 ATL-COM-PHYS-2012-1808
    # formulas_vrss_btag = 'num_b_jets==1 and et_miss_rel>50.0 and abs(m_ll-91.2)>10.0 if is_ee else True) and ((m_ll<90.0 or m_ll>120) if is_mumu else True)'
    # formulas['vrss_btag'] = formulas_vrss_btag

    formulas['vr_emu_razor_ss'] = '(is_emu or is_mue) and mdr>20.0 and '+ss_expr
    formulas['vr_ee_razor_ss'] = 'is_ee and mdr>20.0 and '+ss_expr
    formulas['vr_mumu_razor_ss'] = 'is_mumu and mdr>20.0 and '+ss_expr
    # fake extraction formulas; trig_match is already applied when making the ntuples
    formulas['ext_mumu_ss'] = 'is_mumu and is_same_sign'
    formulas['ext_emu_mue_ss'] = '(is_emu or is_mue) and is_same_sign'
    formulas['ext_emu_pt0_40_ss']= 'is_emu and is_same_sign and l0_pt>40.0'
    formulas['ext_mue_pt0_40_ss']= 'is_mue and is_same_sign and l0_pt>40.0'
    formulas['ext_mumu_pt0_40_ss']= 'is_mumu and is_same_sign and l0_pt>40.0'
    # dbg low pt
    formulas['sr_mue_os_low_pt1_15'] = (formulas['sr_mue_os']+' and l1_pt<15.0')
    return formulas
#___________________________________________________________
def book_histograms(sample_name='', variables=[], systematics=[], selections=[]) :
    "book a dict of histograms with keys [systematics][selection][var]"
    histoName = systUtils.BaseSampleGroup.histoname
    def histo(variable, sam, sys, sel) :
        twopi = +2.0*math.pi
        mljjLab = 'm_{lj}' if '1j' in sel else 'm_{ljj}'
        h = None
        if   v=='onebin'   : h = r.TH1F(histoName(sam, sys, sel, v), ';; entries',                               1, 0.5,   1.5)
        elif v=='njets'    : h = r.TH1F(histoName(sam, sys, sel, v), ';N_{jets}; entries',                      10,-0.5,   9.5)
        elif v=='pt0'      : h = r.TH1F(histoName(sam, sys, sel, v), ';p_{T,l0} [GeV]; entries/bin',            48, 0.0, 240.0)
        elif v=='pt1'      : h = r.TH1F(histoName(sam, sys, sel, v), ';p_{T,l1} [GeV]; entries/bin',            48, 0.0, 240.0)
        elif v=='d_pt0_pt1': h = r.TH1F(histoName(sam, sys, sel, v), ';p_{T,l0}-p_{T,l1} [GeV]; entries/bin',   24, 0.0, 120.0)
        elif v=='eta0'     : h = r.TH1F(histoName(sam, sys, sel, v), ';#eta_{l0}; entries/bin',                 26,-2.6,  +2.6)
        elif v=='eta1'     : h = r.TH1F(histoName(sam, sys, sel, v), ';#eta_{l1}; entries/bin',                 26,-2.6,  +2.6)
        elif v=='phi0'     : h = r.TH1F(histoName(sam, sys, sel, v), ';#phi_{l0} [rad]; entries/bin',           10, 0.0, twopi)
        elif v=='phi1'     : h = r.TH1F(histoName(sam, sys, sel, v), ';#phi_{l1} [rad]; entries/bin',           10, 0.0, twopi)
        elif v=='mcoll'    : h = r.TH1F(histoName(sam, sys, sel, v), ';m_{coll,l0,l1} [GeV]; entries/bin',      40, 0.0, 400.0)
        elif v=='mll'      : h = r.TH1F(histoName(sam, sys, sel, v), ';m_{l0,l1} [GeV]; entries/bin',           24, 0.0, 240.0)
        elif v=='ptll'     : h = r.TH1F(histoName(sam, sys, sel, v), ';p_{T,l0+l1} [GeV]; entries/bin',         24, 0.0, 240.0)
        elif v=='met'      : h = r.TH1F(histoName(sam, sys, sel, v), ';MET [GeV]; entries/bin',                 24, 0.0, 240.0)
        elif v=='dphil0met': h = r.TH1F(histoName(sam, sys, sel, v), ';#Delta#phi(l0, met) [rad]; entries/bin',  10, 0.0, twopi)
        elif v=='dphil1met': h = r.TH1F(histoName(sam, sys, sel, v), ';#Delta#phi(l1, met) [rad]; entries/bin',  10, 0.0, twopi)
        elif v=='nsj'      : h = r.TH1F(histoName(sam, sys, sel, v), ';N_{jets,20<pt<30};entries/bin',           10,-0.5,  9.5)
        elif v=='drl0csj'  : h = r.TH1F(histoName(sam, sys, sel, v), ';#DeltaR(l0, j_{close,soft});entries/bin',10, 0.0,   2.0)
        elif v=='drl1csj'  : h = r.TH1F(histoName(sam, sys, sel, v), ';#DeltaR(l1, j_{close,soft});entries/bin',10, 0.0,   2.0)
        elif v=='mcoll_vs_pt1'     : h = r.TH2F(histoName(sam, sys, sel, v), '; p_{T,l1} [GeV]; m_{coll,l0,l1} [GeV]',      48, 0.0, 240.0, 40, 0.0, 400.0)
        elif v=='pt0_vs_pt1'       : h = r.TH2F(histoName(sam, sys, sel, v), '; p_{T,l1} [GeV]; p_{T,l0} [GeV] [GeV]',      48, 0.0, 240.0, 48, 0.0, 240.0)
        elif v=='met_vs_pt1'       : h = r.TH2F(histoName(sam, sys, sel, v), '; p_{T,l1} [GeV]; MET [GeV]',                 48, 0.0, 240.0, 24, 0.0, 240.0)
        elif v=='dphil0met_vs_pt1' : h = r.TH2F(histoName(sam, sys, sel, v), '; p_{T,l1} [GeV]; #Delta#phi(l0, met) [rad]', 48, 0.0, 240.0, 10, 0.0, twopi)
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
                 ('fake',r.kGray), ('qflip', r.kYellow-9)]
    newColors = [] #[('signal',r.kMagenta), ('WW',r.kAzure-9), ('Higgs',r.kYellow-9)]
    colors = dict((g,c) for g,c in  oldColors + newColors)
    return colors[g]

def regions_to_plot():
    # return ['vr_emu_mue_ss'] # test to debug fake
    # return ['vr_emu_ss_razor']
    # return [k for k in selection_formulas().keys() if ('vr' in k and 'ss' in k)] # test to debug fake
    # return [k for k in selection_formulas().keys() if 'vr' not in k] # tmp until I have vrs
    # return [k for k in selection_formulas().keys() if 'sr' in k] # tmp dbg
    return ['sr_emu_os', 'sr_mue_os', 'vr_emu_os', 'vr_mue_os',
            'sr_emu_ss', 'sr_mue_ss', 'vr_emu_ss', 'vr_mue_ss',
            'ext_mumu_ss', 'ext_emu_mue_ss', 'ext_emu_pt0_40_ss', 'ext_mue_pt0_40_ss', 'ext_mumu_pt0_40_ss',
            'sr_mue_os_low_pt1_15'
            ]
    return selection_formulas().keys()

def variables_to_plot():
    return ['onebin', 'njets', 'pt0', 'pt1', 'd_pt0_pt1', 'eta0', 'eta1', 'phi0', 'phi1', 'mcoll',
            'mll', 'ptll', 'met', 'dphil0met', 'dphil1met',
            'drl0csj', 'drl1csj',
            'nsj',
            ]
def variables_to_fill():
    "do not plot 2d variables, but still fill the corresponding histograms"
    return variables_to_plot() + ['mcoll_vs_pt1', 'pt0_vs_pt1', 'met_vs_pt1', 'dphil0met_vs_pt1']

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
    if verbose : print "plotting ",padMaster.GetName(),' (',padMaster.GetEntries(),' entries)'
    can = r.TCanvas(canvasName, padMaster.GetTitle(), 800, 800)
    botPad, topPad = buildBotTopPads(can, squeezeMargins=False)
    can.cd()
    topPad.Draw()
    # draw top
    topPad.cd()
    topPad._hists = [padMaster]
    padMaster.Draw('axis')
    topPad.Update() # necessary to fool root's dumb object ownership of the stack
    stack = r.THStack('stack_'+padMaster.GetName(), '')
    r.SetOwnership(stack, False)
    topPad._hists.append(stack)
    leg = topRightLegend(can, 0.225, 0.325)
    topPad._leg = leg
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
        topPad._hists.append(histo)
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
    topPad.Update()
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
        topPad.SetTopMargin(2.0*topPad.GetTopMargin())
    atlasLabel = drawAtlasLabel(can, xpos=0.125, align=13, scale=0.75)
    if topLabel : topRightLabel(can, topLabel, ypos=1.0)
    yMin, yMax = getMinMax([histoData, dataGraph, histoTotBkg, histoSignal, totErrBand])
    padMaster.SetMinimum(0.0)
    padMaster.SetMaximum(1.1 * yMax)
    padMaster.GetXaxis().SetLabelSize(0)
    padMaster.GetXaxis().SetTitleSize(0)
    increaseAxisFont(padMaster.GetXaxis())
    increaseAxisFont(padMaster.GetYaxis())
    topPad.RedrawAxis()
    topPad.Update() # force stack to create padMaster
    # draw bot (ratio)
    can.cd()
    botPad.Draw()
    botPad.cd()
    ratio = buildRatioHistogram(histoData, histoTotBkg)
    yMin, yMax = 0.0, 2.0
    ratio.SetMinimum(yMin)
    ratio.SetMaximum(yMax)
    ratio.SetStats(0)
    ratio.Draw('axis')
    x_lo, x_hi = getXrange(ratio)
    refLines = [referenceLine(x_lo, x_hi, y, y) for y in [0.5, 1.0, 1.5]]
    for l in refLines : l.Draw()
    err_band_r = systUtils.buildErrBandRatioGraph(totErrBand)
    err_band_r.Draw('E2 same')
    ratio.Draw('ep same')
    xA, yA = ratio.GetXaxis(), ratio.GetYaxis()
    textScaleUp = 0.75*1.0/botPad.GetHNDC()
    # if xaxis_title : xA.SetTitle(xaxis_title)
    yA.SetNdivisions(-104)
    yA.SetTitle('Data/SM')
    yA.CenterTitle()
    yA.SetTitleOffset(yA.GetTitleOffset()/textScaleUp)
    xA.SetTitleSize(yA.GetTitleSize()) # was set to 0 for padmaster, restore it
    xA.SetLabelSize(yA.GetLabelSize())
    for a in [xA, yA] :
        a.SetLabelSize(a.GetLabelSize()*textScaleUp)
        a.SetTitleSize(a.GetTitleSize()*textScaleUp)
    botPad._graphical_objects = [ratio, err_band_r] + refLines # avoid garbage collection
    botPad.Update()
    can.Update()
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
