#!/bin/env python

# script to plot (from the fake ntuples) the variables (iso, pv) used to define the tight leptons
#
# davide.gerbaudo@gmail.com
# May 2014

import array
import collections
import glob
import kin
import math
import optparse
import os
import pprint
import time
from utils import (dictSum
                   ,first
                   ,mkdirIfNeeded
                   ,rmIfExists
                   )
import rootUtils
from rootUtils import (drawLegendWithDictKeys
                       ,buildErrBandGraph
                       ,computeStatErr2
                       ,getBinContents
                       ,getBinErrors
                       ,getMinMax
                       ,importRoot
                       ,importRootCorePackages
                       ,summedHisto
                       ,topRightLabel
                       ,rightLegend
                       ,writeObjectsToFile)
r = rootUtils.importRoot()
r.gROOT.SetStyle('Plain')
r.gStyle.SetPadTickX(1)
r.gStyle.SetPadTickY(1)
rootUtils.importRootCorePackages()
import dataset
import fakeUtils as fakeu
from indexed_chain import IndexedChain
import settings
import utils
from kin import addTlv, computeCollinearMassLepTau, computeRazor, computeMt, selection_formulas

usage="""
Example usage:
%prog \\
 --verbose  \\
 --tag ${TAG} \\
 --output-dir ./out/fakerate/el_sf_${TAG}
 >& log/fakerate/el_sf_${TAG}.log

 TODO
"""


'use input out/selection_tuple/Jan_22'

def main():
    parser = optparse.OptionParser(usage=usage)
    parser.add_option('-g', '--group', help='group to be processed (used only in fill mode)')
    parser.add_option('-i', '--input-dir', default='./out/fakerate')
    parser.add_option('-o', '--output-dir', default='./out/tight_variables_plots', help='dir for plots')
    parser.add_option('-l', '--lepton', default='el', help='either el or mu')
    parser.add_option('--samples-dir', default='samples/', help='directory with the list of samples; default ./samples/')
    parser.add_option('-f', '--fill-histos', action='store_true', default=False, help='force fill (default only if needed)')
    parser.add_option('--regions', default=None, help='comma-separated list of regions to consider')
    parser.add_option('--include-regions', default='.*', help='regexp to filter regions')
    parser.add_option('--exclude-regions', default=None, help='regext to exclude regions')
    # reminder: submit_batch_fill_job_per_group expects argument-less opt to default to False
    parser.add_option('--debug', action='store_true')
    parser.add_option('--verbose', action='store_true')
    parser.add_option('-b', '--batch',  action='store_true', help='submit to batch (used in fill mode)')
    parser.add_option('--list-all-regions', action='store_true', help='list all possible regions')
    parser.add_option('--require-tight-tight', action='store_true', help='fill histos only when both leps are tight')
    parser.add_option('--quick-test', action='store_true', help='run a quick test and fill only 1% of the events')
    parser.add_option('--disable-cache', action='store_true', help='disable the entry cache')

    (opts, args) = parser.parse_args()
    inputDir  = opts.input_dir
    lepton    = opts.lepton
    regions   = opts.regions
    verbose   = opts.verbose

    if opts.list_all_regions:
        print "All regions:\n\t%s"%'\n\t'.join(sorted(selection_formulas().keys()))
        return

    if lepton not in ['el', 'mu'] : parser.error("invalid lepton '%s'"%lepton)
    if opts.verbose : utils.print_running_conditions(parser, opts)

    runFill(opts)
    runPlot(opts)
    print 'todo: fix errorband'

#___________________________________________________________

def regions_to_plot(include='.*', exclude=None, regions=None):
    selected_regions = selection_formulas().keys()
    if regions:
        selected_regions = [r for r in selected_regions if r in regions.split(',')]
    selected_regions = utils.filterWithRegexp(selected_regions, include)
    selected_regions = utils.excludeWithRegexp(selected_regions, exclude) if exclude else selected_regions
    return selected_regions

def variables_to_plot():
    return ['onebin', 'njets', 'nsj', 'pt0', 'pt1', 'eta0', 'eta1',
            ]
def variables_to_fill():
    return variables_to_plot()

def runFill(opts):
    lepton    = opts.lepton
    batchMode = opts.batch
    inputDir  = opts.input_dir
    outputDir = opts.output_dir
    verbose   = opts.verbose
    debug     = opts.debug

    dataset.Dataset.verbose_parsing = True if debug else False
    groups = dataset.DatasetGroup.build_groups_from_files_in_dir(opts.samples_dir)
    if opts.group : groups = [g for g in groups if g.name==opts.group]
    if verbose : print '\n'.join("group {0} : {1} samples".format(g.name, len(g.datasets)) for g in groups)
    if debug :
        print '\n'.join("group {0} : {1} samples: {2}".format(g.name,
                                                              len(g.datasets),
                                                              '\n\t'+'\n\t'.join(d.name for d in g.datasets))
                        for g in groups)
    if verbose : print "filling histos"
    outputDir = outputDir+'/'+lepton+'/histos'
    mkdirIfNeeded(outputDir)
    print
    if batchMode:
        for group in groups:
            submit_batch_fill_job_per_group(group, opts)
    else:
        for group in groups:
            tree_name = 'hlfv_tuple'
            chain = IndexedChain(tree_name)
            for ds in group.datasets:
                chain.Add(os.path.join(inputDir, ds.name+'.root'))
            if opts.verbose:
                print "{0} : {1} entries from {2} samples".format(group.name,
                                                                  chain.GetEntries(),
                                                                  len(group.datasets))
            chain.cache_directory = os.path.abspath('./selection_cache/'+group.name+'/')
            tcuts = [r.TCut(reg, selection_formulas()[reg])
                     for reg in regions_to_plot(opts.include_regions, opts.exclude_regions, opts.regions)]
            chain.retrieve_entrylists(tcuts)
            counters_pre, histos_pre = dict(), dict()
            counters_npre, histos_npre = dict(), dict()
            cached_tcuts = [] if opts.disable_cache else chain.tcuts_with_existing_list()
            uncached_tcuts = tcuts if opts.disable_cache else chain.tcuts_without_existing_list()
            print 'todo: skip cuts for which the histo files are there'
            if verbose : print 'filling cached cuts: ',' '.join([c.GetName() for c in cached_tcuts])
            for cut in cached_tcuts:
                chain.preselect(cut)
                c_pre, h_pre = count_and_fill(chain=chain, opts=opts,
                                              group=group,
                                              cached_cut=cut)
                counters_pre = dictSum(counters_pre, c_pre)
                histos_pre = dictSum(histos_pre, h_pre)
            if verbose : print 'filling uncached cuts: ',' '.join([c.GetName() for c in uncached_tcuts])
            if uncached_tcuts:
                counters_npre, histos_npre = count_and_fill(chain=chain, opts=opts,
                                                            group=group,
                                                            noncached_cuts=uncached_tcuts)
                chain.save_lists()
            all_histos = dictSum(histos_pre, histos_npre)
            for sel, histos in all_histos.iteritems():
                # write histos for each sel to a separate file (finer granularity, better caching)
                out_filename = os.path.join(outputDir, group.name+'_'+sel+'.root')
                print out_filename
                writeObjectsToFile(out_filename, histos, verbose)

def runPlot(opts):
    lepton    = opts.lepton
    batchMode = opts.batch
    inputDir  = opts.input_dir
    outputDir = opts.output_dir
    verbose   = opts.verbose
    debug     = opts.debug
    dataset.Dataset.verbose_parsing = True if debug else False
    groups = dataset.DatasetGroup.build_groups_from_files_in_dir(opts.samples_dir)
    regions = regions_to_plot(opts.include_regions, opts.exclude_regions, opts.regions)

    inputDir  = outputDir+'/'+lepton+'/histos'
    outputDir = outputDir+'/'+lepton+'/plots'
    mkdirIfNeeded(outputDir)
    histonames = dict((g.name, histonamesOneSample(g.name, variables_to_plot(), regions, leptonSources))
                      for g in groups)
    groups_to_stack = [g.name for g in groups if not g.is_data]
    if verbose:
        print 'groups being included in the compositions: ',groups_to_stack
    for region in regions:
        all_histos = dict([(g.name,
                            rootUtils.fetchObjectsFromFile(os.path.join(inputDir, g.name+'_'+region+'.root'),
                                                           histonames[g.name][region],
                                                           verbose))
                           for g in groups])
        for v in variables_to_plot():
            histos = dict()
            for s in leptonSources:
                histos[s] = summedHisto(histos=[all_histos[g][v][s] for g in groups_to_stack],
                                        label='')
            plotStackedHistos(histos=histos, datakey='data', stackkeys=leptonSources,
                              outputDir=outputDir+'/'+region, region=region,
                              colors=fakeu.colorsFillSources(), verbose=verbose)
    return

#___________________________________________________

allLeptonSources = fakeu.allLeptonSources()
leptonSources = fakeu.leptonSources()
colorsFillSources = fakeu.colorsFillSources()
colorsLineSources = fakeu.colorsLineSources()
markersSources = fakeu.markersSources()
enum2source = fakeu.enum2source

def histoName(sample, var, selection, source):
    return 'h_'+sample+'_'+var+'_'+selection+'_'+source

def histonamesOneSample(sample_name, variables, selections, sources):
    "dict of histogram names with keys [sel][var][source]"
    hn = histoName
    return dict([(se,
                  dict([(v,
                         dict([(so, hn(sample_name, v, se, so))
                               for so in sources]))
                        for v in variables]))
                 for se in selections])

def book_histograms(sample_name, variables, selections, sources) :
    "dict of histograms with keys [sel][var][source]"
    h = settings.histogram
    hn = histoName
    return dict([(se,
                  dict([(v,
                         dict([(so, h(v, hn(sample_name, v, se, so)))
                               for so in sources]))
                        for v in variables]))
                 for se in selections])

def count_and_fill(chain,  opts, group=dataset.DatasetGroup('foo'),
                   cached_cut=None, noncached_cuts=[]):
    """
    count and fill for one sample (or group)
    """
    is_data = group.is_data
    is_mc = not is_data
    is_qflip_sample = group.is_qflip
    assert bool(cached_cut) != bool(noncached_cuts),"must choose either cached selection or non-cached selections: {}, {}".format(cached_cut, noncached_cuts)
    cuts = [cached_cut] if cached_cut else noncached_cuts
    if noncached_cuts:
        chain.preselect(None)
    selections = [c.GetName() for c in cuts]
    counters = dict((sel, 0) for sel in selections)
    quicktest = opts.quick_test
    verbose = opts.verbose

    histos = book_histograms(sample_name=group.name, variables=variables_to_fill(),
                             selections=selections, sources=leptonSources)
    weight_expr = 'event.pars.weight'
    start_time = time.clock()
    num_total_entries = chain.GetEntries()
    num_processed_entries = 0
    print 'todo: select right lepton'
    for iEntry, event in enumerate(chain):
        if quicktest and 100*iEntry > num_total_entries: break
        run_num = event.pars.runNumber
        evt_num = event.pars.eventNumber
        l0 = addTlv(event.l0)
        l1 = addTlv(event.l1)
        met = addTlv(event.met)
        l0_is_el, l0_is_mu = l0.isEl, l0.isMu
        l1_is_el, l1_is_mu = l1.isEl, l1.isMu
        is_emu = int(l0_is_el and l1_is_mu)
        is_mue = int(l0_is_mu and l1_is_el)
        is_mumu = int(l0_is_mu and l1_is_mu)
        is_ee = int(l0_is_el and l1_is_el)
        is_same_sign = int((l0.charge * l1.charge)>0)
        is_opp_sign  = not is_same_sign
        is_qflippable = is_opp_sign and (l0_is_el or l1_is_el) and is_mc
        weight = eval(weight_expr)
        l0_pt, l1_pt = l0.p4.Pt(), l1.p4.Pt()
        d_pt0_pt1 = l0_pt - l1_pt
        l0_eta, l1_eta = abs(l0.p4.Eta()), abs(l1.p4.Eta())
        dphi_l0_met = abs(l0.p4.DeltaPhi(met.p4))
        dphi_l1_met = abs(l1.p4.DeltaPhi(met.p4))
        dphi_l0_l1 = abs(l0.p4.DeltaPhi(l1.p4))
        dpt_l0_l1 = l0.p4.Pt()-l1.p4.Pt()
        # l0_source = None if is_data else enum2source(l0)
        # l1_source = None if is_data else enum2source(l1)
        l0_source = enum2source(l0)
        l1_source = enum2source(l1)
        def jet_pt2(j) : return j.px*j.px+j.py*j.py
        n_cl_jets = sum(1 for j in event.jets if jet_pt2(j)>30.*30.)
        n_jets = n_cl_jets + event.pars.numFjets + event.pars.numBjets
        pass_sels = {}
        for cut in cuts:
            sel = cut.GetName()
            sel_expr = cut.GetTitle()
            pass_sel = eval(sel_expr)
            pass_sels[sel] = pass_sel
            if not pass_sel : continue
            if pass_sel and not cached_cut : chain.add_entry_to_list(cut, iEntry)
            fill_weight = (weight)
            histos[sel]['onebin'][l0_source].Fill(1.0, fill_weight)
            histos[sel]['pt0'   ][l0_source].Fill(l0_pt, fill_weight)
            histos[sel]['pt1'   ][l1_source].Fill(l1_pt, fill_weight)
            histos[sel]['eta0'  ][l0_source].Fill(l0_eta, fill_weight)
            histos[sel]['eta1'  ][l1_source].Fill(l1_eta, fill_weight)
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
        for sel in selections:
            tot_integral = sum(h.Integral() for source, h in histos[sel]['onebin'].iteritems())
            tot_entries = sum(h.GetEntries() for source, h in histos[sel]['onebin'].iteritems())
            print "{0}: integral {1}, entries {2}".format(sel+'_onebin', tot_integral, tot_entries)
    return counters, histos

def plotStackedHistos(histos={}, datakey=None, stackkeys=[], outputDir='', region='', colors={}, verbose=False):
    "input: a dictionary of histos[group]"
    mkdirIfNeeded(outputDir)
    bkg_histos = dict([(k,h) for k,h in histos.iteritems() if k in stackkeys])
    tot_bkg = summedHisto(bkg_histos.values(), label='')
    err_band = None # tmp disable
    # err_band = buildErrBandGraph(tot_bkg, computeStatErr2(tot_bkg))
    empty_bkg = tot_bkg.Integral()==0
    if empty_bkg:
        if verbose : print "empty backgrounds, skip %s"%tot_bkg.GetName()
        return
    histoname = tot_bkg.GetName()
    can = r.TCanvas('c_'+histoname, histoname, 800, 600)
    can.cd()
    pm = tot_bkg # pad master
    pm.SetStats(False)
    pm.Draw('axis')
    can.Update() # necessary to fool root's dumb object ownership
    stack = r.THStack('stack_'+tot_bkg.GetName(),'')
    can.Update()
    r.SetOwnership(stack, False)
    for s, h in bkg_histos.iteritems() :
            h.SetFillColor(colors[s] if s in colors else r.kOrange)
            h.SetDrawOption('bar')
            h.SetDirectory(0)
            stack.Add(h)
    stack.Draw('hist same')
    # err_band.Draw('E2 same')
    data = histos[datakey] if datakey and datakey in histos else None
    if data and data.GetEntries():
        data.SetMarkerStyle(r.kFullDotLarge)
        data.Draw('p same')
    yMin, yMax = getMinMax([h for h in [tot_bkg, data, err_band] if h])
    pm.SetMinimum(0.0)
    pm.SetMaximum(1.1*yMax)
    can.Update()
    topRightLabel(can, "#splitline{%s}{%s}"%(histoname, region), xpos=0.125, align=13)
    drawLegendWithDictKeys(can, dictSum(bkg_histos, {'stat err':err_band}), opt='f')
    can.RedrawAxis()
    can._stack = stack
    can._histos = [h for h in stack.GetHists()]+[data]
    can.Update()
    if verbose : print os.path.join(outputDir, histoname+'.png')
    can.SaveAs(os.path.join(outputDir, histoname+'.png'))


if __name__=='__main__':
    main()
