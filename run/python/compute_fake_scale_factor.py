#!/bin/env python

# script to compute the electron (and muon) data/MC scale factor for the fake rate

# davide.gerbaudo@gmail.com
# April 2014

import array
import collections
import glob
import math
import optparse
import os
import pprint
import time
from utils import (dictSum
                   ,first
                   ,mkdirIfNeeded
                   )
import rootUtils
from rootUtils import (drawLegendWithDictKeys
                       ,getBinContents
                       ,getBinErrors
                       ,getMinMax
                       ,importRoot
                       ,importRootCorePackages
                       ,summedHisto
                       ,topRightLabel)
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
import kin
import systUtils
from kin import addTlv, computeCollinearMassLepTau, computeRazor, computeMt, selection_formulas
from plot_emu import getGroupColor

usage="""
Example usage:
%prog \\
 --verbose  \\
 --lepton el \\
 --output-dir ./out/fake_scale_factor/
 >& out/fake_scale_factor/
"""
def main():
    parser = optparse.OptionParser(usage=usage)
    parser.add_option('-g', '--group', help='group to be processed (used only in fill mode)')
    parser.add_option('-i', '--input-dir', default='./out/fakerate')
    parser.add_option('-o', '--output-dir', default='./out/fake_scale_factor')
    parser.add_option('-l', '--lepton', default='el', help='either el or mu')
    parser.add_option('-r', '--region', help='one of the regions for which we saved the fake ntuples')
    parser.add_option('--samples-dir', default='samples/', help='directory with the list of samples; default ./samples/')
    parser.add_option('-T', '--tight-def', help='on-the-fly tight def, one of defs in fakeUtils.py: fakeu.lepIsTight_std, etc.')
    parser.add_option('-f', '--fill-histos', action='store_true', default=False, help='force fill (default only if needed)')
    parser.add_option('--keep-real', action='store_true', default=False, help='do not subtract real (to get real lep efficiency)')
    parser.add_option('--debug', action='store_true')
    parser.add_option('--verbose', action='store_true')
    parser.add_option('--disable-cache', action='store_true', help='disable the entry cache')
    (options, args) = parser.parse_args()
    inputDir  = options.input_dir
    outputDir = options.output_dir
    lepton    = options.lepton
    region    = options.region
    keepreal  = options.keep_real
    debug     = options.debug
    verbose   = options.verbose
    if lepton not in ['el', 'mu'] : parser.error("invalid lepton '%s'"%lepton)
    regions = kin.selection_formulas().keys()
    assert region in regions,"invalid region '%s', must be one of %s"%(region, str(sorted(regions)))
    regions = [region]

    dataset.Dataset.verbose_parsing = True if debug else False
    groups = dataset.DatasetGroup.build_groups_from_files_in_dir(options.samples_dir)
    if options.group : groups = [g for g in groups if g.name==options.group]
    group_names = [g.name for g in groups]

    outputDir = outputDir+'/'+region+'/'+lepton # split the output in subdirectories, so we don't overwrite things
    mkdirIfNeeded(outputDir)
    templateOutputFilename = "scale_factor_{0}.root".format(lepton)
    outputFileName = os.path.join(outputDir, templateOutputFilename)
    cacheFileName = outputFileName.replace('.root', '_cache.root')
    doFillHistograms = options.fill_histos or not os.path.exists(cacheFileName)
    onthefly_tight_def = eval(options.tight_def) if options.tight_def else None # eval will take care of aborting on typos
    if verbose : utils.print_running_conditions(parser, options)
    vars = ['mt0', 'mt1', 'pt0', 'pt1', 'eta1', 'pt1_eta1']
    #fill histos
    if doFillHistograms :
        start_time = time.clock()
        num_processed_entries = 0
        histosPerGroup = bookHistos(vars, group_names, region=region)
        histosPerSource = bookHistosPerSource(vars, leptonSources, region=region)
        histosPerGroupPerSource = bookHistosPerSamplePerSource(vars, group_names, leptonSources, region=region)
        for group in groups:
            tree_name = 'ss3l_tuple'
            chain = IndexedChain(tree_name)
            for ds in group.datasets:
                fname = os.path.join(inputDir, ds.name+'.root')
                if os.path.exists(fname):
                    chain.Add(fname)
            if verbose:
                print "{0} : {1} entries from {2} samples".format(group.name, chain.GetEntries(), len(group.datasets))
            chain.cache_directory = os.path.abspath('./selection_cache/'+group.name+'/')
            tcuts = [r.TCut(reg, selection_formulas()[reg]) for reg in regions]
            print 'tcuts ',[c.GetName() for c in tcuts]
            chain.retrieve_entrylists(tcuts)
            counters_pre, histos_pre = dict(), dict()
            counters_npre, histos_npre = dict(), dict()
            print 'tcuts_with_existing_list ',str([c.GetName() for c in chain.tcuts_with_existing_list()])
            print 'tcuts_without_existing_list ',str([c.GetName() for c in chain.tcuts_without_existing_list()])
            cached_tcuts = [] if options.disable_cache else chain.tcuts_with_existing_list()
            print 'cached_tcuts ',[c.GetName() for c in cached_tcuts]
            uncached_tcuts = tcuts if options.disable_cache else chain.tcuts_without_existing_list()
            print 'todo: skip cuts for which the histo files are there'
            if verbose:
                print " --- group : {0} ---".format(group.name)
                print '\n\t'.join(chain.filenames)
            if verbose : print 'filling cached cuts: ',' '.join([c.GetName() for c in cached_tcuts])
            if verbose: print "%s : %d entries"%(group.name, chain.GetEntries())
            histosThisGroup = histosPerGroup[group.name]
            histosThisGroupPerSource = dict((v, histosPerGroupPerSource[v][group.name]) for v in histosPerGroupPerSource.keys())
            for cut in cached_tcuts:
                print 'cached_tcut ',cut
                chain.preselect(cut)
                num_processed_entries += fillHistos(chain, histosThisGroup, histosPerSource,
                                                    histosThisGroupPerSource,
                                                    lepton, group,
                                                    cut, cut_is_cached=True,
                                                    onthefly_tight_def=onthefly_tight_def,
                                                    verbose=verbose)
            if verbose : print 'filling uncached cuts: ',' '.join([c.GetName() for c in uncached_tcuts])
            if uncached_tcuts:
                assert len(uncached_tcuts)==1, "expecting only one cut, got {}".format(len(uncached_tcuts))
                cut = uncached_tcuts[0]
                chain.preselect(None)
                num_processed_entries += fillHistos(chain, histosThisGroup, histosPerSource,
                                                    histosThisGroupPerSource,
                                                    lepton, group,
                                                    cut, cut_is_cached=False,
                                                    onthefly_tight_def=onthefly_tight_def,
                                                    verbose=verbose)
                chain.save_lists()

        writeHistos(cacheFileName, histosPerGroup, histosPerSource, histosPerGroupPerSource, verbose)
        end_time = time.clock()
        delta_time = end_time - start_time
        if verbose:
            print ("processed {0:d} entries ".format(num_processed_entries)
                   +"in "+("{0:d} min ".format(int(delta_time/60)) if delta_time>60 else
                           "{0:.1f} s ".format(delta_time))
                   +"({0:.1f} kHz)".format(num_processed_entries/delta_time))
    # return
    # compute scale factors
    histosPerGroup = fetchHistos(cacheFileName, histoNames(vars, group_names, region), verbose)
    histosPerSource = fetchHistos(cacheFileName, histoNamesPerSource(vars, leptonSources, region), verbose)
    histosPerSamplePerSource = fetchHistos(cacheFileName, histoNamesPerSamplePerSource(vars, group_names, leptonSources, region), verbose)
    plotStackedHistos(histosPerGroup, outputDir+'/by_group', region, verbose)
    plotStackedHistosSources(histosPerSource, outputDir+'/by_source', region, verbose)
    plotPerSourceEff(histosPerVar=histosPerSource, outputDir=outputDir+'/by_source', lepton=lepton, region=region, verbose=verbose)
    for g in group_names:
        hps = dict((v, histosPerSamplePerSource[v][g])for v in vars)
        plotPerSourceEff(histosPerVar=hps, outputDir=outputDir, lepton=lepton, region=region, sample=g, verbose=verbose)


    hn_sf_eta = histoname_sf_vs_eta           (lepton)
    hn_sf_pt  = histoname_sf_vs_pt            (lepton)
    hn_da_eta = histoname_data_fake_eff_vs_eta(lepton)
    hn_da_pt  = histoname_data_fake_eff_vs_pt (lepton)
    subtractReal = not keepreal
    objs_eta = subtractRealAndComputeScaleFactor(histosPerGroup, 'eta1', hn_sf_eta, hn_da_eta, outputDir, region, subtractReal, verbose)
    objs_pt  = subtractRealAndComputeScaleFactor(histosPerGroup, 'pt1',  hn_sf_pt,  hn_da_pt,  outputDir, region, subtractReal, verbose)
    objs_pt_eta  = subtractRealAndComputeScaleFactor(histosPerGroup, 'pt1_eta1',
                                                     histoname_sf_vs_pt_eta(lepton),
                                                     histoname_data_fake_eff_vs_pt_eta(lepton),
                                                     outputDir, region, subtractReal, verbose)
    rootUtils.writeObjectsToFile(outputFileName, dictSum(dictSum(objs_eta, objs_pt), objs_pt_eta), verbose)
    if verbose : print "saved scale factors to %s" % outputFileName

#___________________________________________________

leptonTypes = fakeu.leptonTypes()
allLeptonSources = fakeu.allLeptonSources()
leptonSources = fakeu.leptonSources()
colorsFillSources = fakeu.colorsFillSources()
colorsLineSources = fakeu.colorsLineSources()
markersSources = fakeu.markersSources()
enum2source = fakeu.enum2source

def histoname_sf_vs_eta(l) : return 'sf_'+l+'_vs_eta'
def histoname_sf_vs_pt (l) : return 'sf_'+l+'_vs_pt'
def histoname_sf_vs_pt_eta(l) : return 'sf_'+l+'_vs_pt_vs_eta'
def histoname_data_fake_eff_vs_eta(l) : return l+'_fake_rate_data_vs_eta'
def histoname_data_fake_eff_vs_pt (l) : return l+'_fake_rate_data_vs_pt'
def histoname_data_fake_eff_vs_pt_eta(l) : return l+'_fake_rate_data_vs_pt_vs_eta'

def fillHistos(chain, histosThisGroup, histosPerSource, histosThisGroupPerSource,
               lepton, group, cut, cut_is_cached, onthefly_tight_def=None, verbose=False):
    nLoose, nTight = 0, 0
    totWeightLoose, totWeightTight = 0.0, 0.0
    normFactor = 1.0 if group=='heavyflavor' else 1.0 # bb/cc hand-waving normalization factor, see notes 2014-04-17
    addTlv, computeMt = kin.addTlv, kin.computeMt
    region = cut.GetName()
    sel = cut.GetName()
    sel_expr = cut.GetTitle()
    print group
    isData = group.name=='data'
    isHflf = region=='hflf'
    print 'TODO: check isHflf'
    print 'TODO: check hasTrigmatch (already in ntuple creation?)'
    isConversion = region=='conv'
    if group=='heavyflavor':
        if lepton=='el' and not isConversion : normFactor = 1.6
        if lepton=='mu' :                      normFactor = 0.87
    num_processed_entries = 0
    sourceReal = r.ss3l.LeptonTruthType.Prompt
    for iEvent, event in enumerate(chain) :
        num_processed_entries += 1
        pars = event.pars
        weight, evtN, runN = pars.weight, pars.eventNumber, pars.runNumber
        hasTrigmatch = True # pars.has2ltrigmatch==1
        weight = weight*normFactor
        tag, probe, met = addTlv(event.l0), addTlv(event.l1), addTlv(event.met)
        if not (tag.isMu and probe.isMu) : continue
        isSameSign = tag.charge*probe.charge > 0.
        isRightProbe = probe.isEl if lepton=='el' else probe.isMu if lepton=='mu' else False
        if not isRightProbe and lepton=='el':
            tag, probe = probe, tag # if we are  picking emu instead of mue
            isRightProbe = probe.isEl if lepton=='el' else False

        tagIsTight = onthefly_tight_def(tag) if onthefly_tight_def else tag.isTight
        isTight = onthefly_tight_def(probe) if onthefly_tight_def else probe.isTight # todo: rename to probeIsTight
        probeSource = probe.source
        isReal = probeSource==sourceReal and not isData
        isFake = not isReal and not isData
        jets = event.jets
        jets = [addTlv(j) for j in jets] # only if needed
        def isBjet(j, mv1_80=0.3511) : return j.mv1 > mv1_80 # see SusyDefs.h
        hasBjets = any(isBjet(j) for j in jets) # compute only if necessary
        hasFjets = any(abs(j.p4.Eta())>2.4 and j.p4.Pt()>30. for j in jets)
        hasCLjets = any(abs(j.p4.Eta())<2.4 and j.p4.Pt()>30 and not isBjet(j) for j in jets)
        hasJ30jets = any(j.p4.Pt()>30. for j in jets)
        hasJets = hasBjets or hasFjets or hasCLjets
        tag4m, probe4m, met4m = r.TLorentzVector(), r.TLorentzVector(), r.TLorentzVector()
        tag4m.SetPxPyPzE(tag.px, tag.py, tag.pz, tag.E)
        probe4m.SetPxPyPzE(probe.px, probe.py, probe.pz, probe.E)
        met4m.SetPxPyPzE(met.px, met.py, met.pz, met.E)
        pt = probe4m.Pt()
        eta = abs(probe4m.Eta())
        mt0 = computeMt(tag4m, met4m)
        mt1 = computeMt(probe4m, met4m)
        pt0 = tag4m.Pt()
        pt1 = probe4m.Pt()
        passTrigBias =  True
        if   isHflf       : passTrigBias = pt0>20.0 and pt1>20.0
        elif isConversion : passTrigBias = pt1>20.0

        is_mumu = tag.isMu and probe.isMu
        is_same_sign = isSameSign
        is_opp_sign = not is_same_sign
        tag_is_mu = tag.isMu
        tag_is_tight = tagIsTight
        probe_is_tight = isTight
        probe_is_mu = probe.isMu
        probe_is_el = probe.isEl
        l0_is_el, l0_is_mu = event.l0.isEl, event.l0.isMu
        l1_is_el, l1_is_mu = event.l1.isEl, event.l1.isMu
        is_ee   = l0_is_el and l1_is_el
        is_emu  = l0_is_el and l1_is_mu
        is_mue  = l0_is_mu and l1_is_el
        is_mumu = l0_is_mu and l1_is_mu
        m_ll = (tag4m+probe4m).M()
        pass_sel = eval(sel_expr)
        if pass_sel and not cut_is_cached : chain.add_entry_to_list(cut, iEvent)
        # if tag.isMu and isRightProbe and isSameSign : # test 1 : no jet req
        # if tag.isMu and tag.isTight and isRightProbe and isSameSign : # test 2 : reproduce counts from ss3l
        def fillHistosBySource(probe):
            leptonSource = enum2source(probe)
            # if leptonSource=='Unknown':
            #     print 'unknown lep from ',group.name,' pt ',pt,' eta ',eta,' file ',chain.GetCurrentFile().GetName()
            def fillPerSource(tightOrLoose):
                histosPerSource         ['mt1' ][leptonSource][tightOrLoose].Fill(mt1, weight)
                histosPerSource         ['pt1' ][leptonSource][tightOrLoose].Fill(pt,  weight)
                histosPerSource         ['eta1'][leptonSource][tightOrLoose].Fill(eta, weight)
                histosThisGroupPerSource['mt1' ][leptonSource][tightOrLoose].Fill(mt1, weight)
                histosThisGroupPerSource['pt1' ][leptonSource][tightOrLoose].Fill(pt,  weight)
                histosThisGroupPerSource['eta1'][leptonSource][tightOrLoose].Fill(eta, weight)
                histosThisGroupPerSource['pt1_eta1'][leptonSource][tightOrLoose].Fill(pt, eta, weight)
            fillPerSource('loose')
            if isTight :
                fillPerSource('tight')
        def fill(lepType=''):
            histosThisGroup['mt0' ][lepType].Fill(mt0, weight)
            histosThisGroup['pt0' ][lepType].Fill(pt0, weight)
            histosThisGroup['mt1' ][lepType].Fill(mt1, weight)
            histosThisGroup['pt1' ][lepType].Fill(pt, weight)
            histosThisGroup['eta1'][lepType].Fill(eta, weight)
            histosThisGroup['pt1_eta1'][lepType].Fill(pt, eta, weight)

        sourceIsKnown = not isData
        nLoose         += 1
        totWeightLoose += weight
        if isTight:
            nTight         += 1
            totWeightTight += weight

        fill(lepType='loose')
        if sourceIsKnown:
            fillHistosBySource(probe)
        if isTight            : fill(lepType='tight')
        if isReal             : fill(lepType='real_loose')
        if isFake             : fill(lepType='fake_loose')
        if isReal and isTight : fill(lepType='real_tight')
        if isFake and isTight : fill(lepType='fake_tight')
    if verbose:
        counterNames = ['nLoose', 'nTight', 'totWeightLoose', 'totWeightTight']
        print ', '.join(["%s : %.1f"%(c, eval(c)) for c in counterNames])
    return num_processed_entries

def histoNamePerSample(var, sample, tightOrLoose, region) : return 'h_'+var+'_'+sample+'_'+tightOrLoose+'_'+region
def histoNamePerSource(var, leptonSource, tightOrLoose, region) : return 'h_'+var+'_'+leptonSource+'_'+tightOrLoose+'_'+region
def histoNamePerSamplePerSource(var, sample, leptonSource, tightOrLoose, region) : return 'h_'+var+'_'+sample+'_'+leptonSource+'_'+tightOrLoose+'_'+region

def bookHistos(variables, samples, leptonTypes=leptonTypes, region='') :
    "book a dict of histograms with keys [sample][var][tight, loose, real_tight, real_loose]"
    def histo(variable, hname):
        h = None
        v = variable
        mtBinEdges = fakeu.mtBinEdges()
        ptBinEdges = fakeu.ptBinEdges()
        etaBinEdges = fakeu.etaBinEdges()
        if   v=='mt0'     : h = r.TH1F(hname, ';m_{T}(tag,MET) [GeV]; entries/bin',   len(mtBinEdges)-1,  mtBinEdges)
        elif v=='mt1'     : h = r.TH1F(hname, ';m_{T}(probe,MET) [GeV]; entries/bin', len(mtBinEdges)-1,  mtBinEdges)
        elif v=='pt0'     : h = r.TH1F(hname, ';p_{T,l0} [GeV]; entries/bin',   len(ptBinEdges)-1,  ptBinEdges)
        elif v=='pt1'     : h = r.TH1F(hname, ';p_{T,l1} [GeV]; entries/bin',   len(ptBinEdges)-1,  ptBinEdges)
        elif v=='eta1'    : h = r.TH1F(hname, ';#eta_{l1}; entries/bin',        len(etaBinEdges)-1, etaBinEdges)
        elif v=='pt1_eta1': h = r.TH2F(hname, ';p_{T} [GeV]; #eta',
                                       len(ptBinEdges)-1,  ptBinEdges,
                                       len(etaBinEdges)-1, etaBinEdges)
        else : print "unknown variable %s"%v
        h.SetDirectory(0)
        h.Sumw2()
        return h
    return dict([(s,
                  dict([(v,
                         dict([(lt, histo(variable=v, hname=histoNamePerSample(v, s, lt, region)))
                               for lt in leptonTypes]))
                        for v in variables]))
                 for s in samples])

def bookHistosPerSource(variables, sources, region=''):
    "book a dict of histograms with keys [var][lepton_source][tight, loose]"
    def histo(variable, hname):
        h = None
        mtBinEdges = fakeu.mtBinEdges()
        ptBinEdges = fakeu.ptBinEdges()
        etaBinEdges = fakeu.etaBinEdges()
        if   variable=='mt0'     : h = r.TH1F(hname, ';m_{T}(tag,MET) [GeV]; entries/bin', len(mtBinEdges)-1,   mtBinEdges)
        elif variable=='mt1'     : h = r.TH1F(hname, ';m_{T}(probe,MET) [GeV]; entries/bin', len(mtBinEdges)-1, mtBinEdges)
        elif variable=='pt0'     : h = r.TH1F(hname, ';p_{T,l0} [GeV]; entries/bin',   len(ptBinEdges)-1,  ptBinEdges)
        elif variable=='pt1'     : h = r.TH1F(hname, ';p_{T,l1} [GeV]; entries/bin',   len(ptBinEdges)-1,  ptBinEdges)
        elif variable=='eta1'    : h = r.TH1F(hname, ';#eta_{l1}; entries/bin',        len(etaBinEdges)-1, etaBinEdges)
        elif v=='pt1_eta1': h = r.TH2F(hname, ';p_{T} [GeV]; #eta',
                                       len(ptBinEdges)-1,  ptBinEdges,
                                       len(etaBinEdges)-1, etaBinEdges)
        else : print "unknown variable %s"%v
        h.SetDirectory(0)
        h.Sumw2()
        return h
    return dict([(v,
                  dict([(s,
                         {'loose' : histo(variable=v, hname=histoNamePerSource(v, s, 'loose', region)),
                          'tight' : histo(variable=v, hname=histoNamePerSource(v, s, 'tight', region)),
                          })
                        for s in leptonSources]))
                 for v in variables])

def bookHistosPerSamplePerSource(variables, samples, sources, region=''):
    "book a dict of histograms with keys [var][sample][lepton_source][tight, loose]"
    def histo(variable, hname):
        h = None
        mtBinEdges = fakeu.mtBinEdges()
        ptBinEdges = fakeu.ptBinEdges()
        etaBinEdges = fakeu.etaBinEdges()
        if   variable=='mt0'     : h = r.TH1F(hname, ';m_{T}(tag,MET) [GeV]; entries/bin', len(mtBinEdges)-1,  mtBinEdges)
        elif variable=='mt1'     : h = r.TH1F(hname, ';m_{T}(probe,MET) [GeV]; entries/bin', len(mtBinEdges)-1,  mtBinEdges)
        elif variable=='pt0'     : h = r.TH1F(hname, ';p_{T,l0} [GeV]; entries/bin',   len(ptBinEdges)-1,  ptBinEdges)
        elif variable=='pt1'     : h = r.TH1F(hname, ';p_{T,l1} [GeV]; entries/bin',   len(ptBinEdges)-1,  ptBinEdges)
        elif variable=='eta1'    : h = r.TH1F(hname, ';#eta_{l1}; entries/bin',        len(etaBinEdges)-1, etaBinEdges)
        elif v=='pt1_eta1': h = r.TH2F(hname, ';p_{T} [GeV]; #eta',
                                       len(ptBinEdges)-1,  ptBinEdges,
                                       len(etaBinEdges)-1, etaBinEdges)
        else : print "unknown variable %s"%v
        h.SetDirectory(0)
        h.Sumw2()
        return h
    return dict([(v,
                  dict([(g,
                         dict([(s,
                                {'loose' : histo(variable=v, hname=histoNamePerSamplePerSource(v, g, s, 'loose', region)),
                                 'tight' : histo(variable=v, hname=histoNamePerSamplePerSource(v, g, s, 'tight', region)),
                                 })
                               for s in leptonSources]))
                        for g in samples]))
                 for v in variables])

def extractName(dictOrHist):
    "input must be either a dict or something with 'GetName'"
    isDict = type(dictOrHist) is dict
    return dict([(k, extractName(v)) for k,v in dictOrHist.iteritems()]) if isDict else dictOrHist.GetName()
def histoNames(variables, samples, region) :
    return extractName(bookHistos(variables, samples, region=region))
def histoNamesPerSource(variables, samples, region) :
    return extractName(bookHistosPerSource(variables, leptonSources, region=region))
def histoNamesPerSamplePerSource(variables, samples, leptonSources, region) :
    return extractName(bookHistosPerSamplePerSource(variables, samples, leptonSources, region=region))

def writeHistos(outputFileName='', histosPerGroup={}, histosPerSource={}, histosPerSamplePerGroup={}, verbose=False):
    rootUtils.writeObjectsToFile(outputFileName, [histosPerGroup, histosPerSource, histosPerSamplePerGroup], verbose)
def fetchHistos(fileName='', histoNames={}, verbose=False):
    return rootUtils.fetchObjectsFromFile(fileName, histoNames, verbose)

def plotStackedHistos(histosPerGroup={}, outputDir='', region='', verbose=False):
    groups = histosPerGroup.keys()
    variables = first(histosPerGroup).keys()
    leptonTypes = first(first(histosPerGroup)).keys()
    colors = getGroupColor()
    mkdirIfNeeded(outputDir)
    histosPerName = dict([(region+'_'+var+'_'+lt, # one canvas for each histo, so key with histoname w/out group
                           dict([(g, histosPerGroup[g][var][lt]) for g in groups]))
                          for var in variables for lt in leptonTypes])
    for histoname, histosPerGroup in histosPerName.iteritems():
        missingGroups = [g for g, h in histosPerGroup.iteritems() if not h]
        if missingGroups:
            if verbose : print "skip %s, missing histos for %s"%(histoname, str(missingGroups))
            continue
        bkgHistos = dict([(g, h) for g, h in histosPerGroup.iteritems() if g not in ['data', 'signal']])
        totBkg = summedHisto(bkgHistos.values())
        err_band = None # buildErrBandGraph(totBkg, computeStatErr2(totBkg))
        emptyBkg = totBkg.Integral()==0
        if emptyBkg:
            if verbose : print "empty backgrounds, skip %s"%histoname
            continue
        can = r.TCanvas('c_'+histoname, histoname, 800, 600)
        can.cd()
        pm = totBkg # pad master
        pm.SetStats(False)
        pm.Draw('axis')
        can.Update() # necessary to fool root's dumb object ownership
        stack = r.THStack('stack_'+histoname,'')
        can.Update()
        r.SetOwnership(stack, False)
        for s, h in bkgHistos.iteritems() :
            h.SetFillColor(colors[s] if s in colors else r.kOrange)
            h.SetDrawOption('bar')
            h.SetDirectory(0)
            stack.Add(h)
        stack.Draw('hist same')
        # err_band.Draw('E2 same')
        data = histosPerGroup['data']
        if data and data.GetEntries():
            data.SetMarkerStyle(r.kFullDotLarge)
            data.Draw('p same')
        # yMin, yMax = getMinMax([h for h in [totBkg, data, err_band] if h]) # fixme with err_band
        yMin, yMax = 0.0, data.GetMaximum()
        pm.SetMinimum(0.0)
        pm.SetMaximum(1.1*yMax)
        can.Update()
        topRightLabel(can, histoname, xpos=0.125, align=13)
        # drawLegendWithDictKeys(can, dictSum(bkgHistos, {'stat err':err_band}), opt='f')
        drawLegendWithDictKeys(can, bkgHistos, opt='f')
        can.RedrawAxis()
        can._stack = stack
        can._histos = [h for h in stack.GetHists()]+[data]
        can.Update()
        outFname = os.path.join(outputDir, histoname+'.png')
        utils.rmIfExists(outFname)
        can.SaveAs(outFname)
def plotStackedHistosSources(histosPerVar={}, outputDir='', region='', verbose=False):
    variables = histosPerVar.keys()
    sources = first(histosPerVar).keys()
    colors = colorsFillSources
    mkdirIfNeeded(outputDir)
    for var in variables:
        for lOrT in ['loose', 'tight']:
            histos = dict((s, histosPerVar[var][s][lOrT]) for s in sources)
            canvasBasename = region+'_region_'+var+'_'+lOrT
            missingSources = [s for s, h in histos.iteritems() if not h]
            if missingSources:
                if verbose : print "skip %s, missing histos for %s"%(var, str(missingSources))
                continue
            totBkg = summedHisto(histos.values())
            err_band = None # buildErrBandGraph(totBkg, computeStatErr2(totBkg))
            emptyBkg = totBkg.Integral()==0
            if emptyBkg:
                if verbose : print "empty backgrounds, skip %s"%canvasBasename
                continue
            can = r.TCanvas('c_'+canvasBasename, canvasBasename, 800, 600)
            can.cd()
            pm = totBkg # pad master
            pm.SetStats(False)
            pm.Draw('axis')
            can.Update() # necessary to fool root's dumb object ownership
            stack = r.THStack('stack_'+canvasBasename,'')
            can.Update()
            r.SetOwnership(stack, False)
            for s, h in histos.iteritems() :
                h.SetFillColor(colors[s] if s in colors else r.kOrange)
                h.SetDrawOption('bar')
                h.SetDirectory(0)
                stack.Add(h)
            stack.Draw('hist same')
            # err_band.Draw('E2 same')
            yMin, yMax = getMinMax([h for h in [totBkg, err_band] if h is not None])
            pm.SetMinimum(0.0)
            pm.SetMaximum(1.1*yMax)
            can.Update()
            topRightLabel(can, canvasBasename, xpos=0.125, align=13)
            # drawLegendWithDictKeys(can, dictSum(histos, {'stat err':err_band}), opt='f')
            drawLegendWithDictKeys(can, histos, opt='f')
            can.RedrawAxis()
            can._stack = stack
            can._histos = [h for h in stack.GetHists()]
            can.Update()
            outFname = os.path.join(outputDir, canvasBasename+'.png')
            utils.rmIfExists(outFname)
            can.SaveAs(outFname)

def plotPerSourceEff(histosPerVar={}, outputDir='', lepton='', region='', sample='', verbose=False, zoomIn=True):
    "plot efficiency for each source (and 'anysource') as a function of each var; expect histos[var][source][loose,tight]"
    variables = histosPerVar.keys()
    sources = [s for s in first(histosPerVar).keys() if s!='real'] # only fake eff really need a scale factor
    colors = colorsLineSources
    mkdirIfNeeded(outputDir)
    for var in filter(lambda x : x in ['pt1', 'eta1'], histosPerVar.keys()):
        histosPerSource = dict((s, histosPerVar[var][s]) for s in sources)
        canvasBasename = region+'_efficiency_'+lepton+'_'+var+("_%s"%sample if sample else '')
        missingSources = [s for s, h in histosPerSource.iteritems() if not h['loose'] or not h['tight']]
        if missingSources:
            if verbose : print "skip %s, missing histos for %s"%(var, str(missingSources))
            continue
        anySourceLoose = summedHisto([h['loose'] for h in histosPerSource.values()])
        anySourceTight = summedHisto([h['tight'] for h in histosPerSource.values()])
        anySourceLoose.SetName(histoNamePerSource(var, 'any', 'loose', region))
        anySourceTight.SetName(histoNamePerSource(var, 'any', 'tight', region))
        histosPerSource['any'] = { 'loose' : anySourceLoose, 'tight' : anySourceTight }
        emptyBkg = anySourceLoose.Integral()==0 or anySourceTight.Integral()==0
        if emptyBkg:
            if verbose : print "empty backgrounds, skip %s"%canvasBasename
            continue
        def computeEfficiencies(histosPerSource={}) :
            sources = histosPerSource.keys()
            num = dict((s, histosPerSource[s]['tight']) for s in sources)
            den = dict((s, histosPerSource[s]['loose']) for s in sources)
            eff = dict((s, h.Clone(h.GetName().replace('tight', 'tight_over_loose')))
                       for s, h in num.iteritems())
            [eff[s].Divide(den[s]) for s in sources]
            return eff

        effs = computeEfficiencies(histosPerSource)
        can = r.TCanvas('c_'+canvasBasename, canvasBasename, 800, 600)
        can.cd()
        pm = first(effs) # pad master
        pm.SetStats(False)
        pm.Draw('axis')
        can.Update()
        for s, h in effs.iteritems() :
            h.SetMarkerColor(colors[s] if s in colors else r.kBlack)
            h.SetLineColor(h.GetMarkerColor())
            h.SetLineWidth(2*h.GetLineWidth())
            h.SetMarkerStyle(markersSources[s] if s in markersSources else r.kDot)
            h.Draw('ep same')
            h.SetDirectory(0)
        #pprint.pprint(effs)
        yMin, yMax = getMinMax(effs.values())
        pm.SetMinimum(0.0)
        pm.SetMaximum(0.25 if yMax < 0.5 and zoomIn else 1.1)
        can.Update()
        topRightLabel(can, canvasBasename, xpos=0.125, align=13)
        drawLegendWithDictKeys(can, effs, opt='lp')
        can.RedrawAxis()
        can._histos = effs
        can.Update()
        outFname = os.path.join(outputDir, canvasBasename+'.png')
        utils.rmIfExists(outFname)
        can.SaveAs(outFname)

def subtractRealAndComputeScaleFactor(histosPerGroup={}, variable='', outRatiohistoname='',outDataeffhistoname='',
                                      outputDir='./', region='', subtractReal=True, verbose=False):
    "efficiency scale factor"
    groups = histosPerGroup.keys()
    mkdirIfNeeded(outputDir)
    histosPerType = dict([(lt,
                           dict([(g,
                                  histosPerGroup[g][variable][lt])
                                 for g in groups]))
                          for lt in leptonTypes])
    for lt in leptonTypes :
        histosPerType[lt]['totSimBkg'] = summedHisto([histo for group,histo in histosPerType[lt].iteritems()
                                                      if group not in ['data', 'signal']])

    simuTight = histosPerType['fake_tight']['totSimBkg']
    simuLoose = histosPerType['fake_loose']['totSimBkg']
    dataTight = histosPerType['tight'     ]['data'     ]
    dataLoose = histosPerType['loose'     ]['data'     ]
    # subtract real contribution from data
    # _Note to self_: currently estimating the real contr from MC; in
    # the past also used iterative corr, which might be more
    # appropriate in cases like here, where the normalization is
    # so-so.  Todo: investigate the normalization.
    dataSubTight = dataTight.Clone(dataTight.GetName().replace('data_tight','data_minus_prompt_tight'))
    dataSubLoose = dataLoose.Clone(dataLoose.GetName().replace('data_loose','data_minus_prompt_loose'))
    dataSubTight.SetDirectory(0)
    dataSubLoose.SetDirectory(0)
    dataSubTight.Add(histosPerType['real_tight']['totSimBkg'], -1.0 if subtractReal else 0.0)
    dataSubLoose.Add(histosPerType['real_loose']['totSimBkg'], -1.0 if subtractReal else 0.0)
    effData = dataSubTight.Clone(outDataeffhistoname)
    effData.SetDirectory(0)
    effData.Divide(dataSubLoose)
    effSimu = simuTight.Clone(simuTight.GetName().replace('fake_tight','fake_eff'))
    effSimu.SetDirectory(0)
    effSimu.Divide(simuLoose)
    print "eff(T|L) vs. ",variable
    def formatFloat(floats): return ["%.4f"%f for f in floats]
    print "efficiency data : ",formatFloat(getBinContents(effData))
    print "efficiency simu : ",formatFloat(getBinContents(effSimu))
    ratio = effData.Clone(outRatiohistoname)
    ratio.SetDirectory(0)
    ratio.Divide(effSimu)
    print "sf    data/simu : ",formatFloat(getBinContents(ratio))
    print "            +/- : ",formatFloat(getBinErrors(ratio))
    can = r.TCanvas('c_'+outRatiohistoname, outRatiohistoname, 800, 600)
    botPad, topPad = rootUtils.buildBotTopPads(can)
    can.cd()
    topPad.Draw()
    topPad.cd()
    pm = effData
    pm.SetStats(0)
    pm.Draw('axis')
    xAx, yAx = pm.GetXaxis(), pm.GetYaxis()
    xAx.SetTitle('')
    xAx.SetLabelSize(0)
    yAx.SetRangeUser(0.0, 0.25)
    textScaleUp = 1.0/topPad.GetHNDC()
    yAx.SetLabelSize(textScaleUp*0.04)
    yAx.SetTitleSize(textScaleUp*0.04)
    yAx.SetTitle('#epsilon(T|L)')
    yAx.SetTitleOffset(yAx.GetTitleOffset()/textScaleUp)
    effSimu.SetLineColor(r.kRed)
    effSimu.SetMarkerStyle(r.kOpenCross)
    effSimu.SetMarkerColor(effSimu.GetLineColor())
    effData.Draw('same')
    effSimu.Draw('same')
    leg = drawLegendWithDictKeys(topPad, {'data':effData, 'simulation':simuTight}, legWidth=0.4)
    leg.SetHeader('scale factor '+region+' '+('electron' if '_el_'in outRatiohistoname else
                                              'muon' if '_mu_' in outRatiohistoname else ''))
    can.cd()
    botPad.Draw()
    botPad.cd()
    ratio.SetStats(0)
    ratio.Draw()
    textScaleUp = 1.0/botPad.GetHNDC()
    xAx, yAx = ratio.GetXaxis(), ratio.GetYaxis()
    yAx.SetRangeUser(0.0, 2.0)
    xAx.SetTitle({'pt1':'p_{T}', 'eta1':'|#eta|', 'pt1_eta1':'p_{T}'}[variable])
    yAx.SetNdivisions(-202)
    yAx.SetTitle('Data/Sim')
    yAx.CenterTitle()
    xAx.SetLabelSize(textScaleUp*0.04)
    xAx.SetTitleSize(textScaleUp*0.04)
    yAx.SetLabelSize(textScaleUp*0.04)
    yAx.SetTitleSize(textScaleUp*0.04)
    refLine = rootUtils.referenceLine(xAx.GetXmin(), xAx.GetXmax())
    refLine.Draw()
    can.Update()
    outFname = os.path.join(outputDir, region+'_'+outRatiohistoname)
    for ext in ['.eps','.png']:
        utils.rmIfExists(outFname+ext)
        can.SaveAs(outFname+ext)
    return {outRatiohistoname : ratio,
            outDataeffhistoname : effData,
            outDataeffhistoname.replace('_fake_rate_data_', '_tight_data_minus_prompt') : dataSubTight,
            outDataeffhistoname.replace('_fake_rate_data_', '_loose_data_minus_prompt') : dataSubLoose
            }

if __name__=='__main__':
    main()
