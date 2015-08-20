#!/bin/env python

# plot the signal acceptance for our hlfv signal
#
# davide.gerbaudo@gmail.com
# July 2014

import glob
import os
import utils
ntutils = utils.import_susyntutils()
ntcutflow = utils.import_susyntcutflow()

r = ntutils.import_root()

def main():
    ntutils.load_packages()
    ntutils.generate_dicts()
    ntutils.import_SUSYDefs_enums()
    input_dir = '/gdata/atlas/gerbaudo/trashbin/hlfv_signal_shikma/'
    input_files = {
        'ggHtauEl' : "%s/HiggsTauMu_ggH.root"%input_dir,
        'ggHtauMu' : "%s/HiggsTauEl_ggH.root"%input_dir,
        'vbfHtauEl' : "%s/HiggsTauEl_VBF.root"%input_dir,
        'vbfHtauMu' : "%s/HiggsTauMu_VBF.root"%input_dir
        }
    histos = {}
    for label, input_file in input_files.iteritems():
        chain = r.TChain('susyNt')
        chain.Add(input_file)
        num_entries = chain.GetEntries()
        num_entries_to_process = num_entries if num_entries<1e4 else int(1e4)
        print "About to loop on %d entries"%num_entries_to_process
        histos[label] = run_with_chain(chain, num_entries_to_process)
    plot_histos(histos)

def run_with_chain(tree, n_max_entries=-1):
    verbose = False
    nttool = r.SusyNtTools()
    m_entry = r.Long(-1)
    ntevent = r.Susy.SusyNtObject(m_entry)
    ntevent.ReadFrom(tree)
    isSimplifiedModel = False
    nttool.buildSumwMap(tree, isSimplifiedModel)
    period, useRewUtils = 'Moriond', False
    trig_logic = r.DilTrigLogic(period, useRewUtils)
    mcweighter = r.MCWeighter(tree)
    mcweighter.parseAdditionalXsecFile('${ROOTCOREBIN}/data/susynt-ss2lb/LFV.txt', verbose)
    n_entries_to_print = 4
    sys = ntutils.SusyNtSys.NtSys_NOM
    tauId = ntutils.TauID
    tauJetId, tauEleId, tauMuoId = tauId.TauID_loose, tauId.TauID_medium, tauId.TauID_medium
    cutflow = ntcutflow.Cutflow()
    SkipEvent = ntcutflow.SkipEvent
    for iEntry, entry in enumerate(tree):
        m_entry = iEntry
        if n_max_entries>0 and m_entry >= n_max_entries : break
        if iEntry < n_entries_to_print : print 'run ', ntevent.evt().run,' event ',ntevent.evt().event
        pre_elecs  = nttool.getPreElectrons(ntevent, sys)
        pre_muons  = nttool.getPreMuons(ntevent, sys)
        pre_taus   = nttool.getPreTaus(ntevent, sys)
        pre_jets   = nttool.getPreJets(ntevent, sys)
        nttool.performOverlap(pre_elecs, pre_muons, pre_taus, pre_jets) # pre_ is baseline_ after removal
        nttool.removeSFOSPair(pre_elecs, 12.0)
        nttool.removeSFOSPair(pre_muons, 12.0)
        rmLepsFromIso = False
        n_vertices = ntevent.evt().nVtx
        is_mc = ntevent.evt().isMC
        sig_elecs = nttool.getSignalElectrons(pre_elecs, pre_muons, n_vertices, is_mc, rmLepsFromIso)
        sig_muons = nttool.getSignalMuons(pre_muons, pre_elecs, n_vertices, is_mc, rmLepsFromIso)
        sig_taus = nttool.getSignalTaus(pre_taus, tauJetId, tauEleId, tauMuoId)
        sig_jets = nttool.getSignalJets(pre_jets, sys)
        sig_jets2l = nttool.getSignalJets2Lep(pre_jets, sys)
        met = nttool.getMet(ntevent, sys)
        pre_lep, sig_lep = r.LeptonVector(), r.LeptonVector()
        nttool.buildLeptons(pre_lep, pre_elecs, pre_muons)
        nttool.buildLeptons(sig_lep, sig_elecs, sig_muons)
        if iEntry<n_entries_to_print:
            print 'pre_lep:\n','\n'.join(["[%d] %s (eta,phi,pt) = (%.3f, %.3f, %.3f)"
                                          %
                                          (iL, "mu" if l.isMu() else "el", l.Eta(), l.Phi(), l.Pt())
                                          for iL, l in enumerate(pre_lep)])
        event_flag = ntevent.evt().cutFlags[0]
        def mll(leps) : return (leps[0] + leps[1]).M()
        def is_emu_or_mue(leps):
            l0, l1 = leps[0], leps[1]
            return ((l0.isEle() and l1.isMu()) or (l0.isMu() and l1.isEle()))
        def pt3020(leps):
            l0, l1 = leps[0], leps[1]
            return ((l0.Pt()>30.0 and l1.Pt()>20.0) or (l0.Pt()>20.0 and l1.Pt()>30.0))
        killHfor = 4
        histos = book_histos(label)
        try:
            cutflow.cut_if(False, 'input')
            cutflow.cut_if(not nttool.passGRL(event_flag), 'grl')
            cutflow.cut_if(not nttool.passLarErr(event_flag), 'larErr')
            cutflow.cut_if(not nttool.passTileErr(event_flag), 'tileErr')
            cutflow.cut_if(not nttool.passTTCVeto(event_flag), 'ttcVeto')
            cutflow.cut_if(not nttool.passGoodVtx(event_flag), 'goodVtx')
            cutflow.cut_if(not nttool.passTileTripCut(event_flag), 'tileTrip')
            cutflow.cut_if(not nttool.passLAr(event_flag), 'lar')
            cutflow.cut_if(not nttool.passBadJet(event_flag), 'bad_jet')
            cutflow.cut_if(not nttool.passBadMuon(event_flag), 'bad_mu')
            cutflow.cut_if(not nttool.passCosmic(event_flag), 'cosmic')
            cutflow.cut_if(not ntevent.evt().hfor==killHfor, 'hfor')
            cutflow.cut_if(not pre_lep.size()==2, '2lep')
            cutflow.cut_if(not sig_taus.size()==0, 'tauveto')
            # passDilTrig includes both trig & match
            cutflow.cut_if(not trig_logic.passDilTrig(pre_lep, met.Et, ntevent.evt()), 'trigger')
            cutflow.cut_if(not mll(pre_lep)>20.0, 'mll20')
            cutflow.cut_if(not is_emu_or_mue(pre_lep), 'emu')
            cutflow.cut_if(not pt3020(pre_lep), 'pt3020')
            fill_histos(histos)
        except SkipEvent:
            continue
    print '\n'+8*'-'+' cutflow '+8*'-'
    print cutflow
    return histos

def book_histos(label):
    "book histos[selection][sample][variable]"
    histos = {}
    selections = ['0j', '1j', '2j', 'vbf2j']
    samples = ['emu', 'mue']
    variables = ['onebin', 'pt0', 'pt1', 'mll', 'mcoll']
    def histo(variable, sam, sel) :
        twopi = +2.0*math.pi
        h = None
        if   v=='onebin'  : h = r.TH1F(histoName(sam, sel, 'onebin' ), ';; entries',                             1, 0.5,   1.5)
        elif v=='pt0'     : h = r.TH1F(histoName(sam, sel, 'pt0'    ), ';p_{T,l0} [GeV]; entries/bin',          12, 0.0, 240.0)
        elif v=='pt1'     : h = r.TH1F(histoName(sam, sel, 'pt1'    ), ';p_{T,l1} [GeV]; entries/bin',          12, 0.0, 240.0)
        elif v=='mll'     : h = r.TH1F(histoName(sam, sel, 'mll'    ), ';m_{l0,l1} [GeV]; entries/bin',         12, 0.0, 240.0)
        elif v=='mcoll'   : h = r.TH1F(histoName(sam, sel, 'mcoll'  ), ';m^{coll}_{l0,l1} [GeV]; entries/bin',  12, 0.0, 240.0)
        elif v=='ptll'    : h = r.TH1F(histoName(sam, sel, 'ptll'   ), ';p_{T,l0+l1} [GeV]; entries/bin',       12, 0.0, 240.0)
        elif v=='dphill'  : h = r.TH1F(histoName(sam, sel, 'dphill' ), ';#Delta#phi(l, l) [rad]; entries/bin',  10, 0.0, twopi)
        elif v=='detall'  : h = r.TH1F(histoName(sam, sel, 'detall' ), ';#Delta#eta(l, l); entries/bin',        10, 0.0, +3.0 )
        elif v=='dphil0met': h= r.TH1F(histoName(sam, sel, 'dphil0met'),';#Delta#phi(l0, met) [rad]; entries/bin',  10, 0.0, twopi)
        elif v=='dphimumet': h= r.TH1F(histoName(sam, sel, 'dphimumet'),';#Delta#phi(#mu, met) [rad]; entries/bin', 10, 0.0, twopi)
        elif v=='mt2j'    : h = r.TH1F(histoName(sam, sel, 'mt2j'   ), ';m^{J}_{T2} [GeV]; entries/bin',        12, 0.0, 480.0)
        elif v=='mljj'    : h = r.TH1F(histoName(sam, sel, 'mljj'   ), ';'+mljjLab+' [GeV]; entries/30GeV',     12, 0.0, 480.0)
        elif v=='dphijj'  : h = r.TH1F(histoName(sam, sel, 'dphijj' ), ';#Delta#phi(j, j) [rad]; entries/bin',  10, 0.0, twopi)
        elif v=='detajj'  : h = r.TH1F(histoName(sam, sel, 'detajj' ), ';#Delta#eta(j, j); entries/bin',        10, 0.0, +3.0 )
        else : print "unknown variable %s"%v
        h.Sumw2()
        h.SetDirectory(0)
        return h
    return dict([(sam, dict([(sel, dict([(v, histo(v, sam, sel)) for v in variables]))
                         for sel in selections]))
                 for sam in samples])

if __name__=='__main__':
    main()
