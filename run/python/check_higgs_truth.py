#!/bin/env python
#
# Check whether we're saving the H from the truth record
#
# davide.gerbaudo@gmail.com
# October 2014

import glob
import os
import utils
susynt = utils.import_susyntutils()
r = susynt.import_root()


def main():
    susynt.load_packages()
    susynt.generate_dicts()
    susynt.import_SUSYDefs_enums()

    sample_name = 'Sherpa_CT10_lllnu_WZ'
    input_dir = '/var/tmp/susynt_dev/data/ntup_susy/'
    #input_dir = '/var/tmp/susynt_dev/data/ntup_common/'
    input_files = glob.glob(os.path.join(input_dir, '*.root*'))


    input_dir=('/gdata/atlas/ucintprod/SusyNt/susy_n0155/'
               'group.phys-susy.mc12_8TeV.169670.'
               'MadGraphPythia8_AU2CTEQ6L1_ggH125_taumu.SusyNt.e1903_s1581_s1586_r3658_r3549_p1512_n0155b/')
    input_file = input_dir+'group.phys-susy.049421._000001.susyNt.root'
    chain = r.TChain('susyNt')
    chain.Add(input_file)
    num_entries = chain.GetEntries()
    num_entries_to_process = num_entries if num_entries<1e1 else int(1e1)
    print "About to loop on %d entries"%num_entries_to_process
    run_with_chain(chain, num_entries_to_process)

def run_with_chain(tree, n_max_entries=-1):
    nttool = r.SusyNtTools()
    m_entry = r.Long(-1)
    ntevent = r.Susy.SusyNtObject(m_entry)
    ntevent.ReadFrom(tree)
    isSimplifiedModel = False
    nttool.buildSumwMap(tree, isSimplifiedModel)
    period, useRewUtils = 'Moriond', False
    trig_logic = r.DilTrigLogic(period, useRewUtils)
    n_entries_to_print = 4
    sys = susynt.SusyNtSys.NtSys_NOM
    tauId = susynt.TauID
    tauJetId, tauEleId, tauMuoId = tauId.TauID_loose, tauId.TauID_medium, tauId.TauID_medium
    for iEntry, entry in enumerate(tree):
        m_entry = iEntry
        if n_max_entries>0 and m_entry >= n_max_entries : break
        if iEntry < n_entries_to_print : print 'run ', ntevent.evt().run,' event ',ntevent.evt().event
        for iP, particle in enumerate(ntevent.tpr()):
            print "[{0}] pdg: {1}, pt: {2}".format(iP, particle.pdgId, particle.Pt())
        emu = [p for p in ntevent.tpr() if abs(p.pdgId) in [13,15]]
        if len(emu)==2:
            m_h = (emu[0] + emu[1]).M()
            print "m_h ",m_h

if __name__=='__main__':
    main()
