#!/bin/env python

import array
import math
import utils
susyntutils = utils.import_susyntutils()
r = susyntutils.import_root()
susyntutils.load_packages()

import glob

def main():
    test_first_1k_events = True
    print_events_as_txt = False # needed to evaluate the matrix on stdin
    base_dir = '/gdata/atlas/gerbaudo/hlfv/take0/SusyntHlfv/run/out/matrix_prediction'
    dir1 = base_dir+'/Jul_26/'
    dir2 = base_dir+'/May_20/'
    tree_name = 'hlfv_tuple'
    chain1 = r.TChain(tree_name)
    for f in sorted(glob.glob(dir1+'/*.root*')) : chain1.Add(f)
    chain2 = r.TChain(tree_name)
    for f in sorted(glob.glob(dir2+'/*.root*')) : chain2.Add(f)
    entries1 =chain1.GetEntries()
    entries2 =chain2.GetEntries()
    print 'chain1 ',entries1
    print 'chain2 ',entries2
    # assert entries1==entries2,"the two chains must have the same entries (matching events), {0}!={1}".format(entries1, entries2)
    min_weight, max_weight = 1.0e6, -1.0e6
    h_w1_vs_w2   = r.TH2F('h_w2_vs_w1',   'fake weight: Jul_26 vs. May_20; w Jul_26; w May_20', 100, -0.75, +0.75, 100, -0.75, +0.75)
    h_w1_vs_w2_s = r.TH2F('h_w2_vs_w1_s', 'fake weight: Jul_26 vs. May_20; w Jul_26; w May_20', 100, -0.75, +0.75, 100, -0.75, +0.75)
    chain1.SetAlias("weight1", "weight")
    chain2.SetAlias("weight2", "weight")
    print 'chain1 ',dir1
    chain1.Scan("runNumber:eventNumber:weight", "", "", 10)
    print 'chain2 ',dir2
    chain2.Scan("runNumber:eventNumber:weight", "", "", 10) 
    chain1.AddFriend(chain2, "May_20")
    chain1.Draw("weight1:May_20.weight2 >> "+h_w1_vs_w2_s.GetName())
    out_file = r.TFile.Open('fake_weight_comparison.root', 'recreate')
    out_file.cd()
    out_tree = r.TTree('fake_weight_tree',
                       'tree to compare fake weights: '
                       +"w1 from {0}".format(dir1)
                       +"w2 from {0}".format(dir2))
    weights = array.array('d', 2*[0.0])
    pts = array.array('d', 2*[0.0])
    isEmu = array.array( 'i', 1*[0] )
    isSameSign = array.array( 'i', 1*[0] )
    out_tree.Branch('weights', weights, 'weights[%d]/D'%len(weights))
    out_tree.Branch('pts', pts, 'pts[%d]/D'%len(pts))
    out_tree.Branch('isEmu', isEmu, 'isEmu[%d]/I'%len(isEmu))
    out_tree.Branch('isSameSign', isSameSign, 'isSameSign[%d]/I'%len(isSameSign))

    for iEntry in xrange(entries1):
        if test_first_1k_events and iEntry>1000 : break
        chain1.GetEntry(iEntry)
        chain2.GetEntry(iEntry)
        r1, e1, w1 = chain1.pars.runNumber, chain1.pars.eventNumber, chain1.pars.weight
        r2, e2, w2 = chain2.pars.runNumber, chain2.pars.eventNumber, chain2.pars.weight
        if r1!=r2 or e1!=e2:
            print "1) run {0}, event {1}, weight {2}".format(r1, e1, w1)
            print "2) run {0}, event {1}, weight {2}".format(r2, e2, w2)
        assert r1==r2,"different run number {0}!={1}".format(r1, r2)
        assert e1==e2,"different evt number {0}!={1}".format(e1, e2)
        min_weight = min_weight if min_weight<min([w1, w2]) else min([w1, w2])
        max_weight = max_weight if max_weight>max([w1, w2]) else max([w1, w2])
        h_w1_vs_w2.Fill(w1, w2)
        weights[0], weights[1] = w1, w2
        l0, l1 = chain1.l0, chain1.l1
        l0 = addTlv(l0)
        l1 = addTlv(l1)
        pts[0], pts[1] = l0.p4.Pt(), l1.p4.Pt()
        isEl0, isMu0 = l0.isEl, l0.isMu
        isEl1, isMu1 = l1.isEl, l1.isMu
        isEmu[0] = int((isEl0 and isMu1) or (isMu1 and isMu0))
        isSameSign[0] = int((l0.charge * l1.charge)>0)
        values = {'run':r1, 'event':e1,
                  'weight':w1,
                  'isEl1': 1 if isEl0 else 0,
                  'isEl2': 1 if isEl1 else 0,
                  'isTight1': 1 if l0.isTight else 0,
                  'isTight2': 1 if l1.isTight else 0,
                  'pt1': pts[0], 'pt2': pts[1],
                  'eta1': l0.p4.Eta(), 'eta2': l1.p4.Eta()
                  }
        if print_events_as_txt and isEmu[0]:
            print "{run} {event} {weight} {isEl1} {isEl2} {isTight1} {isTight2} {pt1} {pt2} {eta1} {eta2}".format(**values)

        out_tree.Fill()
    print "weight: min {0}, max {1}".format(min_weight, max_weight)
    print "output entries: ",out_tree.GetEntries()
    out_tree.Write()
    out_file.Close()
    
    c = r.TCanvas('c_w1_vs_w2_draw','fake weight: Jul_26 vs. May_20')
    c.cd()
    h_w1_vs_w2.Draw('scat')
    c.SaveAs(h_w1_vs_w2.GetName()+'.png')
    c.SaveAs(h_w1_vs_w2.GetName()+'.root')
    
    c = r.TCanvas('c_w1_vs_w2_draw','fake weight: Jul_26 vs. May_20')
    c.cd()
    c.Clear()
    h_w1_vs_w2_s.Draw('scat')
    c.SaveAs(h_w1_vs_w2_s.GetName()+'.png')
    c.SaveAs(h_w1_vs_w2_s.GetName()+'.root')


tlv = r.TLorentzVector
def FourMom2TLorentzVector(fm) :
    l = tlv()
    l.SetPxPyPzE(fm.px, fm.py, fm.pz, fm.E)
    return l
def addTlv(l) :
    if not hasattr(l, 'p4') : l.p4 = FourMom2TLorentzVector(l)
    return l

if __name__=='__main__':
    main()
