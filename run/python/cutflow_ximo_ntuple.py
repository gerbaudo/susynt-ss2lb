#!/bin/env python


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
from utils import (commonPrefix
                   ,commonSuffix
                   ,dictSum
                   ,first
                   ,getCommandOutput
                   ,mkdirIfNeeded
                   ,filterWithRegexp
                   ,remove_duplicates
                   ,sortedAs
                   )
from indexed_chain import IndexedChain
import settings
import utils

susyntutils = utils.import_susyntutils()
R = susyntutils.import_root()
susyntutils.load_packages()

from CutflowTable import CutflowTable


def main() :
    usage = """
simple cutflow, run with -h
"""
    parser = optparse.OptionParser(usage=usage)
    parser.add_option('-i', '--input')
    parser.add_option('--debug', action='store_true')
    parser.add_option('--verbose', action='store_true')
    (opts, args) = parser.parse_args()

    treename = 'AnaNtup'
    chain = R.TChain(treename)
    R.ChainHelper.addInput(chain, opts.input, opts.verbose)
    print "Number of events in the input chain: {0}".format(chain.GetEntries())

    eventlist = None
    eventIndices = xrange(chain.GetEntries())
    # if True:
    #     eventlist = R.TEventList("ss3lEventlist", "events selected")
    #     chain.Draw(">>"+eventlist.GetName(),
    #                "RunNb==270806 && "
    #                "(EventNumber==5600684 ||"
    #                " EventNumber==7314718 ||"
    #                " EventNumber==8903766 )",
    #                'groff')
    #     eventIndices = [eventlist.GetList()[i] for i in xrange(eventlist.GetN())]
    # for iEvent, iEntry in enumerate(eventIndices):
    for iEvent, event in enumerate(chain):
        # chain.GetEntry(iEntry)
        # if chain.RunNb!=270806 or chain.EventNumber not in [5600684, 7314718, 8903766]:
        if (chain.RunNb!=270806 or
            chain.EventNumber not in [7684284, 16022745, 4073193, 14768930, 16595938]):
            continue
        print "run {0} event {1}".format(chain.RunNb, chain.EventNumber)
        electrons = Electron.build_electrons(chain)
        muons = Muon.build_muons(chain)
        jets = Jet.build_jets(chain)

        print_muons(muons)
        print_electrons(electrons)
        print_jets(jets)

        # if iEvent>10:
        #     break


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

# # for reference: branch names v18 (the ones that are not in Electron, Muon, Jet)
# L1_2EM15
# L1_2EM7
# L1_2MU4
# L1_2MU6
# L1_EM12
# L1_EM15
# L1_EM15_MU4
# L1_EM7_MU10
# L1_MU10
# L1_MU15
# L1_MU4
# L1_XE35
# HLT_e24_lhmedium_iloose_L1EM18VH
# HLT_e24_lhmedium_iloose_L1EM20VH
# HLT_e24_lhmedium_L1EM18VH
# HLT_e24_lhtight_iloose
# HLT_e24_lhtight_iloose_L1EM20VH
# HLT_e24_medium_iloose_L1EM18VH
# HLT_e24_medium_iloose_L1EM20VH
# HLT_e24_medium_L1EM18VH
# HLT_e24_tight_iloose
# HLT_e24_tight_iloose_L1EM20VH
# HLT_e26_lhtight_iloose
# HLT_e26_tight_iloose
# HLT_e60_lhmedium
# HLT_e60_medium
# HLT_j100_xe80
# HLT_j80_xe80
# HLT_mu18
# HLT_mu20_iloose_L1MU15
# HLT_mu24_iloose_L1MU15
# HLT_mu24_imedium
# HLT_mu26_imedium
# HLT_mu50
# HLT_x700_pufit_wEFMu
# HLT_xe100
# HLT_xe100_mht
# HLT_xe100_mht_wEFMu
# HLT_xe100_pueta
# HLT_xe100_pueta_wEFMu
# HLT_xe100_pufit
# HLT_xe100_pufit_wEFMu
# HLT_xe100_wEFMu
# HLT_xe35
# HLT_xe35_mht
# HLT_xe35_mht_wEFMu
# HLT_xe35_pueta
# HLT_xe35_pueta_wEFMu
# HLT_xe35_pufit
# HLT_xe35_pufit_wEFMu
# HLT_xe35_wEFMu
# HLT_xe70
# HLT_xe70_mht
# HLT_xe70_mht_wEFMu
# HLT_xe70_pueta
# HLT_xe70_pueta_wEFMu
# HLT_xe70_pufit
# HLT_xe70_wEFMu
# HLT_2e12_lhloose_cutd0dphideta_L12EM10VH
# HLT_2e12_lhloose_L12EM10VH
# HLT_2e12_lhloose_nod0_L12EM10VH
# HLT_2e12_lhloose_nodeta_L12EM10VH
# HLT_2e12_lhloose_nodphires_L12EM10VH
# HLT_2e12_loose_L12EM10VH
# HLT_2e17_lhloose
# HLT_2e17_loose
# HLT_mu18_mu8noL1
# HLT_mu20_mu8noL1
# HLT_mu22_mu8noL1
# HLT_mu24_mu8noL1
# HLT_e17_lhloose
# HLT_e17_loose
# HLT_e17_lhloose_mu14
# HLT_e17_loose_mu14
# HLT_e17_lhloose_nod0_mu14
# HLT_mu40
# HLT_2mu10
# HLT_2mu14
# HLT_2e15_lhloose_L12EM13VH
# EventNumber
# ChannelNumber
# AvgMu
# EventWeight
# bcid
# LB
# passGRL
# RunNb
# DetError
# Etmiss_CST_Etx
# Etmiss_CST_Ety
# Etmiss_CST_Et
# Etmiss_TST_Etx
# Etmiss_TST_Ety
# Etmiss_TST_Et
# PV_z



#-----------------------------------------------------------
class Electron(object):
    "wraps electron from flat ntuple"

    def __init__(self, event, index):
        for attribute in Electron.attribute_names():
            prefix = Electron.prefix()
            setattr(self, attribute, getattr(event, prefix+attribute)[index])
        MeV2GeV = 1.0e-3
        self.pT *= MeV2GeV
        self.E *= MeV2GeV

    @staticmethod
    def size_leaf():
        return 'NEl'

    @staticmethod
    def prefix():
        return 'El_'

    @staticmethod
    def attribute_names():
        return [
            'pT'
            ,'eta'
            ,'etaclus'
            ,'phi'
            ,'E'
            ,'charge'
            ,'sigd0'
            # ,'sigd0old'
            ,'z0pvtx'
            # ,'d0pvtx'
            # ,'SFwMedium'
            # ,'SFwMediumLH'
            # ,'SFwTight'
            # ,'SFwTightLH'
            ,'isLooseLH'
            ,'isMediumLH'
            ,'isTightLH'
            # ,'isLoose'
            # ,'isMedium'
            # ,'isTight'
            # ,'nBLayerHits'
            # ,'expectBLayerHit'
            # ,'ptcone20'
            # ,'ptcone30'
            # ,'ptcone40'
            # ,'ptvarcone20'
            # ,'ptvarcone30'
            # ,'ptvarcone40'
            # ,'topoetcone20'
            # ,'topoetcone30'
            # ,'topoetcone40'
            # ,'passIsoLooseTO'
            # ,'passIsoLoose'
            # ,'passIsoTight'
            # ,'passIsoGrad'
            # ,'passIsoGradLoose'
        ]
    @staticmethod
    def build_electrons(event):
        n_electrons = getattr(event, Electron.size_leaf())
        print 'n_electrons ',n_electrons
        return [Electron(event, i) for i in xrange(n_electrons)]
    @staticmethod
    def column_width():
        return 5

    @staticmethod
    def header():
        cell = '%'+str(Electron.column_width())+'s'
        return ' '.join([cell % a for a in Electron.attribute_names()])

    def str(self):
        cell = '%'+str(Electron.column_width())+'.2f'
        return ' '.join([cell % getattr(self, a) for a in Electron.attribute_names()])

def print_electrons(electrons):
    print 'electrons: '+('--' if not electrons else '')
    if electrons:
        w = Electron.column_width()
        print w*' '+Electron.header()
        template_line = '[%'+str(w)+'s] %s'
        print '\n'.join([template_line % (i, m.str()) for i, m in enumerate(electrons)])

#-----------------------------------------------------------
class Muon(object):
    "wraps muons from flat ntuple"

    def __init__(self, event, index):
        for attribute in Muon.attribute_names():
            prefix = Muon.prefix()
            setattr(self, attribute, getattr(event, prefix+attribute)[index])
        MeV2GeV = 1.0e-3
        self.pT *= MeV2GeV

    @staticmethod
    def size_leaf():
        return 'NMu'

    @staticmethod
    def prefix():
        return 'Mu_'

    @staticmethod
    def attribute_names():
        return [
            'pT'
            ,'eta'
            ,'phi'
            # ,'SFw'
            # ,'charge'
            # ,'d0pvtx'
            ,'sigd0'
            # ,'sigd0old'
            ,'z0pvtx'
            # ,'isBad'
            # ,'isCosmic'
            # ,'ptcone20'
            # ,'ptcone30'
            # ,'ptcone40'
            # ,'ptvarcone20'
            # ,'ptvarcone30'
            # ,'ptvarcone40'
            # ,'topoetcone20'
            # ,'topoetcone30'
            # ,'topoetcone40'
            # ,'passIsoLooseTO'
            # ,'passIsoLoose'
            # ,'passIsoTight'
            # ,'passIsoGrad'
            # ,'passIsoGradLoose'
            ]
    @staticmethod
    def build_muons(event):
        n_muons = getattr(event, Muon.size_leaf())
        print 'n_muons ',n_muons
        return [Muon(event, i) for i in xrange(n_muons)]

    @staticmethod
    def column_width():
        return 5

    @staticmethod
    def header():
        cell = '%'+str(Muon.column_width())+'s'
        return ' '.join([cell % a for a in Muon.attribute_names()])

    def str(self):
        cell = '%'+str(Muon.column_width())+'.2f'
        return ' '.join([cell % getattr(self, a) for a in Muon.attribute_names()])

def print_muons(muons):
    print 'muons: '+('--' if not muons else '')
    if muons:
        w = Muon.column_width()
        print w*' '+Muon.header()
        template_line = '[%'+str(w)+'s] %s'
        print '\n'.join([template_line % (i, m.str()) for i, m in enumerate(muons)])
#-----------------------------------------------------------
class Jet(object):
    "wraps jets from flat ntuple"

    def __init__(self, event, index):
        for attribute in Jet.attribute_names():
            prefix = Jet.prefix()
            setattr(self, attribute, getattr(event, prefix+attribute)[index])
        MeV2GeV = 1.0e-3
        self.pT *= MeV2GeV
        self.E *= MeV2GeV

    @staticmethod
    def size_leaf():
        return 'NJet'

    @staticmethod
    def prefix():
        return 'Jet_'

    @staticmethod
    def attribute_names():
        return [
            'eta'
            ,'phi'
            ,'pT'
            ,'E'
            # ,'quality'
            # ,'JVF'
            # ,'JVT'
            # ,'MV2c20'
            # ,'SFw'
            ,'nTrk'
            ]

    @staticmethod
    def build_jets(event):
        n_jets = getattr(event, Jet.size_leaf())
        print 'n_jets ',n_jets
        return [Jet(event, i) for i in xrange(n_jets)]

    @staticmethod
    def column_width():
        return 5

    @staticmethod
    def header():
        cell = '%'+str(Jet.column_width())+'s'
        return ' '.join([cell % a for a in Jet.attribute_names()])

    def str(self):
        cell = '%'+str(Jet.column_width())+'.2f'
        return ' '.join([cell % getattr(self, a) for a in Jet.attribute_names()])

def print_jets(jets):
    print 'jets: '+('--' if not jets else '')
    if jets:
        w = Jet.column_width()
        print w*' '+Jet.header()
        template_line = '[%'+str(w)+'s] %s'
        print '\n'.join([template_line % (i, m.str()) for i, m in enumerate(jets)])


if __name__=='__main__' :
    main()
