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
r = importRoot()
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
    for iEvent, event in enumerate(chain):
        print "run {0} event {1}".format(event.RunNb, event.EventNumber)
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

if __name__=='__main__' :
    main()

# # for reference: branch names v18
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
# NMu
# Mu_eta
# Mu_phi
# Mu_pT
# Mu_SFw
# Mu_charge
# Mu_d0pvtx
# Mu_sigd0
# Mu_sigd0old
# Mu_z0pvtx
# Mu_isBad
# Mu_isCosmic
# Mu_ptcone20
# Mu_ptcone30
# Mu_ptcone40
# Mu_ptvarcone20
# Mu_ptvarcone30
# Mu_ptvarcone40
# Mu_topoetcone20
# Mu_topoetcone30
# Mu_topoetcone40
# Mu_passIsoLooseTO
# Mu_passIsoLoose
# Mu_passIsoTight
# Mu_passIsoGrad
# Mu_passIsoGradLoose
# MuTrigSF_2mu14
# MuTrigSF_mu14
# NEl
# El_eta
# El_etaclus
# El_phi
# El_pT
# El_E
# El_charge
# El_sigd0
# El_sigd0old
# El_z0pvtx
# El_d0pvtx
# El_SFwMedium
# El_SFwMediumLH
# El_SFwTight
# El_SFwTightLH
# El_isLooseLH
# El_isMediumLH
# El_isTightLH
# El_isLoose
# El_isMedium
# El_isTight
# El_nBLayerHits
# El_expectBLayerHit
# El_ptcone20
# El_ptcone30
# El_ptcone40
# El_ptvarcone20
# El_ptvarcone30
# El_ptvarcone40
# El_topoetcone20
# El_topoetcone30
# El_topoetcone40
# El_passIsoLooseTO
# El_passIsoLoose
# El_passIsoTight
# El_passIsoGrad
# El_passIsoGradLoose
# NJet
# Jet_eta
# Jet_phi
# Jet_pT
# Jet_E
# Jet_quality
# Jet_JVF
# Jet_JVT
# Jet_MV2c20
# Jet_SFw
# Jet_nTrk
# Etmiss_CST_Etx
# Etmiss_CST_Ety
# Etmiss_CST_Et
# Etmiss_TST_Etx
# Etmiss_TST_Ety
# Etmiss_TST_Et
# PV_z
