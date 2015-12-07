#!/bin/env python

# loop over one of Ximo's nutples and count the electrons as they're being classified with Julien's function.

# davide.gerbaudo@gmail.com
# Dec 0215

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
    def int_attribute_names():
        return ['type', 'origin', 'chFlip']

    @staticmethod
    def float_attribute_names():
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
            # ,'isLooseLH'
            # ,'isMediumLH'
            # ,'isTightLH'
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
    @classmethod
    def attribute_names(cls):
        return cls.int_attribute_names() + cls.float_attribute_names()

    @staticmethod
    def build_electrons(event):
        n_electrons = getattr(event, Electron.size_leaf())
        return [Electron(event, i) for i in xrange(n_electrons)]
    @staticmethod
    def column_width():
        return 5

    @staticmethod
    def header():
        cell = '%'+str(Electron.column_width())+'s'
        return ' '.join([cell % a for a in Electron.attribute_names()])

    def str(self):
        cellf = '%'+str(Electron.column_width())+'.2f'
        celli = '%'+str(Electron.column_width())+'.d'
        float_attributes = Electron.float_attribute_names()
        int_attributes = Electron.int_attribute_names()
        return ' '.join([cellf % getattr(self, a) for a in float_attributes]+
                        [celli % getattr(self, a) for a in int_attributes])

def print_electrons(electrons):
    print 'electrons: '+('--' if not electrons else '')
    if electrons:
        w = Electron.column_width()
        print w*' '+Electron.header()
        template_line = '[%'+str(w)+'s] %s'
        print '\n'.join([template_line % (i, m.str()) for i, m in enumerate(electrons)])

#-----------------------------------------------------------

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
    counter_per_type_origin_flag = collections.defaultdict(int)
    for iEvent, event in enumerate(chain):
        electrons = Electron.build_electrons(chain)
        if iEvent<10:
            print "run {0} event {1}".format(chain.RunNb, chain.EventNumber)
            print_electrons(electrons)
        for electron in electrons:
            key = (electron.type, electron.origin, electron.chFlip)
            counter_per_type_origin_flag[key] += 1
    print "Summary of electron classification:"
    print '\t'.join(['type', 'origin', 'qflipflag'])+" : counts"
    summary_lines = ["%3d / %3d / %3d : %d" % (el_type, el_origin, el_qflipflag, counts)
                     for (el_type, el_origin, el_qflipflag),counts in counter_per_type_origin_flag.iteritems()]
    summary_lines = sorted(summary_lines)
    print '\n'.join(summary_lines)


if __name__=='__main__':
    main()
