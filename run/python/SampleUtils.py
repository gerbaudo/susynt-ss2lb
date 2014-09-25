#!/bin/env python

# Utility classes and functions to manage samples and subsamples
#
# davide.gerbaudo@gmail.com
# Jan 2013

import glob, os, re, unittest
from datasets import datasets, allGroups, allDatasets, activeDatasets
from rootUtils import importRoot
r = importRoot()
import utils

colors = {
    'ttbar'       : r.kRed+1,
    'zjets'       : r.kOrange-2,
    'wjets'       : r.kBlue-2,
    'diboson'     : r.kSpring+2,
    'singletop'   : r.kAzure-4,
    'multijet'    : r.kGray,
    'fake'        : r.kGray, # just another name for the same thing
    'heavyflavor' : r.kViolet+1
    }
markers = {
    'ttbar'       : r.kFullTriangleDown,
    'zjets'       : r.kOpenSquare,
    'wjets'       : r.kFullTriangleUp,
    'diboson'     : r.kFullDiamond,
    'singletop'   : r.kOpenCross,
    'multijet'    : r.kFullCross,
    'fake'        : r.kFullCircle, # just another name for the same thing
    'heavyflavor' : r.kFullSquare
    }

allGroups, allDatasets, datasets = allGroups(datasets), allDatasets(datasets), activeDatasets(datasets)

def fastSamplesFromFilenames(filenames=[], verbose=False):
    """Assume a list of filenames with some common prefix (basepath)
    and suffix (tag): strip away the prefix&suffix, and return a list
    of Dataset objects obtained by matching the name to the dataset
    """
    prefix = utils.commonPrefix(filenames)
    suffix = utils.commonSuffix(filenames)
    if verbose : print "fastSamplesFromFilenames: prefix '%s' suffix '%s'"%(prefix, suffix)
    if not prefix or not suffix : print "fastSamplesFromFilenames: warning prefix '%s' suffix '%s'"%(prefix, suffix)
    filenames = [f.replace(prefix, '').replace(suffix, '') for f in filenames]
    dsets = dict([(d.name, d) for d in datasets])
    return [dsets[f] for f in filenames if f in dsets]

def guessGroupFromFilename(filename='') :
    "Guess group from filename, either merged or un-merged"
    print 'guessGroupFromFilename obsolete, use guessSampleFromFilename'
    group = next((g for g in allGroups if any(e in filename for e in [g+'_', g+'.'])), None)
    group = group if group else next((d.group
                                      for d in datasets
                                      if d.name+'.' in filename
                                      or d.name+'_' in filename), None)
    return group
def guessSampleFromFilename(filename='') :
    "Guess sample from filename"
    def longest(lst=[]) : return max(lst, key=len) if lst else None
    def filenameMatchesSample(fn, smp) : return re.search('([_.])'+smp+'([_.])', fn)
    return longest([s for s in allDatasets if filenameMatchesSample(filename, s)])


def isDataSample(samplename) : return 'data' in samplename or 'period' in samplename
def isSigSample(samplename) : return 'WH_' in samplename
def isBkgSample(samplename) : return not isDataSample(samplename) and not isSigSample(samplename)
def guessReqidFromFilename(filename='', verbose=False) :
    match = re.search('mc12\_8TeV\.(\d+)\.', filename)
    if verbose and not match : print "'%s' does not contain mc12\_8TeV\.(\d+)\."%filename
    return match.group(1) if match else None 


def basePathArea() :
    path = os.path.realpath(__file__)
    return path[:path.rfind('SusyTest0')]
def xsReaderDataDir(basePath='') :
    relPath = '/SusyXSReader/data'
    return (basePath if basePath else basePathArea()) + relPath

class ModeAWhDbPar :
    fields = ['ds', 'mc1', 'mn1', 'xsec', 'xsecSys']
    class Entry:
        def __init__(self, line) :
            """parse a line that is expected to be formatted as follow:
            DS MC1 MN1 xsec xsec_sys
            and store what is necessary.
            """
            line = line.strip()
            words = line.split()
            for a,w in zip(ModeAWhDbPar.fields, words) : setattr(self, a, w)
        def valid(self) :
            return all([hasattr(self, a) for a in ModeAWhDbPar.fields]) and self.ds.isdigit()
    def __init__(self) :
        filenames = [xsReaderDataDir()+'/'+'modeA_WH_MC1eqMN2.txt',
                     xsReaderDataDir()+'/'+'modeA_WH_notauhad_MC1eqMN2_DiagonalMatrix.txt'
                     ]
        self.entries = [e for e in [ModeAWhDbPar.Entry(l)
                                    for fn in filenames
                                    for l in open(fn).readlines()]]
        self.entries = filter(lambda e: e.valid(), self.entries)
    def mc1Mn1ByReqid(self, reqid) :
        entry = next(e for e in self.entries if e.ds == reqid)
        return float(entry.mc1), float(entry.mn1)
    def allMc1(self) : return [e.mc1 for e in self.entries]
    def allMn1(self) : return [e.mn1 for e in self.entries]

class ModeAWhDbReqid :
    "Using the filelists, map reqids to samplenames"
    def __init__(self, filenames = []) :
        self.entries = {}
        filelistDir = basePathArea()+'/SusyTest0/run/filelist/'
        filenames = (filenames if filenames
                     else
        glob.glob(filelistDir+'Herwigpp_simplifiedModel_wA_noslep_WH_*Lep_*.txt')
        + glob.glob(filelistDir+'Herwigpp_sM_wA_noslep_notauhad_WH_2Lep_*.txt'))
        for f in filenames :
            rootfile = open(f).read()
            if not rootfile :
                print "warning, emtpy filelist %s"%f
                continue
            reqid  = guessReqidFromFilename(rootfile)
            sample = guessSampleFromFilename(rootfile)
            if reqid and sample :
                assert sample not in self.entries, "Multiple reqids for one sample : %s, %s"%(sample, str([reqid, self.entries[sample]]))
                self.entries[sample] = reqid
            else :
                print "skipping invalid entry reqid='%s', sample='%s' from '%s'"%(reqid, sample, f)
    def reqidBySample(self, sample) :
        return self.entries[sample]
    def sampleByReqid(self, reqid) :
        return next((s for s,r in self.entries.iteritems() if r == reqid), None)

class ModeAWhDbMergedFake2Lreqid :
    "provide fake request ids to emulate merged samples in the 2L (mc1,mn1) plane"
    def __init__(self) :
        pass
    def reqidByMc1Mn1(self, mc1, mn1) :
        x, y = float(mc1), float(mn1)
        if   y<(-x+210) : return 1765700
        elif y<(-x+290) : return 1765800 if y>(x-210) else 1765900
        elif y<(-x+360) : return 1766000 if y>(x-210) else 1766100 if y>(x-290) else 1766200
        else            : return 1766300
        #if y<(x-100) else None # this is messing up also the bottom half? later on...

#
# testing
#
class KnownReqidModeAWhDb(unittest.TestCase) :
    def testMatchAllAvailabeAttrs(self) :
        knownValues = [ ('176574', (130.0 , 0.0))
                       ,('176575', (140.0, 10.0))
                        ]
        db = ModeAWhDbPar()
        for reqid, valuePair in knownValues :
            v1T, v2T = valuePair[0], valuePair[1]
            v1, v2  = db.mc1Mn1ByReqid(reqid)
            self.assertEqual(v1, v1T)
            self.assertEqual(v2, v2T)

class KnownEntriesModeAWhDbReqid(unittest.TestCase) :
    def testMatchAllAvailabeAttrs(self) :
        knownValues = [ ('176584', 'WH_2Lep_11')
                       ,('176581', 'WH_2Lep_8')
                        ]
        db = ModeAWhDbReqid()
        for reqid, sample in knownValues :
            self.assertEqual(reqid, db.reqidBySample(sample))


if __name__ == "__main__":
    unittest.main()
