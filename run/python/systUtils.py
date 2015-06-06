# Utility functions to handle systematic variations
#
# davide.gerbaudo@gmail.com
# March 2014

# todo : replace the various range(...) with rootUtils.getBinIndices()

import collections
import math
import os
from rootUtils import importRoot
r = importRoot()

import rootUtils
from rootUtils import integralAndError, getBinContents

try:
    import numpy as np
except ImportError:
    np = rootUtils.np

from utils import filterWithRegexp, first, remove_duplicates, sortedAs

def fakeSystWeightLeaves() :
    """
    Names of the leaves used to store the syst variations for the
    fake.  See WeightVariations.h and TupleMaker::initTreeBranches()
    The names of the systematics are the ones specified in
    susy::fake::Systematic.
    """
    branch = 'relWeights.'
    return {'EL_RE_UP'    : branch+'fakeElRealUp',
            'EL_RE_DOWN'  : branch+'fakeElRealDo',
            'MU_RE_UP'    : branch+'fakeElFakeUp',
            'MU_RE_DOWN'  : branch+'fakeElFakeDo',
            'EL_FR_UP'    : branch+'fakeMuRealUp',
            'EL_FR_DOWN'  : branch+'fakeMuRealDo',
            'MU_FR_UP'    : branch+'fakeMuFakeUp',
            'MU_FR_DOWN'  : branch+'fakeMuFakeDo',
            'EL_FRAC_DOWN': branch+'fakeElFracDo',
            'EL_FRAC_UP'  : branch+'fakeElFracUp',
            'MU_FRAC_DOWN': branch+'fakeMuFracDo',
            'MU_FRAC_UP'  : branch+'fakeMuFracUp'
            }
def fakeSystVariations() :
    return fakeSystWeightLeaves().keys()

def mcObjectVariations() :
    "See definitions in SusyDefs.h:SusyNtSystNames, and active list in SusyPlotter::toggleStdSystematics()"
    return ['EERDOWN', 'EERUP',
            'EESLOWDOWN', 'EESLOWUP',
            'EESMATDOWN', 'EESMATUP',
            'EESPSDOWN', 'EESPSUP',
            'EESZDOWN', 'EESZUP',
            'JER',
            'JESDOWN', 'JESUP',
            'MESDOWN', 'MESUP',
            'MIDDOWN', 'MIDUP',
            'RESOST',
            'SCALESTDOWN', 'SCALESTUP',
            ]

def mcWeightLeaves() :
    """
    Leaves where the relative weight variations for weight systematics
    are stored. See WeightVariations.h
    """
    branch = 'relWeights.'
    return {'ETRIGREWUP'  : branch+'elTrigUp',
            'ETRIGREWDOWN': branch+'elTrigDo',
            'MTRIGREWUP'  : branch+'muTrigUp',
            'MTRIGREWDOWN': branch+'muTrigDo',
            'BJETUP'      : branch+'bTagUp',
            'BJETDOWN'    : branch+'bTagDo',
            'CJETUP'      : branch+'cTagUp',
            'CJETDOWN'    : branch+'cTagDo',
            'BMISTAGUP'   : branch+'lTagUp',
            'BMISTAGDOWN' : branch+'lTagDo',
            'ESFUP'       : branch+'elEffUp',
            'ESFDOWN'     : branch+'elEffDo',
            'MEFFUP'      : branch+'muEffUp',
            'MEFFDOWN'    : branch+'muEffDo',
            }

def mcWeightVariations() :
    return mcWeightLeaves().keys()

def getAllVariations() :
    return ['NOM'] + fakeSystVariations() + mcObjectVariations() + mcWeightVariations()

def fetchVariationHistos(input_fake_file=None, nominal_histo=None, variations=fakeSystVariations()) :
    nom_hname = nominal_histo.GetName()
    return dict([(v, input_fake_file.Get(nom_hname.replace('_NONE','_'+v))) for v in variations])
def computeSysErr2(nominal_histo=None, vars_histos={}) :
    """Add in quadrature the the bin-by-bin delta^2 for up&down
    systematic variations.  We do not know, and we do not care,
    whether the ones labeled as 'up' give a larger yield and the
    'down' ones a smaller yield. So we lump all the positive
    variations together, and do the same for negative ones.
    This gives the most conservative error estimate.
    Return the err^2 for up and down.
    """
    def bc(h) : return [h.GetBinContent(b) for b in range(1, 1+h.GetNbinsX())]
    nom_bcs  = bc(nominal_histo)
    vars_bcs = dict([(v, bc(h)) for v, h in vars_histos.iteritems()])
    bins = range(nominal_histo.GetNbinsX())
    deltas = [[vars_bcs[v][b] - nom_bcs[b] for v in vars_bcs.keys()] for b in bins]
    def positive(ll) : return [l if not l<0.0 else 0.0 for l in ll ]
    def negative(ll) : return [l if     l<0.0 else 0.0 for l in ll ]
    def sumquad(ll) : return sum([l*l for l in ll])
    up_e2s = np.array([sumquad(positive(deltas[b])) for b in bins])
    do_e2s = np.array([sumquad(negative(deltas[b])) for b in bins])
    return {'up' : up_e2s, 'down' : do_e2s}
def computeStatErr2(nominal_histo=None) :
    "Compute the bin-by-bin stat err2"
    bins = range(1, 1+nominal_histo.GetNbinsX())
    bes = [nominal_histo.GetBinError(b)   for b in bins]
    be2s = np.array([e*e for e in bes])
    return {'up' : be2s, 'down' : be2s}
def computeFakeSysStatErr2(nominal_histo=None, vars_histos={}) :
    "Compute the bin-by-bin sum2 err including up&down fake systematic variations + stat. unc."
    print "refactor, use computeStatErr2 and computeSysErr2"
    def bc(h) : return [h.GetBinContent(b) for b in range(1, 1+h.GetNbinsX())]
    def be(h) : return [h.GetBinError(b) for b in range(1, 1+h.GetNbinsX())]
    nom_bcs  = bc(nominal_histo)
    nom_be2s = [e*e for e in be(nominal_histo)]
    vars_bcs = dict([(v, bc(h)) for v, h in vars_histos.iteritems()])
    bins = range(nominal_histo.GetNbinsX())
    deltas = [[vars_bcs[v][b] - nom_bcs[b] for v in vars_bcs.keys()] for b in bins]
    def positive(ll) : return [l if not l<0.0 else 0.0 for l in ll ]
    def negative(ll) : return [l if     l<0.0 else 0.0 for l in ll ]
    def sumquad(ll) : return sum([l*l for l in ll])
    up_e2s = np.array([sumquad(positive(deltas[b])) for b in bins]) + np.array(nom_be2s)
    do_e2s = np.array([sumquad(negative(deltas[b])) for b in bins]) + np.array(nom_be2s)
    return {'up' : up_e2s, 'down' : do_e2s}
def fetchFakeSysHistosAndComputeSysErr2(input_fake_file=None, nominal_histo=None) :
    vars_histos = fetchVariationHistos(input_fake_file, nominal_histo)
    return computeFakeSysStatErr2(nominal_histo, vars_histos)
def buildErrBandGraph(histo_tot_bkg, err2s) :
    h = histo_tot_bkg
    bins = range(1, 1+h.GetNbinsX())
    x = np.array([h.GetBinCenter (b) for b in bins])
    y = np.array([h.GetBinContent(b) for b in bins])
    ex_lo = ex_hi = np.array([0.5*h.GetBinWidth(b) for b in bins])
    ey_lo, ey_hi = np.sqrt(err2s['down']), np.sqrt(err2s['up'])
    using_tvectors = hasattr(x, 'ClassName')
    gr = (r.TGraphAsymmErrors(x, y, ex_lo, ex_hi, ey_lo, ey_hi) if using_tvectors else
          r.TGraphAsymmErrors(len(bins), x, y, ex_lo, ex_hi, ey_lo, ey_hi))
    gr.SetMarkerSize(0)
    gr.SetFillStyle(3004)
    gr.SetFillColor(r.kGray+3)
    gr.SetLineWidth(2)
    return gr
def buildErrBandRatioGraph(errband_graph) :
    gr = errband_graph.Clone()
    points = range(gr.GetN())
    xs     = np.array([gr.GetX()[i] for i in points])
    ys     = np.array([gr.GetY()[i] for i in points])
    eys_lo = np.array([abs(gr.GetErrorYlow (i)) for i in points])
    eys_hi = np.array([abs(gr.GetErrorYhigh(i)) for i in points])
    ys_lo  = ys - eys_lo
    ys_hi  = ys + eys_hi
    def absFracDeltaFromUnity(y_nom, y_var) : return abs(y_var/y_nom - 1.0) if y_nom else 0.0
    eys_lo = [absFracDeltaFromUnity(n, v) for n, v in zip(ys, ys_lo)]
    eys_hi = [absFracDeltaFromUnity(n, v) for n, v in zip(ys, ys_hi)]
    for p, x, ey_lo, ey_hi in zip(points, xs, eys_lo, eys_hi) :
        gr.SetPoint(p, x, 1.0) # TGraph does not have a SetPointY, so we need to set both x and y
        gr.SetPointEYlow (p, ey_lo)
        gr.SetPointEYhigh(p, ey_hi)
    return gr

def buildAsymmRatioGraph(num_data, den_reference, name='') :
    """
    Given two TGraphAsymmErrors, provide one that is the bin-by-bin ratio.

    This function is meant to be used when the errors on the input
    histograms are asymmetric (e.g. poisson errors from
    graphWithPoissonError). The main case is for the data points in
    the bottom pad showing the ratio observed/expected.

    For symmetric errors, prefer rootUtils.buildRatioHistogram().

    Note to self: I assume that the denominator is the 'reference'
    (i.e. expected total background) that will be represented as an
    error band centered around 1. Therefore, this function takes into
    account only the error on the numerator.
    """
    num, den = num_data, den_reference
    valid_p_num   = [i for i in range(num.GetN())]
    nonzero_p_den = [i for i in range(den.GetN()) if den.GetX()[i]]
    valid_x_num   = [num.GetX()[i] for i in valid_p_num]
    nonzero_x_den = [den.GetX()[i] for i in nonzero_p_den]
    common_x = sorted(list(set(valid_x_num).intersection(nonzero_x_den)))
    points_num = [i for i,x in zip(valid_p_num, valid_x_num) if x in common_x]
    points_den = [i for i,x in zip(nonzero_p_den, nonzero_x_den) if x in common_x]
    assert len(points_num)==len(points_den),("num and den: unmatched number of points: num[{}]!=den[{}], {}{}"
                                             .format(len(points_num), len(points_den), num.GetName(), den.GetName()))
    gr = num.Clone()
    gr.SetName("{}_over_{}".format(num.GetName(), den.GetName()))
    [gr.RemovePoint(i) for i in range(gr.GetN())] # keep the formatting, then add only the points we want
    xs     = np.array([num.GetX()[i] for i in points_num])
    ys     = np.array([num.GetY()[i] for i in points_num])
    ys_ref = np.array([den.GetY()[i] for i in points_den])
    eys_lo = np.array([abs(num.GetErrorYlow (i)) for i in points_num])
    eys_hi = np.array([abs(num.GetErrorYhigh(i)) for i in points_num])
    ys_lo  = ys - eys_lo
    ys_hi  = ys + eys_hi
    for x, y, y_lo, y_hi, y_ref in zip(xs, ys, ys_lo, ys_hi, ys_ref):
        p = gr.GetN()
        central = y/y_ref
        ey_lo = abs(y_lo/y_ref - central)
        ey_hi = abs(y_hi/y_ref - central)
        gr.SetPoint      (p, x, central) # TGraph does not have a SetPointY, so we need to set both x and y
        gr.SetPointEYlow (p, ey_lo)
        gr.SetPointEYhigh(p, ey_hi)
    return gr

def setHistErrFromErrBand(h, errband_graph):
    "given an histogram and a graph, set the error of the hist to the error from the graph"
    assert h.GetNbinsX()==errband_graph.GetN()
    points = range(errband_graph.GetN())
    # histos don't support asymmetric errors, take the average
    eys_lo = np.array([abs(errband_graph.GetErrorYlow (i)) for i in points])
    eys_hi = np.array([abs(errband_graph.GetErrorYhigh(i)) for i in points])
    def avg(v1, v2) : return np.array([0.5*(e1+e2) for e1, e2 in zip(v1, v2)])
    errs = avg(eys_lo, eys_hi)
    for i in points:
        h.SetBinError(i+1, errs[i])
    return h
#___________________________________________________________
class BaseSampleGroup(object) :
    """
    Base class for a sample or a group of samples.

    Holds the info about the type of sample (data/mc/fake, etc.) and
    is aware of which systematic variation should be processed for
    each sample. Also keep track of the yield variation for each
    systematic.
    """
    def __init__(self, name) :
        self.name = name
        self.setSystNominal()
        self.varCounts = collections.defaultdict(dict)
        self.currentSelection = None # used only for filenameHisto
    @property
    def label(self) : return self.groupname if hasattr(self, 'groupname') else self.name
    @property
    def isFake(self) : return self.label=='fake'
    @property
    def isData(self) : return self.label=='data'
    @property
    def isMc(self) : return not (self.isFake or self.isData)
    @property
    def isSignal(self) : return 'signal' in self.label
    @property
    def isMcBkg(self) : return self.isMc and not self.isSignal
    def isNeededForSys(self, sys) :
        return (sys=='NOM'
                or (self.isMc and sys in mcWeightVariations())
                or (self.isMc and sys in mcObjectVariations())
                or (self.isFake and sys in fakeSystVariations()))
    def setSystNominal(self) : return self.setSyst()
    def setSyst(self, sys='NOM') :
        "Set the syst; if we should not consider this syst for the current sample, set to nominal"
        nominal = 'NOM' # do we have differnt names for nom (mc vs fake)?
        self.isObjSys    = sys in mcObjectVariations()
        self.isWeightSys = sys in mcWeightVariations() + fakeSystVariations()
        self.isFakeSys   = sys in fakeSystVariations()
        def nameObjectSys(s) : return s if self.isMc else nominal
        def nameWeightSys(s) : return s if self.isMc else nominal
        def nameFakeSys(s) : return s if self.isFake else nominal
        def identity(s) : return s
        sysNameFunc = (nameObjectSys if self.isObjSys else
                       nameWeightSys if self.isWeightSys else
                       nameFakeSys if self.isFakeSys else
                       identity)
        self.syst = sysNameFunc(sys)
        self.syst = sys
        return self
    def setCurrentSelection(self, sel=''):
        self.currentSelection = sel
        return self
    def logVariation(self, sys='', selection='', counts=0.0) :
        "log this systematic variation and internally store it as [selection][sys]"
        self.varCounts[selection][sys] = counts
        return self
    def variationsSummary(self) :
        summaries = {} # one summary for each selection
        for selection, sysCounts in self.varCounts.iteritems() :
            nominalCount = sysCounts['NOM']
            summaries[selection] = [(sys, sysCount, (100.0*(sysCount-nominalCount)/nominalCount)
                                     if nominalCount else None)
                                    for sys, sysCount in sortedAs(sysCounts, getAllVariations())]
        return summaries
    def printVariationsSummary(self):
        for selection, summarySel in self.variationsSummary().iteritems() :
            colW = str(12)
            header = ' '.join([('%'+colW+'s')%colName for colName in ['variation', 'yield', 'delta[%]']])
            lineTemplate = '%(sys)'+colW+'s'+'%(counts)'+colW+'s'+'%(delta)'+colW+'s'
            print "---- summary of variations for %s ----" % self.name
            print "---     [selection: %s]            ---" % selection
            print header
            print '\n'.join(lineTemplate%{'sys':s,
                                          'counts':(("%.3f"%c) if type(c) is float
                                                    else (str(c)+str(type(c)))),
                                          'delta' :(("%.3f"%d) if type(d) is float
                                                    else '--' if d==None
                                                    else (str(d)+str(type(d)))) }
                            for s,c,d in summarySel)
    @property
    def weightLeafname(self) :
        leafname = 'event.pars.weight'
        if  self.isWeightSys:
            rel_weight = ""
            if self.isFakeSys: rel_weight = "event.{0}".format(fakeSystWeightLeaves()[self.syst])
            elif self.isWeightSys: rel_weight = "event.{0}".format(mcWeightLeaves()[self.syst])
            else: print "not implemented yet"
            leafname += " * {0}".format(rel_weight)
        return leafname
    @classmethod
    def histoname(cls, sample='', syst='', selection='', variable=''):
        return "h_%s_%s_%s_%s"%(variable, sample, syst, selection)


def findByName(bsgs=[], name='') : return [b for b in bsgs if b.name==name][0]
#___________________________________________________________
class Sample(BaseSampleGroup) :
    "A sample that is aware of its input ntuples and event weight formulas"
    def __init__(self, name, groupname) :
        super(Sample, self).__init__(name) # this is either the name (for data and fake) or the dsid (for mc)
        self.groupname = groupname
        self.setHftInputDir()
    @property
    def filename(self):
        "filename for the ntuple holding the current syst (without directory)"
        return (self.name+'.root' if not self.isObjSys else
                self.name+'_'+self.syst+'.root')
    def setHftInputDir(self, dir='') :
        useDefaults = not dir
        defaultDir = 'out/fakepred' if self.isFake else 'out/susyplot'
        self.hftInputDir = defaultDir if useDefaults else dir
        return self
    @property
    def filenameHftTree(self) :
        def dataFilename(sample, dir, sys) : return "%(dir)s/%(sys)s_%(sam)s.PhysCont.root" % {'dir':dir, 'sam':sample, 'sys':sys}
        def fakeFilename(sample, dir, sys) : return "%(dir)s/%(sys)s_fake.%(sam)s.PhysCont.root" % {'dir':dir, 'sam':sample, 'sys':sys}
        def mcFilename  (sample, dir, sys) : return "%(dir)s/%(sys)s_%(dsid)s.root" % {'dir':dir, 'sys':sys, 'dsid':sample}
        fnameFunc = dataFilename if self.isData else fakeFilename if self.isFake else mcFilename
        sys = (self.syst
               if (self.isMc and self.isObjSys or self.isFake and self.isFakeSys)
               else 'NOM')
        return fnameFunc(self.name, self.hftInputDir, sys)
    @property
    def hftTreename(self) :
        def dataTreename(samplename) : return "id_%(s)s.PhysCont" % {'s' : samplename}
        def fakeTreename(samplename) : return "id_fake.%(s)s.PhysCont"%{'s':samplename}
        def mcTreename(dsid=123456) :  return "id_%d"%dsid
        getTreename = dataTreename if self.isData else fakeTreename if self.isFake else mcTreename
        return getTreename(self.name)
    def hasInputHftFile(self, msg) :
        filename = self.filenameHftTree
        isThere = os.path.exists(filename)
        if not isThere : print msg+"%s %s missing : %s"%(self.groupname, self.name, filename)
        return isThere
    def hasInputHftTree(self, msg='') :
        treeIsThere = False
        if self.hasInputHftFile(msg) :
            filename, treename = self.filenameHftTree, self.hftTreename
            inputFile = r.TFile.Open(filename) if self.hasInputHftFile(msg) else None
            if inputFile :
                if inputFile.Get(treename) : treeIsThere = True
                else : print msg+"%s %s missing tree '%s' from %s"%(self.groupname, self.name, treename, filename)
            inputFile.Close()
        return treeIsThere
    def group(self) :
        return Group(self.groupname).setSyst(self.syst)
#___________________________________________________________
class Group(BaseSampleGroup) :
    "A group of sameples that's aware of the root files with the histograms to be plotted"
    def __init__(self, name) :
        super(Group, self).__init__(name)
        self.setSyst()
        self.setHistosDir()
        self._histoCache = collections.defaultdict(dict) # [syst][histoname]
    def setHistosDir(self, dir='') :
        self.histosDir = dir if dir else 'out/hft'
        return self
    @property
    def filenameHisto(self):
        "file containig the histograms for the current syst"
        fname = "%(dir)s/%(group)s_%(sys)s.root" % {'group':self.name, 'dir':self.histosDir, 'sys':self.syst}
        if self.currentSelection:
            fname = "%(dir)s/%(group)s_%(sys)s_%(sel)s.root" % {'group':self.name, 'dir':self.histosDir,
                                                                'sys':self.syst, 'sel':self.currentSelection}
        return fname
    def exploreAvailableSystematics(self, verbose=False) :
        systs = ['NOM']
        if self.isFake :
            systs += fakeSystVariations()
        elif self.isMc :
            systs += mcObjectVariations()
            systs += mcWeightVariations()
        self.systematics = []
        for sys in systs :
            self.setSyst(sys)
            if os.path.exists(self.filenameHisto) :
                self.systematics.append(sys)
        if verbose : print "%s : found %d variations : %s"%(self.name, len(self.systematics), str(self.systematics))
    def filterAndDropSystematics(self, include='.*', exclude=None, verbose=False) :
        "include and exclude can be either a regex, a single value, or a list"
        nBefore = len(self.systematics)
        def is_regex(exp) : return exp and '*' in exp
        def is_list(exp) : return exp and ',' in exp
        def is_single_value(exp) : return exp and len(exp)
        def str_to_list(exp) : return eval("['{0}']".format(exp))
        toBeIncluded = ([s for s in self.systematics if s in str_to_list(include)] if is_list(include) else
                        filterWithRegexp(self.systematics, include) if is_regex(include) else
                        str_to_list(include) if is_single_value(include) else
                        self.systematics)
        toBeExcluded = ([s for s in self.systematics if s in str_to_list(include)] if is_list(include) else
                        filterWithRegexp(self,systematics, exclude) if is_regex(exclude) else
                        [])
        self.systematics = remove_duplicates([s for s in toBeIncluded if s not in toBeExcluded])
        nAfter = len(self.systematics)
        if verbose : print "%s : dropped %d systematics, left with %s"%(self.name, nBefore-nAfter, str(self.systematics))
        assert self.systematics.count('NOM')==1 or not nBefore, "%s : 'NOM' required %s"%(self.name, str(self.systematics))
    def getHistogram(self, variable, selection, cacheIt=False) :
        hname = BaseSampleGroup.histoname(sample=self.name, syst=self.syst,
                                          selection=selection, variable=variable)
        histo = None
        try :
            histo = self._histoCache[self.syst][hname]
        except KeyError :
            self.setCurrentSelection(selection)
            file = r.TFile.Open(self.filenameHisto)
            if not file : print "missing file %s"%self.filenameHisto
            histo = file.Get(hname)
            if not histo : print "%s : cannot get histo %s from %s"%(self.name, hname, self.filenameHisto)
            elif cacheIt :
                histo.SetDirectory(0)
                self._histoCache[self.syst][hname] = histo
                file.Close()
            else :
                histo.SetDirectory(0)
                file.Close()
        if variable=='onebin' and histo : self.logVariation(self.syst, selection, histo.Integral(0, -1))
        return histo
    def getBinContents(self, variable, selection) :
        return getBinContents(self.getHistogram)
#___________________________________________________________
def buildTotBackgroundHisto(histoFakeBkg=None, histosSimBkgs={}) :
    hTemplate = histoFakeBkg
    if not hTemplate :
        hTemplate = first(histosSimBkgs)
        print "warning, cannot use fake as template; histo name will be wrong, based on %s"%hTemplate.GetName()
    totBkg = hTemplate.Clone(hTemplate.GetName().replace('_fake','_totbkg'))
    totBkg.Reset()
    if not totBkg.GetSumw2N() : totBkg.Sumw2()
    allBackgrounds = dict((group, histo) for group, histo in [('fake',histoFakeBkg)]+[(g,h) for g,h in histosSimBkgs.iteritems()])
    for group, histo in allBackgrounds.iteritems() :
        if histo :
            totBkg.Add(histo)
            integral, error = integralAndError(histo)
            print "adding %.3f +/- %.3f  : %s (%s)" % (integral, error, group, histo.GetName())
        else : print "buildStatisticalErrorBand: missing %s"%group
    integral, error = integralAndError(totBkg)
    print "totBkg %.3f +/- %.3f  : %s" % (integral, error, totBkg.GetName())
    return totBkg
#___________________________________________________________
def buildStatisticalErrorBand(histoTotBkg= None) :
    return buildErrBandGraph(histoTotBkg, computeStatErr2(histoTotBkg))
#___________________________________________________________
def buildSystematicErrorBand(fake=None, simBkgs=[], variable='', selection='',
                             fakeVariations = fakeSystVariations(),
                             mcVariations = mcObjectVariations()+mcWeightVariations(),
                             verbose=False) :
    "build the syst error band accounting for fake, mc-object, and mc-weight systematics"
    if verbose : print "buildSystematicErrorBand(%s, %s)"%(variable, selection)
    for g in [fake] + simBkgs : g.setSystNominal()
    nominalFake   = fake.getHistogram(variable, selection)
    nominalOthers = dict([(g.name, g.getHistogram(variable=variable, selection=selection)) for g in simBkgs])

    fakeSysErrorBand  = buildFakeSystematicErrorBand(fake, nominalOthers, variable, selection, fakeVariations, verbose)
    mcSysErrorBand    = buildMcSystematicErrorBand  (nominalFake, simBkgs, variable, selection, mcVariations, verbose)

    totErrorBand = None
    if fakeSysErrorBand and mcSysErrorBand : totErrorBand = addErrorBandsInQuadrature(fakeSysErrorBand, mcSysErrorBand)
    elif fakeSysErrorBand : totErrorBand = fakeSysErrorBand
    elif mcSysErrorBand : totErrorBand = mcSysErrorBand
    return totErrorBand
#___________________________________________________________
def buildFakeSystematicErrorBand(fake=None, nominalHistosSimBkg={},
                                 variable='', selection='', variations=[], verbose=False) :
    if verbose : print "buildFakeSystematicErrorBand(%s, %s), %s"%(variable, selection, str(variations))
    fake.setSystNominal()
    nominalTotBkg = buildTotBackgroundHisto(fake.getHistogram(variable, selection), nominalHistosSimBkg)
    variedTotBkgs = dict()
    for sys in variations :
        if verbose : print "buildFakeSystematicErrorBand(%s)"%sys
        variedTotBkgs[sys] = buildTotBackgroundHisto(fake.setSyst(sys).getHistogram(variable=variable, selection=selection),
                                                     nominalHistosSimBkg)
    err2 = computeSysErr2(nominal_histo=nominalTotBkg, vars_histos=variedTotBkgs)
    return buildErrBandGraph(nominalTotBkg, err2)
#___________________________________________________________
def buildMcSystematicErrorBand(fakeNominalHisto=None, simulatedBackgrounds=[],
                               variable='', selection='', variations=[], verbose=False) :
    if verbose : print "buildMcSystematicErrorBand(%s, %s), %s"%(variable, selection, str(variations))
    for b in simulatedBackgrounds : b.setSystNominal()
    nominalTotBkg = buildTotBackgroundHisto(fakeNominalHisto, dict((g.name, g.getHistogram(variable, selection))
                                                                   for g in simulatedBackgrounds))
    variedTotBkgs = dict()
    for sys in variations :
        if verbose : print "buildMcSystematicErrorBand(%s)"%sys
        variedMcBkgs = dict((g.name, g.setSyst(sys).getHistogram(variable, selection)) for g in simulatedBackgrounds)
        variedTotBkgs[sys] = buildTotBackgroundHisto(fakeNominalHisto, variedMcBkgs)
    fakeErr2 = computeSysErr2(nominal_histo=nominalTotBkg, vars_histos=variedTotBkgs)
    return buildErrBandGraph(nominalTotBkg, fakeErr2)
#___________________________________________________________
def addErrorBandsInQuadrature(errBand1, errBand2) :
    sqrt = math.sqrt
    totErrBand = None
    if errBand1 and errBand2 :
        totErrBand = errBand1.Clone()
        points = range(totErrBand.GetN())
        eys_stat_lo = np.array([abs(errBand1.GetErrorYlow (i)) for i in points])
        eys_stat_hi = np.array([abs(errBand1.GetErrorYhigh(i)) for i in points])
        eys_syst_lo = np.array([abs(errBand2.GetErrorYlow (i)) for i in points])
        eys_syst_hi = np.array([abs(errBand2.GetErrorYhigh(i)) for i in points])
        eys_lo = np.sqrt(np.square(eys_stat_lo) + np.square(eys_syst_lo))
        eys_hi = np.sqrt(np.square(eys_stat_hi) + np.square(eys_syst_hi))
        for p, ey_lo, ey_hi in zip(points, eys_lo, eys_hi) :
            totErrBand.SetPointEYlow (p, ey_lo)
            totErrBand.SetPointEYhigh(p, ey_hi)
    return totErrBand
#___________________________________________________________
def totalUpDownVariation(errBand) :
    points = range(errBand.GetN())
    up   = [abs(errBand.GetErrorYhigh(p)) for p in points]
    down = [abs(errBand.GetErrorYlow (p)) for p in points]
    return sum(up), sum(down)
#___________________________________________________________
