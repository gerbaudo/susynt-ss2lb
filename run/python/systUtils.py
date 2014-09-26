# Utility functions to handle systematic variations
#
# davide.gerbaudo@gmail.com
# March 2014

# todo : replace the various range(...) with rootUtils.getBinIndices()

import math
from rootUtils import importRoot
r = importRoot()

import rootUtils
from rootUtils import integralAndError

try:
    import numpy as np
except ImportError:
    np = rootUtils.np

def fakeSystVariations() :
    "syst variations for the fake estimate, see DiLeptonMatrixMethod::systematic_names"
    return ['EL_RE_UP', 'EL_RE_DOWN', 'MU_RE_UP', 'MU_RE_DOWN',
            'EL_FR_UP', 'EL_FR_DOWN', 'MU_FR_UP', 'MU_FR_DOWN',
            'EL_FRAC_DO', 'EL_FRAC_UP', 'MU_FRAC_DO', 'MU_FRAC_UP']
def mcObjectVariations() :
    "See definitions in SusyDefs.h:SusyNtSystNames, and active list in SusyPlotter::toggleStdSystematics()"
    return ['EES_Z_UP', 'EES_Z_DN',
            'EES_MAT_UP','EES_MAT_DN',
            'EES_PS_UP', 'EES_PS_DN',
            'EES_LOW_UP', 'EES_LOW_DN',
            'EER_UP', 'EER_DN',
            'MS_UP', 'MS_DN',
            'ID_UP', 'ID_DN',
            'JES_UP', 'JES_DN',
            'JER',
            'SCALEST_UP', 'SCALEST_DN',
            'RESOST',
            ]
def mcWeightVariations() :
    "See list at HftFiller::assignWeightVars()"
    return ['BKGMETHODUP' ,'BKGMETHODDOWN'
            ,'ETRIGREWUP' ,'ETRIGREWDOWN'
            ,'MTRIGREWUP' ,'MTRIGREWDOWN'
            ,'ESFUP'      ,'ESFDOWN'
            ,'MEFFUP'     ,'MEFFDOWN'
            ,'BJETUP'     ,'BJETDOWN'
            ,'CJETUP'     ,'CJETDOWN'
            ,'BMISTAGUP'  ,'BMISTAGDOWN'
            ,'XSUP'       ,'XSDOWN'
            ]
def mcWeightBranchname(mcWeightVariation='') : return 'syst_'+mcWeightVariation
def mcWeightBranches() : return [mcWeightBranchname(v) for v in mcWeightVariations()]

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

#___________________________________________________________
# todo : move Group here (fake and simBkgs are Group objects)
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
