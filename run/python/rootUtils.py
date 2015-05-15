#!/bin/env python

# Root utility functions
#
# davide.gerbaudo@gmail.com
# 2013-08-26

import os
import math

try:
    import numpy as np
except ImportError:
    print "missing numpy: some functions will not be available"
    class HackyNumpy(object):
        "a hack replacement for np; it is only expected to work for basic array operations"
        def __init__(self):
            pass
        def array(self, an_array):
            length = an_array.GetNoElements() if hasattr(an_array, 'GetNoElements') else len(an_array)
            vv = r.TVectorD(length)
            for i in range(length) : vv[i] = an_array[i]
            return vv
        def sqrt(self, an_array):
            return self.array(an_array).Sqrt()
        def square(self, an_array):
            return self.array(an_array).Sqr()
    np = HackyNumpy()

def get_np():
    return np

def importRoot() :
    import ROOT as r
    r.gROOT.SetBatch(True)                     # no windows popping up
    r.PyConfig.IgnoreCommandLineOptions = True # don't let root steal our cmd-line options
    r.gErrorIgnoreLevel = r.kWarning           # avoid annoying 'Info in <TCanvas::Print>'
    return r
r = importRoot()

def importRootCorePackages() :
    "same functionality as RootCore/scripts/load_packages.C"
    rcoreDir = os.environ['ROOTCOREDIR']
    [r.gSystem.Load(l.strip()) for l in open(os.path.join(rcoreDir, 'preload'))]
    [r.gSystem.Load('lib%s'%l.strip()) for l in open(os.path.join(rcoreDir, 'load'))]

from utils import verticalSlice, commonPrefix, commonSuffix

def referenceLine(xmin=0., xmax=100.0, ymin=1.0, ymax=1.0) :
    l1 = r.TLine(xmin, ymin, xmax, ymax)
    l1.SetLineStyle(3)
    l1.SetLineColor(r.kGray+1)
    return l1
def firstHisto(histos) :
    return (histos.itervalues().next() if type(histos) is dict
            else histos[0] if type(histos) is list
            else None)
def unitLineFromFirstHisto(histos) :
    fH = firstHisto(histos)
    xAx = fH.GetXaxis()
    return referenceLine(xAx.GetXmin(), xAx.GetXmax())

def topLeftLegend(pad,  legWidth, legHeight, shift=0.0) :
    rMarg, lMarg, tMarg = pad.GetRightMargin(), pad.GetLeftMargin(), pad.GetTopMargin()
    leg = r.TLegend(0.0 + lMarg + shift,
                    1.0 - tMarg - legHeight + shift,
                    0.0 + rMarg + legWidth + shift,
                    1.0 - tMarg + shift)
    leg.SetBorderSize(1)
    leg.SetFillColor(0)
    leg.SetFillStyle(0)
    pad._leg = leg
    return leg
def topRightLegend(pad,  legWidth, legHeight, shift=0.0) :
    rMarg, lMarg, tMarg = pad.GetRightMargin(), pad.GetLeftMargin(), pad.GetTopMargin()
    leg = r.TLegend(1.0 - rMarg - legWidth + shift,
                    1.0 - tMarg - legHeight + shift,
                    1.0 - rMarg + shift,
                    1.0 - tMarg + shift)
    leg.SetBorderSize(1)
    leg.SetFillColor(0)
    leg.SetFillStyle(0)
    pad._leg = leg
    return leg
def rightLegend(pad) :
    bMarg, rMarg, lMarg, tMarg = pad.GetBottomMargin(), pad.GetRightMargin(), pad.GetLeftMargin(), pad.GetTopMargin()
    leg = r.TLegend(1.0 - rMarg, bMarg, 1.0, 1.0 - tMarg)
    leg.SetBorderSize(1)
    leg.SetFillColor(0)
    leg.SetFillStyle(0)
    pad._leg = leg
    return leg
def drawLegendWithDictKeys(pad, histosDict, legWidth=0.325, legHeight=0.225, opt='p') :
    leg = topRightLegend(pad, legWidth, legHeight)
    for s,h in histosDict.iteritems() :
        leg.AddEntry(h, s, opt)
    leg.Draw()
    pad.Update()
    return leg
def dummyHisto(name='dummy', title='', n=1, x_min=0.0, x_max=1.0):
    "dummy histogram, sometime useful for empty lines in a TLegend"
    h = r.TH1F(name, title, n, x_min, x_max)
    h.SetFillColor(0)
    h.SetLineWidth(0)
    h.SetMarkerSize(0)
    h.SetDirectory(0)
    return h
def getMinMaxFromTGraph(gr) :
    points = range(gr.GetN())
    y = np.array([gr.GetY()[p] for p in points])
    return (min(y), max(y)) if points else (0.0, 0.0)
def getMinMaxFromTGraphAsymmErrors(gr) :
    points = range(gr.GetN())
    y    = np.array([gr.GetY()[p] for p in points])
    y_el = np.array([abs(gr.GetErrorYlow (i)) for i in points])
    y_eh = np.array([abs(gr.GetErrorYhigh(i)) for i in points])
    return (min(y-y_el), max(y+y_eh)) if points else (0.0, 0.0)
def getMinMaxFromTH1(h) :
    bins = range(1, 1+h.GetNbinsX())
    y   = np.array([h.GetBinContent(b) for b in bins])
    y_e = np.array([h.GetBinError(b)   for b in bins])
    return (min(y-y_e), max(y+y_e)) if bins else (0.0, 0.0)
def getMinMax(histosOrGraphs=[]) :
    def mM(obj) :
        cname = obj.Class().GetName()
        if   cname.startswith('TH1') :               return getMinMaxFromTH1(obj)
        elif cname.startswith('TGraphAsymmErrors') : return getMinMaxFromTGraphAsymmErrors(obj)
        elif cname.startswith('TGraph') :            return getMinMaxFromTGraph(obj)
    ms, Ms = verticalSlice([mM(o) for o in histosOrGraphs if o])
    return min(ms), max(Ms)
def buildRatioHistogram(num, den, name='', divide_opt='B') :
    "Given two histograms, provide one that is the bin-by-bin ratio"
    ratio = num.Clone(name if name else num.GetName()+'_over_'+den.GetName())
    ratio.SetDirectory(0) # we usually don't care about the ownership of these temporary objects
    ratio.Reset()
    ratio.Divide(num, den, 1, 1, divide_opt)
    return ratio
def getNumDenHistos(f, baseHname='base_histo_name', suffNum='_num', suffDen='_den') :
    "baseHname is something like 'lep_controlreg_chan_var', see MeasureFakeRate2::initHistos()"
    num = f.Get(baseHname+suffNum)
    den = f.Get(baseHname+suffDen)
    return {'num':num, 'den':den}
def buildBotTopPads(canvas, splitFraction=0.275, squeezeMargins=True) :
    canvas.cd()
    botPad = r.TPad(canvas.GetName()+'_bot', 'bot pad', 0.0, 0.0, 1.0, splitFraction, 0, 0, 0)
    interPadMargin = 0.5*0.05
    botPad.SetTopMargin(interPadMargin)
    botPad.SetBottomMargin(botPad.GetBottomMargin()/splitFraction)
    if squeezeMargins : botPad.SetRightMargin(0.20*botPad.GetRightMargin())
    r.SetOwnership(botPad, False)
    canvas.cd()
    canvas.Update()
    topPad = r.TPad(canvas.GetName()+'_top', 'top pad', 0.0, splitFraction, 1.0, 1.0, 0, 0)
    topPad.SetBottomMargin(interPadMargin)
    if squeezeMargins : topPad.SetTopMargin(0.20*topPad.GetTopMargin())
    if squeezeMargins : topPad.SetRightMargin(0.20*topPad.GetRightMargin())
    r.SetOwnership(topPad, False)
    canvas._pads = [topPad, botPad]
    return botPad, topPad

def summedHisto(histos=[], label='sum', remove_extra_underscore=True):
    "return an histogram that is the sum of the inputs"
    def guess_out_name(input_names=[]):
        pref = commonPrefix(hnames)
        suff = commonSuffix(hnames)
        out_name = (pref+label+suff) if (pref or suff) else (input_names[0]+'_'+label)
        out_name = out_name.replace('__','_') if remove_extra_underscore else out_name
        out_name = out_name[:-1] if remove_extra_underscore and out_name.endswith('_') else out_name
        return out_name
    hnames = [h.GetName() for h in histos]
    htemplate = histos[0]
    hsum = htemplate.Clone(guess_out_name(hnames))
    hsum.Reset()
    for h in histos : hsum.Add(h)
    return hsum

def binContentsWithUoflow(h) :
    nBinsX = h.GetNbinsX()+1
    return [h.GetBinContent(0)] + [h.GetBinContent(i) for i in range(1, nBinsX)] + [h.GetBinContent(nBinsX+1)]

def cloneAndFillHisto(histo, bincontents=[], suffix='', zeroErr=True) :
    h, bc= histo, bincontents
    assert h.GetNbinsX()==len(bc),"%d bincontents for %d bins"%(len(bc), h.GetNbinsX())
    h = h.Clone(h.GetName()+suffix)
    for i, c in enumerate(bc) :
        h.SetBinContent(i+1, c)
        if zeroErr : h.SetBinError(i+1, 0.)
    return h

def cumEffHisto(histoTemplate, bincontents=[], leftToRight=True) :
    h, bc= histoTemplate, bincontents
    assert h.GetNbinsX()==len(bc),"%d bincontents for %d bins"%(len(bc), h.GetNbinsX())
    h = h.Clone(h.GetName()+'_ce')
    tot = bc[-1] if leftToRight else bc[0]
    for i, c in enumerate(bc) :
        h.SetBinContent(i+1, c/tot if tot else 0.)
        h.SetBinError(i+1, 0.)
    h.SetMinimum(0.0)
    h.SetMaximum(1.0)
    h.SetTitle('')
    h.SetFillStyle(0)
    return h
def maxSepVerticalLine(hSig, hBkg, yMin=0.0, yMax=1.0) :
    nxS, nxB = hSig.GetNbinsX(), hBkg.GetNbinsX()
    assert nxS==nxB,"maxSepVerticalLine : histos with differen binning (%d!=%d)"%(nxS,nxB)
    bcS = [hSig.GetBinContent(i) for i in range(1,1+nxS)]
    bcB = [hBkg.GetBinContent(i) for i in range(1,1+nxB)]
    def indexMaxDist(bcS, bcB) :
        return sorted([(i,d) for i,d in enumerate([abs(a-b) for a,b in zip(bcS, bcB)])],
                      key= lambda x : x[1])[-1][0]
    iMax = indexMaxDist(bcS, bcB)
    xPos = hSig.GetBinCenter(iMax+1)
    sep = [abs(a-b) for a,b in zip(bcS, bcB)]
    return r.TLine(xPos, yMin, xPos, yMax)

def topRightLabel(pad, label, xpos=None, ypos=None, align=33, scale=1.0) :
    pad.cd()
    tex = r.TLatex(0.0, 0.0, '')
    tex.SetNDC()
    tex.SetTextAlign(align)
    tex.SetTextSize(scale*tex.GetTextSize())
    tex.DrawLatex((1.0-pad.GetRightMargin()) if not xpos else xpos, (1.0-pad.GetTopMargin()) if not ypos else ypos, label)
    pad._label = tex
    return tex

def drawAtlasLabel(pad, xpos=None, ypos=None, align=33, scale=1.0) :
    label = "#bf{#it{ATLAS}} Internal, #sqrt{s} = 8 TeV, 20.3 fb^{-1}"
    return topRightLabel(pad, label, xpos, ypos, align, scale)

def getBinning(h) :
    cname = h.Class().GetName()
    if   cname.startswith('TH1') : return (h.GetNbinsX())
    elif cname.startswith('TH2') : return (h.GetNbinsX(), h.GetNbinsY())
    elif cname.startswith('TH3') : return (h.GetNbinsX(), h.GetNbinsY(), h.GetNbinsZ())
    else : return None

def getBinIndices(h) :
    "Return a list of the internal indices used by TH1/TH2/TH3; see TH1::GetBin for info on internal mapping"
    cname = h.Class().GetName()
    if   cname.startswith('TH1') :
        return [h.GetBin(i)
                for i in range(1, 1+h.GetNbinsX())]
    elif cname.startswith('TH2') :
        return [h.GetBin(i, j)
                for i in range(1, 1+h.GetNbinsX())
                for j in range(1, 1+h.GetNbinsY())]
    elif cname.startswith('TH3') :
        return [h.GetBin(i, j, k)
                for i in range(1, 1+h.GetNbinsX())
                for j in range(1, 1+h.GetNbinsY())
                for k in range(1, 1+h.GetNbinsZ())]
    else : return []
def getBinCenters(h) :
    return [h.GetBinCenter(b) for b in getBinIndices(h)]
def getBinContents(h) :
    return [h.GetBinContent(b) for b in getBinIndices(h)]
def getBinErrors(h) :
    return [h.GetBinError(b) for b in getBinIndices(h)]


def writeObjectsToFile(outputFileName='', objects={}, verbose=False):
    """
    Objects can either be a list or a dict of TObjects (or a TObject).
    The list/dict will be recursively navigated down to TObjects.
    """
    outputFile = r.TFile.Open(outputFileName, 'recreate')
    outputFile.cd()
    if verbose : print "writing to %s"%outputFile.GetName()
    def write(dictOrObj):
        isDict = type(dictOrObj) is dict
        isList = type(dictOrObj) is list
        if isDict:
            for v in dictOrObj.values():
                write(v)
        elif isList:
            for v in dictOrObj:
                write(v)
        else:
            dictOrObj.Write()
    write(objects)
    outputFile.Close()
def fetchObjectsFromFile(fileName='', objects={}, verbose=False, closeFileOnExit=False):
    """
    Objects is either a name, or a list of names, or a dict of names
    (or a nested dict,list).  This function will return the same
    structure (that is list, dict, nested dict...), but with objects
    instead of names.
    """
    fetched_objects = type(objects)()
    inputFile = r.TFile.Open(fileName)
    def fetch(object):
        isDict = type(object) is dict
        isList = type(object) is list
        if isDict:
            return dict([(k, fetch(v)) for k,v in object.iteritems()])
        elif isList:
            return [fecth(o) for o in object]
        else:
            obj = inputFile.Get(object)
            if verbose and not obj : print "cannot get '%s' from '%s'"%(str(object), inputFile.GetName())
            return obj
    if not inputFile : print "cannot open %s"%fileName
    else:
        if verbose : print "fetching histograms from %s"%inputFile.GetName()
        fetched_objects = fetch(objects)
        if closeFileOnExit : inputFile.Close()
    return fetched_objects
def printHistoIntegrals(histos={}, keyColWidth=20, nameColWidth=40, integralColWidth=20, integralPrecision=1):
    "given a dict of histos, print a table with their integrals"
    templateLine = '{0:<'+str(keyColWidth)+'} {1:<'+str(nameColWidth)+'} : {2:>'+str(integralColWidth)+'.'+str(integralPrecision)+'f}'
    print ('integrals :\n'
           +'\n'.join([templateLine.format(k, h.GetName(), h.Integral())
                       for k, h in histos.iteritems()]))

def integralAndError(h) :
    error = r.Double(0.0)
    integral = h.IntegralAndError(0,-1, error)
    return integral, float(error)

def setAtlasStyle() :
    aStyle = getAtlasStyle()
    r.gROOT.SetStyle("ATLAS")
    r.gROOT.ForceStyle()

def getAtlasStyle() :
    style = r.TStyle('ATLAS', 'Atlas style')
    white = 0
    style.SetFrameBorderMode(white)
    style.SetFrameFillColor(white)
    style.SetCanvasBorderMode(white)
    style.SetCanvasColor(white)
    style.SetPadBorderMode(white)
    style.SetPadColor(white)
    style.SetStatColor(white)
    #style.SetPaperSize(20,26)
    style.SetPadTopMargin(0.05)
    style.SetPadRightMargin(0.05)
    nokidding = 0.75 # limit the exaggerated margins
    style.SetPadBottomMargin(nokidding*0.16)
    style.SetPadLeftMargin(nokidding*0.16)
    style.SetTitleXOffset(nokidding*1.4)
    style.SetTitleYOffset(nokidding*1.4)
    font, fontSize = 42, 0.04 # helvetica large
    style.SetTextFont(font)
#     style.SetTextSize(fontSize)
    style.SetLabelFont(font,"xyz")
    style.SetTitleFont(font,"xyz")
    style.SetPadTickX(1)
    style.SetPadTickY(1)
    style.SetOptStat(0)
    style.SetOptTitle(0)
    style.SetEndErrorSize(0)
    return style

def increaseAxisFont(axis, factorLabel=1.25, factorTitle=1.25) :
    axis.SetLabelSize(factorLabel*axis.GetLabelSize())
    axis.SetTitleSize(factorTitle*axis.GetTitleSize())

def graphWithPoissonError(histo, fillZero=False) :
    "From TGuiUtils.cxx; no idea where this implementation is coming from. Ask Anyes et al."
    gr = r.TGraphAsymmErrors()
    gr.SetLineWidth(histo.GetLineWidth())
    gr.SetLineColor(histo.GetLineColor())
    gr.SetLineStyle(histo.GetLineStyle())
    gr.SetMarkerSize(histo.GetMarkerSize())
    gr.SetMarkerColor(histo.GetMarkerColor())
    gr.SetMarkerStyle(histo.GetMarkerStyle())
    binCenters  = getBinCenters(histo)
    binContents = getBinContents(histo)
    def poissonErr(n) :
        sqrt = math.sqrt
        err_up, err_do = 0.0, 0.0
        if n :
            y1, y2 = n+1.0, n
            d1 = 1.0 - 1.0/(9.0*y1) + 1.0/(3.0*sqrt(y1))
            d2 = 1.0 - 1.0/(9.0*y2) - 1.0/(3.0*sqrt(y2))
            err_up = y1*d1*d1*d1 - n
            err_do = n - y2*d2*d2*d2;
        return err_do, err_up
    xErr = 0.0
    yErrors = [poissonErr(v) for v in binContents]
    for bc, bv, (ed, eu) in zip(binCenters, binContents, yErrors) :
        if bv or fillZero :
            point = gr.GetN()
            gr.SetPoint(point, bc, bv)
            gr.SetPointError(point, xErr, xErr, ed, eu)
    histo._poissonErr = gr # attach to histo for persistency
    return gr

def getXrange(h) :
    nbins = h.GetNbinsX()
    x_lo = h.GetBinCenter(1) - 0.5*h.GetBinWidth(1)
    x_hi = h.GetBinCenter(nbins) + 0.5*h.GetBinWidth(nbins)
    return x_lo, x_hi

def computeStatErr2(nominal_histo=None) :
    "Compute the bin-by-bin err2 (should include also mc syst, but for now it does not)"
    bins = range(1, 1+nominal_histo.GetNbinsX())
    bes = [nominal_histo.GetBinError(b)   for b in bins]
    be2s = np.array([e*e for e in bes])
    return {'up' : be2s, 'down' : be2s}

def buildErrBandGraph(histo_tot_bkg, err2s) :
    h = histo_tot_bkg
    bins = range(1, 1+h.GetNbinsX())
    x = np.array([h.GetBinCenter (b) for b in bins])
    y = np.array([h.GetBinContent(b) for b in bins])
    ex_lo = ex_hi = np.array([0.5*h.GetBinWidth(b) for b in bins])
    ey_lo, ey_hi = np.sqrt(err2s['down']), np.sqrt(err2s['up'])
    gr = r.TGraphAsymmErrors(len(bins), x, y, ex_lo, ex_hi, ey_lo, ey_hi)
    gr.SetMarkerSize(0)
    gr.SetFillStyle(3004)
    gr.SetFillColor(r.kGray+3)
    gr.SetLineWidth(2)
    return gr
