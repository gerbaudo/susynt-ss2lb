#!/bin/env python

# Root utility functions
#
# davide.gerbaudo@gmail.com
# 2013-08-26

import os
try:
    import numpy as np
except ImportError:
    print "missing numpy: some functions will not be available"

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

from utils import verticalSlice

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
def getMinMaxFromTGraph(gr) :
    points = range(gr.GetN())
    y = np.array([gr.GetY()[p] for p in points])
    return min(y), max(y)
def getMinMaxFromTGraphAsymmErrors(gr) :
    points = range(gr.GetN())
    y    = np.array([gr.GetY()[p] for p in points])
    y_el = np.array([abs(gr.GetErrorYlow (i)) for i in points])
    y_eh = np.array([abs(gr.GetErrorYhigh(i)) for i in points])
    return min(y-y_el), max(y+y_eh)
def getMinMaxFromTH1(h) :
    bins = range(1, 1+h.GetNbinsX())
    y   = np.array([h.GetBinContent(b) for b in bins])
    y_e = np.array([h.GetBinError(b)   for b in bins])
    return min(y-y_e), max(y+y_e)
def getMinMax(histosOrGraphs=[]) :
    def mM(obj) :
        cname = obj.Class().GetName()
        if   cname.startswith('TH1') :               return getMinMaxFromTH1(obj)
        elif cname.startswith('TGraphAsymmErrors') : return getMinMaxFromTGraphAsymmErrors(obj)
        elif cname.startswith('TGraph') :            return getMinMaxFromTGraph(obj)
    ms, Ms = verticalSlice([mM(o) for o in histosOrGraphs])
    return min(ms), max(Ms)
def buildRatioHistogram(num, den, name='', divide_opt='B') :
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
def summedHisto(histos) :
    "return an histogram that is the sum of the inputs"
    hsum = histos[0].Clone(histos[0].GetName()+'_sum')
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

def topRightLabel(pad, label, xpos=None, ypos=None, align=33) :
    pad.cd()
    tex = r.TLatex(0.0, 0.0, '')
    tex.SetNDC()
    tex.SetTextAlign(align)
    tex.DrawLatex((1.0-pad.GetRightMargin()) if not xpos else xpos, (1.0-pad.GetTopMargin()) if not ypos else ypos, label)
    pad._label = tex
    return tex
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
