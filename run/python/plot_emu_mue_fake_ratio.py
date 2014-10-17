#!/bin/env python

# plot the emu/mue ratio in m_coll
# required input: fake m_coll distributions from plot_emu.py
#
# davide.gerbaudo@gmail.com
# Oct 2014

import os
import sys

from rootUtils import (buildBotTopPads
                       ,drawAtlasLabel
                       ,getMinMax
                       ,increaseAxisFont
                       ,referenceLine
                       ,topLeftLegend
                       ,importRoot
                       ,setAtlasStyle
                       )
r = importRoot()
setAtlasStyle()
import systUtils
from systUtils import buildTotBackgroundHisto, buildStatisticalErrorBand, buildFakeSystematicErrorBand
import utils

def main():
    if len(sys.argv)!=3:
        print "Usage: {0} inputdir outputdir".format(sys.argv[0])
        return
    inputdir = sys.argv[1]
    outputdir = sys.argv[2]
    verbose = True
    if not os.path.exists(inputdir):
        print "missing input dir {0}".format(inputdir)
        return
    utils.mkdirIfNeeded(outputdir)

    fake = systUtils.Group('fake')
    fake.setHistosDir(inputdir)
    fake.exploreAvailableSystematics(verbose)
    fakeSystematics = [s for s in fake.systematics if s!='NOM']

    c = r.TCanvas('c','')
    var = 'mcoll'

    h_emu = fake.getHistogram(variable=var, selection='sr_emu_os', cacheIt=True)
    h_mue = fake.getHistogram(variable=var, selection='sr_mue_os', cacheIt=True)
    h_ratio = h_mue.Clone(h_mue.GetName().replace('mue', 'mue_over_emu'))
    h_ratio.Divide(h_emu)
    plot_emu_mue_with_ratio(canvas=c, h_mue=h_mue, h_emu=h_emu, h_ratio=h_ratio,
                            filename="{0}/emu_over_mue_wout_sys_err.png".format(outputdir))

    h_with_totErrBand = {} # histo with stat+syst err (to get the correct error in the ratio)
    for sel in regions_to_plot():
        fake.setSystNominal()
        nominalHistoData    = None
        nominalHistoFakeBkg = fake.getHistogram(variable=var, selection=sel, cacheIt=True)
        nominalHistosBkg    = {'fake', nominalHistoFakeBkg}
        nominalHistoTotBkg  = buildTotBackgroundHisto(histoFakeBkg=nominalHistoFakeBkg, histosSimBkgs={})
        statErrBand = buildStatisticalErrorBand(nominalHistoTotBkg)
        systErrBand = buildFakeSystematicErrorBand(fake=fake, nominalHistosSimBkg={}, variable=var, selection=sel,
                                                   variations=fakeSystematics, verbose=verbose)
        totErrBand = systUtils.addErrorBandsInQuadrature(statErrBand, systErrBand)
        c.cd()
        c.Clear()
        nominalHistoFakeBkg.Draw()
        totErrBand.Draw('E2 same')
        totErrBand.SetFillStyle(3005)
        c.SaveAs("{0}/{1}_{2}.png".format(outputdir, sel, var))
        h_with_totErrBand[sel] = systUtils.setHistErrFromErrBand(nominalHistoFakeBkg, totErrBand)
    h_emu = [h for k,h in h_with_totErrBand.iteritems() if 'emu' in k][0]
    h_mue = [h for k,h in h_with_totErrBand.iteritems() if 'mue' in k][0]

    h_ratio = h_mue.Clone(h_mue.GetName().replace('mue', 'mue_over_emu'))
    h_ratio.Divide(h_emu)
    plot_emu_mue_with_ratio(canvas=c, h_mue=h_mue, h_emu=h_emu, h_ratio=h_ratio,
                            filename="{0}/emu_over_mue_with_sys_err.png".format(outputdir))
    return

def plot_emu_mue_with_ratio(canvas=None, h_mue=None, h_emu=None, h_ratio=None, filename=''):
    c = canvas
    c.Clear()
    c._graphics = []
    c.cd()
    botPad, topPad = buildBotTopPads(c, splitFraction=0.375)
    topPad.Draw()
    topPad.cd()
    h_mue.SetStats(0)
    h_emu.SetStats(0)
    h_emu.SetLineColor(r.kBlue)
    h_mue.SetLineColor(r.kRed)
    h_mue.Draw()
    h_emu.Draw('same')
    pm = h_mue
    yMin, yMax = getMinMax([h_mue, h_emu])
    pm.SetMinimum(0.0)
    pm.SetMaximum(1.1 * yMax)
    xax = pm.GetXaxis()
    yax = pm.GetYaxis()
    xax.SetLabelSize(0)
    xax.SetTitleSize(0)
    leg = topLeftLegend(c, 0.275, 0.275, shift=0.05)
    leg.SetBorderSize(0)
    leg.SetHeader('non-prompt prediction')
    leg.AddEntry(h_mue, '#mue', 'L')
    leg.AddEntry(h_emu, 'e#mu', 'L')
    leg.Draw()
    c._graphics.append(leg)
    c.cd()
    botPad.Draw()
    botPad.cd()
    h_ratio.SetLineColor(r.kBlack)
    h_ratio.SetTitle('')
    xax = h_ratio.GetXaxis()
    yax = h_ratio.GetYaxis()
    yax.SetRangeUser(0.0, 2.0)
    yax.SetTitle('#mue / e#mu ratio')
    xA, yA = h_ratio.GetXaxis(), h_ratio.GetYaxis()
    textScaleUp = 1.0/botPad.GetHNDC()
    yA.SetNdivisions(-104)
    yA.CenterTitle()
    yA.SetTitleOffset(yA.GetTitleOffset()/textScaleUp)
    for a in [xA, yA] :
        a.SetLabelSize(a.GetLabelSize()*textScaleUp)
        a.SetTitleSize(a.GetTitleSize()*textScaleUp)
    h_ratio.Draw()
    refline = referenceLine(xA.GetXmin(), xA.GetXmax())
    refline.Draw()
    c._graphics.append(refline)
    c.Update()
    c.SaveAs(filename)

def regions_to_plot() : return ['sr_mue_os', 'sr_emu_os']

if __name__=='__main__':
    main()
