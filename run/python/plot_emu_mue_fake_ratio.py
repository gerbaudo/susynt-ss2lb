#!/bin/env python

# plot the emu/mue ratio in m_coll
# required input: fake m_coll distributions from plot_emu.py
#
# davide.gerbaudo@gmail.com
# Oct 2014

import ROOT as r
r.gROOT.SetBatch(1)
r.gStyle.SetPadTickX(1)
r.gStyle.SetPadTickY(1)

def main():
    fname = '/tmp/NOM_fake.root'
    input_file = r.TFile.Open(fname)
    h_emu = input_file.Get('h_mcoll_fake_sr_mue_os')
    h_mue = input_file.Get('h_mcoll_fake_sr_emu_os')

    c = r.TCanvas('c', '')
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
    xax = pm.GetXaxis()
    yax = pm.GetYaxis()
    xax.SetLabelSize(0)
    xax.SetTitleSize(0)
    leg = topLeftLegend(c, 0.275, 0.275)
    leg.SetBorderSize(0)
    leg.SetHeader('non-prompt prediction')
    leg.AddEntry(h_mue, '#mue', 'L')
    leg.AddEntry(h_emu, 'e#mu', 'L')
    leg.Draw()
    c._graphics.append(leg)
    c.cd()
    botPad.Draw()
    botPad.cd()
    h_r = h_emu.Clone('h_r')
    h_r.Divide(h_mue)
    h_r.SetLineColor(r.kBlack)
    h_r.SetTitle('')
    xax = h_r.GetXaxis()
    yax = h_r.GetYaxis()
    yax.SetRangeUser(0.0, 2.0)
    yax.SetTitle('#mue / e#mu ratio')
    xA, yA = h_r.GetXaxis(), h_r.GetYaxis()
    textScaleUp = 1.0/botPad.GetHNDC()
    yA.SetNdivisions(-104)
    yA.CenterTitle()
    yA.SetTitleOffset(yA.GetTitleOffset()/textScaleUp)
    for a in [xA, yA] :
        a.SetLabelSize(a.GetLabelSize()*textScaleUp)
        a.SetTitleSize(a.GetTitleSize()*textScaleUp)
    h_r.Draw()
    refline = referenceLine(xA.GetXmin(), xA.GetXmax())
    refline.Draw()
    c._graphics.append(refline)
    c.Update()
    c.SaveAs('ratio_emu_mue.png')

    
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
def referenceLine(xmin=0., xmax=100.0, ymin=1.0, ymax=1.0) :
    l1 = r.TLine(xmin, ymin, xmax, ymax)
    l1.SetLineStyle(3)
    l1.SetLineColor(r.kGray+1)
    return l1

if __name__=='__main__':
    main()
