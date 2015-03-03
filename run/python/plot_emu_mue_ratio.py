#!/bin/env python

# plot the emu/mue ratio in pt1, m_coll for data, data-fake, and mc
# required input: fake m_coll distributions from plot_emu.py
#
# davide.gerbaudo@gmail.com
# Dec 2014

import os
import pprint
import sys

import dataset
from rootUtils import (buildBotTopPads
                       ,drawAtlasLabel
                       ,getMinMax
                       ,increaseAxisFont
                       ,referenceLine
                       ,topRightLegend
                       ,importRoot
                       ,setAtlasStyle
                       )
r = importRoot()
setAtlasStyle()
import systUtils
from systUtils import buildTotBackgroundHisto, buildStatisticalErrorBand, buildFakeSystematicErrorBand
import utils
from utils import (first
                   ,mkdirIfNeeded
                   )

def main():
    if len(sys.argv)!=3:
        print "Usage: {0} inputdir outputdir".format(sys.argv[0])
        return
    inputdir = sys.argv[1]
    outputdir = sys.argv[2]
    selections = regions_to_plot()
    verbose = True
    if not os.path.exists(inputdir):
        print "missing input dir {0}".format(inputdir)
        return
    utils.mkdirIfNeeded(outputdir)


    buildTotBkg = systUtils.buildTotBackgroundHisto
    buildStat = systUtils.buildStatisticalErrorBand
    # buildSyst = systUtils.buildSystematicErrorBand
    variables = ['mcoll', 'pt1']


    groups = dataset.DatasetGroup.build_groups_from_files_in_dir('./samples')#opts.samples_dir)
    groups.append(first([g for g in groups if g.is_data]).clone_data_as_fake())
    plot_groups = [systUtils.Group(g.name) for g in groups]
    for group in plot_groups :
        group.setCurrentSelection(first(selections))
        group.setHistosDir(inputdir).setCurrentSelection(first(selections))
        group.exploreAvailableSystematics(verbose)
    available_systematics = sorted(list(set([s for g in plot_groups for s in g.systematics])))
    # systematics_to_use = get_list_of_syst_to_fill(opts)
    # systematics = [s for s in systematics_to_use if s in available_systematics]
    systematics = available_systematics
    if verbose :
        print "using the following systematics : {0}".format(systematics)
    #     print "missing the following systematics : {0}".format([s for s in systematics_to_use if s not in available_systematics])
    # fakeSystematics = [s for s in systematics if s in systUtils.fakeSystVariations()]
    # mcSystematics = [s for s in systematics if s in systUtils.mcObjectVariations() + systUtils.mcWeightVariations()]

    mkdirIfNeeded(outputdir)
    findByName = systUtils.findByName
    simBkgs = [g for g in plot_groups if g.isMcBkg]
    data = findByName(plot_groups, 'data')
    fake = findByName(plot_groups, 'fake')

    fake = systUtils.Group('fake')
    fake.setHistosDir(inputdir)
    fake.exploreAvailableSystematics(verbose)
    fakeSystematics = [s for s in fake.systematics if s!='NOM']

    print '---> ',fake.syst
    fake.setSyst()
    print '---> ',fake.syst
    c = r.TCanvas('c','')
    variables = ['mcoll', 'pt1']

    fpt1_histos = {'data_fake':{}, 'sim_bkg':{}}
    for var in variables:
        print ">>>plotting ",var
        for g in plot_groups : g.setSystNominal()
        # sr
        h_emu = fake.getHistogram(variable=var, selection='sr_emu_os', cacheIt=True)
        h_mue = fake.getHistogram(variable=var, selection='sr_mue_os', cacheIt=True)
        h_ratio = h_emu.Clone(h_emu.GetName().replace('emu', 'emu_over_mue'))
        h_ratio.Divide(h_mue)
        plot_emu_mue_with_ratio(canvas=c, h_mue=h_mue, h_emu=h_emu, h_ratio=h_ratio,
                                filename=outputdir+'/sr_'+var+'_emu_over_mue_fake',
                                label='SR: fake')
        h_emu = data.getHistogram(variable=var, selection='sr_emu_os', cacheIt=True)
        h_mue = data.getHistogram(variable=var, selection='sr_mue_os', cacheIt=True)
        h_ratio = h_emu.Clone(h_emu.GetName().replace('emu', 'emu_over_mue'))
        h_ratio.Divide(h_mue)
        plot_emu_mue_with_ratio(canvas=c, h_mue=h_mue, h_emu=h_emu, h_ratio=h_ratio,
                                filename=outputdir+'/sr_'+var+'_emu_over_mue_data',
                                label='SR: data')

        h_emu = data.getHistogram(variable=var, selection='sr_emu_os', cacheIt=True).Clone('sr_emu_data_minus_fake')
        h_mue = data.getHistogram(variable=var, selection='sr_mue_os', cacheIt=True).Clone('sr_mue_data_minus_fake')
        h_emu.Add(fake.getHistogram(variable=var, selection='sr_emu_os', cacheIt=True), -1.0)
        h_mue.Add(fake.getHistogram(variable=var, selection='sr_mue_os', cacheIt=True), -1.0)
        h_ratio = h_emu.Clone(h_emu.GetName().replace('emu', 'emu_over_mue'))
        h_ratio.Divide(h_mue)
        plot_emu_mue_with_ratio(canvas=c, h_mue=h_mue, h_emu=h_emu, h_ratio=h_ratio,
                                filename=outputdir+'/sr_'+var+'_emu_over_mue_data_minus_fake',
                                label='SR: data-fake')
        fpt1_histos['data_fake']['sr'] = h_ratio

        h_emu  = buildTotBkg(histoFakeBkg=None,
                             histosSimBkgs=dict([(g.name,
                                                  g.getHistogram(variable=var, selection='sr_emu_os', cacheIt=True))
                                                 for g in simBkgs]))
        h_mue  = buildTotBkg(histoFakeBkg=None,
                             histosSimBkgs=dict([(g.name,
                                                  g.getHistogram(variable=var, selection='sr_mue_os', cacheIt=True))
                                                 for g in simBkgs]))
        h_ratio = h_emu.Clone(h_emu.GetName().replace('emu', 'emu_over_mue'))
        h_ratio.Divide(h_mue)
        plot_emu_mue_with_ratio(canvas=c, h_mue=h_mue, h_emu=h_emu, h_ratio=h_ratio,
                                filename=outputdir+'/sr_'+var+'_emu_over_mue_simbkg',
                                label='SR: simbkg')
        fpt1_histos['sim_bkg']['sr'] = h_ratio

        # continue # if you don't have vr
        # vr
        h_emu = fake.getHistogram(variable=var, selection='vr_emu_os', cacheIt=True)
        h_mue = fake.getHistogram(variable=var, selection='vr_mue_os', cacheIt=True)
        h_ratio = h_emu.Clone(h_emu.GetName().replace('emu', 'emu_over_mue'))
        h_ratio.Divide(h_mue)
        plot_emu_mue_with_ratio(canvas=c, h_mue=h_mue, h_emu=h_emu, h_ratio=h_ratio,
                                filename=outputdir+'/vr_'+var+'_emu_over_mue_fake',
                                label='VR: fake')

        h_emu = data.getHistogram(variable=var, selection='vr_emu_os', cacheIt=True)
        h_mue = data.getHistogram(variable=var, selection='vr_mue_os', cacheIt=True)
        h_ratio = h_emu.Clone(h_emu.GetName().replace('emu', 'emu_over_mue'))
        h_ratio.Divide(h_mue)
        plot_emu_mue_with_ratio(canvas=c, h_mue=h_mue, h_emu=h_emu, h_ratio=h_ratio,
                                filename=outputdir+'/vr_'+var+'_emu_over_mue_data',
                                label='VR: data')

        h_emu = data.getHistogram(variable=var, selection='vr_emu_os', cacheIt=True).Clone('vr_emu_data_minus_fake')
        h_mue = data.getHistogram(variable=var, selection='vr_mue_os', cacheIt=True).Clone('vr_mue_data_minus_fake')
        h_emu.Add(fake.getHistogram(variable=var, selection='vr_emu_os', cacheIt=True), -1.0)
        h_mue.Add(fake.getHistogram(variable=var, selection='vr_mue_os', cacheIt=True), -1.0)
        h_ratio = h_emu.Clone(h_emu.GetName().replace('emu', 'emu_over_mue'))
        h_ratio.Divide(h_mue)
        plot_emu_mue_with_ratio(canvas=c, h_mue=h_mue, h_emu=h_emu, h_ratio=h_ratio,
                                filename=outputdir+'/vr_'+var+'_emu_over_mue_data_minus_fake',
                                label='VR: data-fake')
        fpt1_histos['data_fake']['vr'] = h_ratio

        h_emu  = buildTotBkg(histoFakeBkg=None,
                             histosSimBkgs=dict([(g.name,
                                                  g.getHistogram(variable=var, selection='vr_emu_os', cacheIt=True))
                                                 for g in simBkgs]))
        h_mue  = buildTotBkg(histoFakeBkg=None,
                             histosSimBkgs=dict([(g.name,
                                                  g.getHistogram(variable=var, selection='vr_mue_os', cacheIt=True))
                                                 for g in simBkgs]))
        h_ratio = h_emu.Clone(h_emu.GetName().replace('emu', 'emu_over_mue'))
        h_ratio.Divide(h_mue)
        plot_emu_mue_with_ratio(canvas=c, h_mue=h_mue, h_emu=h_emu, h_ratio=h_ratio,
                                filename=outputdir+'/vr_'+var+'_emu_over_mue_simbkg',
                                label='VR: simbkg')
        fpt1_histos['sim_bkg']['vr'] = h_ratio

        print "chi2 test p-value between SR and VR:"
        print "data-fake : ",fpt1_histos['data_fake']['vr'].Chi2Test(fpt1_histos['data_fake']['sr'], 'WW')
        print "sim-bkg   : ",fpt1_histos['sim_bkg'  ]['vr'].Chi2Test(fpt1_histos['sim_bkg'  ]['sr'], 'WW')
        print "chi2 test chi2/NDF between SR and VR:"
        print "data-fake : ",fpt1_histos['data_fake']['vr'].Chi2Test(fpt1_histos['data_fake']['sr'], 'WW CHI2/NDF')
        print "sim-bkg   : ",fpt1_histos['sim_bkg'  ]['vr'].Chi2Test(fpt1_histos['sim_bkg'  ]['sr'], 'WW CHI2/NDF')
    return

def plot_emu_mue_with_ratio(canvas=None, h_mue=None, h_emu=None, h_ratio=None,
                            filename='', figformats=['png', 'eps'], label=''):
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
    leg = topRightLegend(c, 0.275, 0.275, shift=0.05)
    leg.SetBorderSize(0)
    leg.SetHeader(label)
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
    yax.SetTitle('e#mu / #mue ratio')
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
    for ext in figformats:
        c.SaveAs("{0}.{1}".format(filename, ext))
        print "{0}.{1}".format(filename, ext)

def regions_to_plot() : return ['sr_mue_os', 'sr_emu_os']

if __name__=='__main__':
    main()
