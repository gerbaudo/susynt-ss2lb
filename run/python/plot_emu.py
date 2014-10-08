#!/bin/env python

# make basic plots for the emu selection
#
# davide.gerbaudo@gmail.com
# Sep 2014

import collections
import datetime
import math
import optparse
import os
import pprint

import dataset
from rootUtils import (drawAtlasLabel
                       ,getBinContents
                       ,getMinMax
                       ,graphWithPoissonError
                       ,increaseAxisFont
                       ,topRightLegend
                       ,importRoot
                       ,integralAndError
                       ,setAtlasStyle
                       ,writeObjectsToFile
                       )
r = importRoot()
from utils import (first
                   ,getCommandOutput
                   ,mkdirIfNeeded
                   ,filterWithRegexp
                   ,remove_duplicates
                   ,sortedAs
                   )

import utils
susyntutils = utils.import_susyntutils()
r = susyntutils.import_root()
susyntutils.load_packages()

from CutflowTable import CutflowTable

import systUtils # from Susy2014_Nt_dev/SusyTest0/run/python/systUtils.py; needed?

usage="""
This code is used either (1) to fill the histos, or (2) to make plots
and tables. The output of (1) is used as input of (2).

Required inputs: trees produced with `submitJobs.py --seltuple` and
with `--matrix-prediction`.

Example usage ('fill' mode):
%prog \\
 --input-other  out/selection_tuple/Jul_26 \\
 --input-fake out/matrix_prediction/Jul_26/ \\
 --output-dir out/plot_emu/Jul_26/histos \\
 --verbose \\
 2>&1 | tee log/plot_emu/Jul_26/fill.log

Example usage ('plot' mode):
%prog \\
 --input-dir out/plot_emu/Jul_26/histos \\
 --output-dir out/plot_emu/Jul_26/ \\
 --verbose \\
 2>&1 | tee log/plot_emu/Jul_26/plot.log"""

def main() :
    parser = optparse.OptionParser(usage=usage)
    parser.add_option('-b', '--batch',  action='store_true', help='submit to batch (used in fill mode)')
    parser.add_option('-g', '--group', help='group to be processed (used only in fill mode)')
    parser.add_option('-f', '--input-fake', help='location of fake trees')
    parser.add_option('-O', '--input-other', help='location other trees')
    parser.add_option('-i', '--input-dir')
    parser.add_option('-o', '--output-dir')
    parser.add_option('--samples-dir', default='samples/',
                      help='directory with the list of samples; default ./samples/')
    parser.add_option('-v', '--verbose', action='store_true', default=False)
    parser.add_option('--debug', action='store_true', default=False)

    (opts, args) = parser.parse_args()
    inOtherSpecified, inDirSpecified = opts.input_other!=None, opts.input_dir!=None
    eitherMode = inOtherSpecified != inDirSpecified
    if not eitherMode : parser.error("Run either in 'fill' or 'plot' mode")
    mode = 'fill' if inOtherSpecified else 'plot' if inDirSpecified else None
    requiredOptions = (['input_fake', 'input_other', 'output_dir'] if mode=='fill'
                       else ['input_dir', 'output_dir'])
    allOptions = [x.dest for x in parser._get_all_options()[1:]]
    def optIsNotSpecified(o) : return not hasattr(opts, o) or getattr(opts,o) is None
    if any(optIsNotSpecified(o) for o in requiredOptions):
        parser.error('Missing required option\n'
                     +'\n'.join(["%s : %s"%(o, getattr(opts, o)) for o in requiredOptions]))
    if opts.verbose:
        print ('\nUsing the following options:\n'
               +'\n'.join("%s : %s"%(o, str(getattr(opts, o))) for o in allOptions))

    if opts.debug : dataset.Dataset.verbose_parsing = True
    groups = dataset.DatasetGroup.build_groups_from_files_in_dir(opts.samples_dir)
    groups.append(first([g for g in groups if g.is_data]).clone_data_as_fake())
    if opts.group : groups = [g for g in groups if g.name==opts.group]
    print '\n'.join("group {0} : {1} samples".format(g.name, len(g.datasets)) for g in groups)

    if   mode=='fill' : runFill(opts, groups)
    elif mode=='plot' : runPlot(opts, groups)

def runFill(opts, groups) :
    batchMode    = opts.batch
    inputFakeDir = opts.input_fake
    inputGenDir  = opts.input_other
    outputDir    = opts.output_dir
    verbose      = opts.verbose

    if verbose : print "filling histos"
    mkdirIfNeeded(outputDir)
    systematics = ['NOM']
    if verbose : print "about to loop over these systematics:\n %s"%str(systematics)
    if batchMode:
        for syst in systematics: # todo : probably not needed, leave it in for now
            for group in groups :
                submit_batch_fill_job_per_group(group, opts)
    else:
        selections = regions_to_plot()
        variables = variables_to_plot()
        group_names = [g.name for g in groups]
        counters_all_groups = bookCounters(group_names, selections)
        histos_all_groups = bookHistos(variables, group_names, selections)
        for group in groups:
            tree_name = 'hlfv_tuple'
            chain = r.TChain(tree_name)
            input_dir = opts.input_fake if group.name=='fake' else opts.input_other
            for dataset in group.datasets:
                chain.Add(os.path.join(input_dir, dataset.name+'.root'))
            if opts.verbose:
                print "{0} : {1} entries from {2} samples".format(group.name,
                                                                  chain.GetEntries(),
                                                                  len(group.datasets))
            selection = 'emu'
            histos = histos_all_groups[group.name][selection]
            counters = counters_all_groups[group.name][selection]
            for iEntry, event in enumerate(chain):
                run_num = event.pars.runNumber
                evt_num = event.pars.eventNumber
                weight =  event.pars.weight
                l0 = addTlv(event.l0)
                l1 = addTlv(event.l1)
                isEl0, isMu0 = l0.isEl, l0.isMu
                isEl1, isMu1 = l1.isEl, l1.isMu
                isEmu = int((isEl0 and isMu1) or (isMu1 and isMu0))
                isSameSign = int((l0.charge * l1.charge)>0)
                if l0.p4.Pt()<45.0 or not isSameSign : continue
                histos['onebin'].Fill(1.0, weight)
                histos['pt0'].Fill(l0.p4.Pt(), weight)
                histos['pt1'].Fill(l1.p4.Pt(), weight)
                counters += (weight) # if passSels[sel] else 0.0)

            for v in ['onebin', 'pt0', 'pt1']:
                h = histos[v]
                print "{0}: integral {1}, entries {2}".format(h.GetName(), h.Integral(), h.GetEntries())
        plotting_groups = dict([(g.name, Group(g.name)) for g in groups])
        saveHistos(plotting_groups, histos_all_groups, outputDir, opts.verbose)
        # print counters

def runPlot(opts, groups) :
    inputDir     = opts.input_dir
    outputDir    = opts.output_dir
    sysOption    = 'NOM' #opts.syst
    excludedSyst = None #opts.exclude
    verbose      = opts.verbose
    mkdirIfNeeded(outputDir)
    buildTotBkg = systUtils.buildTotBackgroundHisto
    buildStat = systUtils.buildStatisticalErrorBand
    buildSyst = systUtils.buildSystematicErrorBand

    selections = regions_to_plot()
    variables = variables_to_plot()
    plot_groups = [Group(g.name) for g in groups]
    for group in plot_groups :
        group.setHistosDir(inputDir)
        group.exploreAvailableSystematics(verbose)
        group.filterAndDropSystematics(sysOption, excludedSyst, verbose)

    mkdirIfNeeded(outputDir)
    systematics = ['NOM']
    anySys = sysOption==None
    if sysOption=='fake'   or anySys : systematics += systUtils.fakeSystVariations()
    if sysOption=='object' or anySys : systematics += systUtils.mcObjectVariations()
    if sysOption=='weight' or anySys : systematics += systUtils.mcWeightVariations()
    if sysOption and sysOption.count(','):
        systematics = [s for s in systUtils.getAllVariations() if s in sysOption.split(',')]
    elif sysOption in systUtils.getAllVariations() : systematics = [sysOption]
    if not anySys and len(systematics)==1 and sysOption!='NOM' : raise ValueError("Invalid syst %s"%str(sysOption))
    if excludedSyst : systematics = [s for s in systematics if s not in filterWithRegexp(systematics, excludedSyst)]
    if verbose : print "using the following systematics : %s"%str(systematics)

    fakeSystematics = [s for s in systematics if s in systUtils.fakeSystVariations()]
    mcSystematics = [s for s in systematics if s in systUtils.mcObjectVariations() + systUtils.mcWeightVariations()]

    simBkgs = [g for g in plot_groups if g.isMcBkg]
    data = findByName(plot_groups, 'data')
    fake = findByName(plot_groups, 'fake')
    signal = findByName(plot_groups, 'signal')
    print 'names_stacked_groups to be improved'
    names_stacked_groups = [g.name for g in groups if g.name not in ['data', 'signal']]
    for sel in selections :
        if verbose : print '-- plotting ',sel
        for var in variables :
            if verbose : print '---- plotting ',var
            for g in plot_groups : g.setSystNominal()
            nominalHistoData    = data.getHistogram(variable=var, selection=sel, cacheIt=True)
            nominalHistoSign    = signal.getHistogram(variable=var, selection=sel, cacheIt=True)
            nominalHistoFakeBkg = fake.getHistogram(variable=var, selection=sel, cacheIt=True)
            nominalHistosSimBkg = dict([(g.name, g.getHistogram(variable=var, selection=sel, cacheIt=True))
                                        for g in simBkgs])
            nominalHistosBkg    = dict([('fake', nominalHistoFakeBkg)] +
                                       [(g, h) for g, h in nominalHistosSimBkg.iteritems()])
            nominalHistoTotBkg  = buildTotBkg(histoFakeBkg=nominalHistoFakeBkg,
                                              histosSimBkgs=nominalHistosSimBkg)
            statErrBand = buildStat(nominalHistoTotBkg)
            systErrBand = buildSyst(fake=fake, simBkgs=simBkgs, variable=var, selection=sel,
                                    fakeVariations=fakeSystematics, mcVariations=mcSystematics,
                                    verbose=verbose)

            plotHistos(histoData=nominalHistoData, histoSignal=nominalHistoSign,
                       histoTotBkg=nominalHistoTotBkg, histosBkg=nominalHistosBkg,
                       statErrBand=statErrBand, systErrBand=systErrBand,
                       stack_order=names_stacked_groups,
                       canvasName=(sel+'_'+var), outdir=outputDir, verbose=verbose)
    for group in plot_groups :
        group.printVariationsSummary()

def submit_batch_fill_job_per_group(group, opts):
    verbose = opts.verbose

    group_name = group.name if hasattr(group, 'name') else group
    newOptions  = " --input-other %s" % opts.input_other
    newOptions += " --input-fake %s" % opts.input_fake
    newOptions += " --output-dir %s" % opts.output_dir
    newOptions += " --group %s" % group_name
    newOptions += (" --verbose " if opts.verbose else '')
    template = 'batch/templates/plot_emu.sh'
    log_dir = mkdirIfNeeded('log/plot_emu')
    script_dir = mkdirIfNeeded('batch/plot_emu')
    print "trying to join '{0}' with '{1}'".format(script_dir, group_name+'.sh')
    script_name = os.path.join(script_dir, group_name+'.sh')
    script_file = open(script_name, 'w')
    script_file.write(open(template).read()
                      .replace('%(opt)s', newOptions)
                      .replace('%(logfile)s', log_dir+'/'+group_name+'.log')
                      .replace('%(jobname)s', group_name))
    script_file.close()
    cmd = "sbatch %s"%script_name
    if verbose : print cmd
    out = getCommandOutput(cmd)
    if verbose : print out['stdout']
    if out['stderr'] : print  out['stderr']

def countAndFillHistos(samplesPerGroup={}, syst='', verbose=False, outdir='./') :

    selections = allRegions()
    variables = variablesToPlot()

    mcGroups, fakeGroups = mcDatasetids().keys(), ['fake']
    objVariations, weightVariations, fakeVariations = systUtils.mcObjectVariations(), systUtils.mcWeightVariations(), systUtils.fakeSystVariations()
    def groupIsRelevantForSys(g, s) :
        isRelevant = (s=='NOM' or (g in mcGroups and s in objVariations+weightVariations) or (g in fakeGroups and s in fakeVariations))
        if verbose and not isRelevant : print "skipping %s for %s"%(g, s)
        return isRelevant
    def dropIrrelevantGroupsForThisSys(groups, sys) : return dict((g, samples) for g, samples in groups.iteritems() if groupIsRelevantForSys(g, syst))
    def dropSamplesWithoutTree(samples) : return [s for s in samples if s.hasInputHftTree(msg='Warning! ')]
    def dropGroupsWithoutSamples(groups) : return dict((g, samples) for g, samples in groups.iteritems() if len(samples))
    samplesPerGroup = dropIrrelevantGroupsForThisSys(samplesPerGroup, syst)
    samplesPerGroup = dict((g, dropSamplesWithoutTree(samples)) for g, samples in samplesPerGroup.iteritems())
    samplesPerGroup = dropGroupsWithoutSamples(samplesPerGroup)

    groups = samplesPerGroup.keys()
    counters = bookCounters(groups, selections)
    histos = bookHistos(variables, groups, selections)
    for group, samplesGroup in samplesPerGroup.iteritems() :
        logLine = "---->"
        if verbose : print 1*' ',group
        histosGroup = histos  [group]
        countsGroup = counters[group]
        for sample in samplesGroup :
            if verbose : logLine +=" %s"%sample.name
            fillAndCount(histosGroup, countsGroup, sample, blind=False)
        if verbose : print logLine
    if verbose : print 'done'
    return counters, histos

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
#___________________________________________________________
class BaseSampleGroup(object) :
    def __init__(self, name) :
        self.name = name
        self.setSystNominal()
        self.varCounts = collections.defaultdict(dict)
    @property
    def label(self) : return self.groupname if hasattr(self, 'groupname') else self.name
    @property
    def isFake(self) : return self.label=='fake'
    @property
    def isData(self) : return self.label=='data'
    @property
    def isMc(self) : return not (self.isFake or self.isData)
    @property
    def isSignal(self) : return self.label=='signal'
    @property
    def isMcBkg(self) : return self.isMc and not self.isSignal
    def isNeededForSys(self, sys) :
        return (sys=='NOM'
                or (self.isMc and sys in systUtils.mcWeightVariations())
                or (self.isMc and sys in systUtils.mcObjectVariations())
                or (self.isFake and sys in systUtils.fakeSystVariations()))
    def setSystNominal(self) : return self.setSyst()
    def setSyst(self, sys='NOM') :
        nominal = 'NOM' # do we have differnt names for nom (mc vs fake)?
        self.isObjSys    = sys in systUtils.mcObjectVariations()
        self.isWeightSys = sys in systUtils.mcWeightVariations()
        self.isFakeSys   = sys in systUtils.fakeSystVariations()
        def nameObjectSys(s) : return s if self.isMc else nominal
        def nameWeightSys(s) : return s if self.isMc else nominal
        def nameFakeSys(s) : return s if self.isFake else nominal
        def identity(s) : return s
        sysNameFunc = nameObjectSys if self.isObjSys else nameWeightSys if self.isWeightSys else nameFakeSys if self.isFakeSys else identity
        self.syst = sysNameFunc(sys)
        self.syst = sys
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
                                    for sys, sysCount in sortedAs(sysCounts, systUtils.getAllVariations())]
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


def findByName(bsgs=[], name='') : return [b for b in bsgs if b.name==name][0]
#___________________________________________________________
class Sample(BaseSampleGroup) :
    def __init__(self, name, groupname) :
        super(Sample, self).__init__(name) # this is either the name (for data and fake) or the dsid (for mc)
        self.groupname = groupname
        self.setHftInputDir()
    def setHftInputDir(self, dir='') :
        useDefaults = not dir
        defaultDir = 'out/fakepred' if self.isFake else 'out/susyplot'
        self.hftInputDir = defaultDir if useDefaults else dir
        return self
    @property
    def weightLeafname(self) :
        leafname = 'eventweight'
        if  self.isWeightSys : leafname += " * %s"%systUtils.mcWeightBranchname(self.syst)
        return leafname
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
    def __init__(self, name) :
        super(Group, self).__init__(name)
        self.setSyst()
        self.setHistosDir()
        self._histoCache = collections.defaultdict(dict) # [syst][histoname]
    def setHistosDir(self, dir='') :
        self.histosDir = dir if dir else 'out/hft'
        return self
    @property
    def filenameHisto(self) :
        "file containig the histograms for the current syst"
        fname = "%(dir)s/%(sys)s_%(group)s.root" % {'group':self.name, 'dir':self.histosDir, 'sys':self.syst}
        # do we still need the lines below?
        print 'filenameHisto, cleanup'
        return fname
        def dataFilename(group, dir, sys) : return "%(dir)s/%(sys)s_%(gr)s.PhysCont.root" % {'dir':dir, 'gr':group, 'sys':sys}
        def fakeFilename(group, dir, sys) : return "%(dir)s/%(sys)s_fake.%(gr)s.PhysCont.root" % {'dir':dir, 'gr':group, 'sys':sys}
        def mcFilename  (group, dir, sys) : return "%(dir)s/%(sys)s_%(gr)s.root" % {'dir':dir, 'sys':sys, 'gr':group}
        fnameFunc = dataFilename if self.isData else fakeFilename if self.isFake else mcFilename
        return fnameFunc(self.name, self.histosDir, self.syst)
    def exploreAvailableSystematics(self, verbose=False) :
        systs = ['NOM']
        if self.isFake :
            systs += systUtils.fakeSystVariations()
        elif self.isMc :
            systs += systUtils.mcObjectVariations()
            systs += systUtils.mcWeightVariations()
        self.systematics = []
        for sys in systs :
            self.setSyst(sys)
            if os.path.exists(self.filenameHisto) :
                self.systematics.append(sys)
        if verbose : print "%s : found %d variations : %s"%(self.name, len(self.systematics), str(self.systematics))
    def filterAndDropSystematics(self, include='.*', exclude=None, verbose=False) :
        nBefore = len(self.systematics)
        anyFilter = include or exclude
        toBeExcluded = filter(self,systematics, exclude) if exclude else []
        systs = ['NOM'] if 'NOM' in self.systematics else []
        if include : systs += filterWithRegexp(self.systematics, include)
        if exclude : systs  = [s for s in systs if toBeExcluded and s not in toBeExcluded]
        self.systematics = systs if anyFilter else self.systematics
        self.systematics = remove_duplicates(self.systematics)
        nAfter = len(self.systematics)
        if verbose : print "%s : dropped %d systematics, left with %s"%(self.name, nBefore-nAfter, str(self.systematics))
        assert self.systematics.count('NOM')==1 or not nBefore, "%s : 'NOM' required %s"%(self.name, str(self.systematics))
    def getHistogram(self, variable, selection, cacheIt=False) :
        hname = histoName(sample=self.name, selection=selection, variable=variable)
        histo = None
        try :
            histo = self._histoCache[self.syst][hname]
        except KeyError :
            file = r.TFile.Open(self.filenameHisto)
            if not file : print "missing file %s"%self.filenameHisto
            hname = histoName(sample=self.name, selection=selection, variable=variable)
            histo = file.Get(hname)
            if not histo : print "%s : cannot get histo %s"%(self.name, hname)
            elif cacheIt :
                histo.SetDirectory(0)
                self._histoCache[self.syst][hname] = histo
            else :
                histo.SetDirectory(0)
                file.Close()
        if variable=='onebin' and histo : self.logVariation(self.syst, selection, histo.Integral(0, -1))
        return histo
    def getBinContents(self, variable, selection) :
        return getBinContents(self.getHistogram)
#___________________________________________________________
def allGroups(noData=False, noSignal=True) :
    return ([k for k in mcDatasetids().keys() if k!='signal' or not noSignal]
            + ([] if noData else ['data'])
            + ['fake']
            )

def selectionFormulas(sel) :
    ee, em, mm = 'isEE', 'isEMU', 'isMUMU'
    pt32  = '(lept1Pt>30000.0 && lept2Pt>20000.0)'
    pt33  = '(lept1Pt>30000.0 && lept2Pt>30000.0)'
    j1    = '(L2nCentralLightJets==1)'
    j23   = '(L2nCentralLightJets==2 || L2nCentralLightJets==3)'
    vetoZ = '(L2Mll<(91200.0-10000.0) || L2Mll>(91200.0+10000.))'
    dEll  = '(TMath::Abs(deltaEtaLl)<1.5)'
    ss    = '(!isOS || L2qFlipWeight!=1.0)' # ssOrQflip
    mlj1  = 'mlj < 90000.0'
    mlj2  = 'mljj<120000.0'
    formulas = {
        'eeSR1jet'   : '('+ee+' && '+ss+' && '+j1 +' && '+pt32+' && '+vetoZ+' && '+mlj1+' && L2METrel>55000.0 && Ht>200000.0)',
        'eeSR23jets' : '('+ee+' && '+ss+' && '+j23+' && '+pt32+' && '+vetoZ+' && '+mlj2+' && L2METrel>30000.0 &&                mtmax>110000.0)',
        'mmSR1jet'   : '('+mm+' && '+ss+' && '+j1 +' && '+pt32+' && '+dEll +' && '+mlj1+' &&                     Ht>200000.0 && mtmax>110000.0)',
        'mmSR23jets' : '('+mm+' && '+ss+' && '+j23+' && '+pt33+' && '+dEll +' && '+mlj2+' &&                     Ht>200000.0)',
        'emSR1jet'   : '('+em+' && '+ss+' && '+j1 +' && '+pt33+' && '+dEll +' && '+mlj1+' &&                     Ht>200000.0 && mtmax>110000.0)',
        'emSR23jets' : '('+em+' && '+ss+' && '+j23+' && '+pt33+' && '+dEll +' && '+mlj2+' &&                     Ht>200000.0 && mtmax>110000.0)',
        }
    for f in formulas.keys() :
        formulas['pre'+f] = formulas[f].replace(mlj1, '1').replace(mlj2, '1')
    mlj1Not, mlj2Not = mlj1.replace('<','>'), mlj2.replace('<','>')
    for f in formulas.keys() :
        formulas['bld'+f] = formulas[f].replace(mlj1, mlj1Not).replace(mlj2, mlj2Not)
    return formulas[sel]

def fillAndCount(histos, counters, sample, blind=True) :
    group    = sample.group
    filename = sample.filenameHftTree
    treename = sample.hftTreename
    file = r.TFile.Open(filename)
    tree = file.Get(treename)
    selections = allRegions()
    selWeights = dict((s, r.TTreeFormula(s, selectionFormulas(s), tree)) for s in selections)
    weightFormula = r.TTreeFormula('weightFormula', sample.weightLeafname, tree)
    l1 = r.TLorentzVector()
    l2 = r.TLorentzVector()
    met = r.TLorentzVector()
    for iEvent, event in enumerate(tree) :
        weight = weightFormula.EvalInstance()
        passSels = dict((s, selWeights[s].EvalInstance()) for s in selections)
        for sel in selections : counters[sel] += (weight if passSels[sel] else 0.0)
        for sel in selections :
            fillHisto = passSels[sel]
            if blind and sample.isData :
                if sel in signalRegions() : fillHisto = False
                else : fillHisto = passSels[blindRegionFromAnyRegion(sel)] and not passSels[signalRegionFromAnyRegion(sel)]
            oneJet = event.L2nCentralLightJets==1
            mev2gev = 1.0e-3
            mljj = mev2gev*(event.mlj if oneJet else event.mljj)
            l1.SetPtEtaPhiM(event.lept1Pt*mev2gev, event.lept1Eta, event.lept1Phi, 0.0) # massless here is good enough
            l2.SetPtEtaPhiM(event.lept2Pt*mev2gev, event.lept2Eta, event.lept2Phi, 0.0)
            met.SetPtEtaPhiM(event.met*mev2gev,               0.0, event.metPhi,   0.0)
            ll = l1+l2
            ptll = ll.Pt()
            mll = ll.M()
            l1IsMu = event.lept1Flav==1
            l2IsMu = event.lept2Flav==1
            dphil0met = abs(l1.DeltaPhi(met)) if l1.Pt()>l2.Pt() else abs(l2.DeltaPhi(met))
            if fillHisto :
                histos[sel]['mll'   ].Fill(mll, weight)
                histos[sel]['mljj'  ].Fill(mljj, weight)
                histos[sel]['ptll'  ].Fill(ptll, weight)
                histos[sel]['onebin'].Fill(1.0,  weight)
                histos[sel]['dphil0met'].Fill(dphil0met, weight)
                if l1IsMu or l2IsMu:
                    dphimumet = abs(l1.DeltaPhi(met)) if l1IsMu else abs(l2.DeltaPhi(met))
                    histos[sel]['dphimumet'].Fill(dphimumet, weight)
            # checks
            if (True and fillHisto
                and sel in signalRegions()
                and (sample.isData or sample.isFake)) :
                channel = 'ee' if event.isEE else 'mm' if event.isMUMU else 'em'
                dataOrFake = 'data' if sample.isData else 'fake' if sample.isFake else 'other'
                print "ev %d run %d channel %s sel %s sample %s weight %f"%(event.runNumber, event.eventNumber, channel, sel, dataOrFake, weight)
    file.Close()

def dataSampleNames() :
    return ["period%(period)s.physics_%(stream)s"%{'period':p, 'stream':s}
            for p in ['A','B','C','D','E','G','H','I','J','L']
            for s in ['Egamma','Muons']]
def mcDatasetids() :
    print 'todo: now it is an attribute Dataset.dsid parsed from the dataset name'

def allSamplesAllGroups() :
    asg = dict( [(group, [Sample(groupname=group, name=dsid) for dsid in dsids]) for group, dsids in mcDatasetids().iteritems()]
               +[('data', [Sample(groupname='data', name=s) for s in dataSampleNames()])]
               +[('fake', [Sample(groupname='fake', name=s) for s in dataSampleNames()])])
    return asg
def allGroups() :
    return [Group(g) for g in mcDatasetids().keys()+['data']+['fake']]

def stackedGroups(groups) :
    return [g for g in allSamplesAllGroups().keys() if g not in ['data', 'signal']]

def variablesToPlot() :
    return ['onebin','mljj', 'ptll', 'mll', 'dphil0met', 'dphimumet']
    return ['pt0','pt1','mll','mtmin','mtmax','mtllmet','ht','metrel','dphill','detall',
            'mt2j','mljj','dphijj','detajj']
def histoName(sample, selection, variable) : return "h_%s_%s_%s"%(variable, sample, selection)
def bookHistos(variables, samples, selections) :
    "book a dict of histograms with keys [sample][selection][var]"
    def histo(variable, sam, sel) :
        twopi = +2.0*math.pi
        mljjLab = 'm_{lj}' if '1j' in sel else 'm_{ljj}'
        h = None
        if   v=='onebin'  : h = r.TH1F(histoName(sam, sel, 'onebin' ), ';; entries',                             1, 0.5,   1.5)
        elif v=='pt0'     : h = r.TH1F(histoName(sam, sel, 'pt0'    ), ';p_{T,l0} [GeV]; entries/bin',          12, 0.0, 240.0)
        elif v=='pt1'     : h = r.TH1F(histoName(sam, sel, 'pt1'    ), ';p_{T,l1} [GeV]; entries/bin',          12, 0.0, 240.0)
        elif v=='mll'     : h = r.TH1F(histoName(sam, sel, 'mll'    ), ';m_{l0,l1} [GeV]; entries/bin',         12, 0.0, 240.0)
        elif v=='ptll'    : h = r.TH1F(histoName(sam, sel, 'ptll'   ), ';p_{T,l0+l1} [GeV]; entries/bin',       12, 0.0, 240.0)
        elif v=='dphil0met': h= r.TH1F(histoName(sam, sel, 'dphil0met'),';#Delta#phi(l0, met) [rad]; entries/bin',  10, 0.0, twopi)
        elif v=='dphil1met': h= r.TH1F(histoName(sam, sel, 'dphil1met'),';#Delta#phi(l1, met) [rad]; entries/bin',  10, 0.0, twopi)
        else : print "unknown variable %s"%v
        h.Sumw2()
        h.SetDirectory(0)
        return h
    return dict([(sam, dict([(sel, dict([(v, histo(v, sam, sel)) for v in variables]))
                         for sel in selections]))
                 for sam in samples])
def bookCounters(samples, selections) :
    "book a dict of counters with keys [sample][selection]"
    return dict((s, dict((sel, 0.0) for sel in selections)) for s in samples)
def countTotalBkg(counters={'sample' : {'sel':0.0}}) :
    backgrounds = [g for g in counters.keys() if g!='signal' and g!='data']
    selections = first(counters).keys()
    counters['totBkg'] = dict((s, sum(counters[b][s] for b in backgrounds)) for s in selections)
def getGroupColor(g) :
    oldColors = [('data', r.kBlack), ('diboson',r.kSpring+2), ('higgs',r.kAzure-4),
                 ('signal',r.kMagenta), ('top', r.kRed+1), ('zjets', r.kOrange-2),
                 ('fake',r.kGray)]
    newColors = [] #[('signal',r.kMagenta), ('WW',r.kAzure-9), ('Higgs',r.kYellow-9)]
    colors = dict((g,c) for g,c in  oldColors + newColors)
    return colors[g]

def regions_to_plot():
    return ['emu']
def variables_to_plot():
    return ['onebin', 'pt0', 'pt1']

def plotHistos(histoData=None, histoSignal=None, histoTotBkg=None, histosBkg={},
               statErrBand=None, systErrBand=None, # these are TGraphAsymmErrors
               canvasName='canvas', outdir='./', verbose=False,
               stack_order=[],
               drawStatErr=False, drawSystErr=False,
               drawYieldAndError=False) :
    "Note: blinding can be required for only a subrange of the histo, so it is taken care of when filling"
    setAtlasStyle()
    padMaster = histoData
    if verbose : print "plotting ",padMaster.GetName()
    can = r.TCanvas(canvasName, padMaster.GetTitle(), 800, 600)
    can.cd()
    can._hists = [padMaster]
    padMaster.Draw('axis')
    can.Update() # necessary to fool root's dumb object ownership of the stack
    stack = r.THStack('stack_'+padMaster.GetName(), '')
    r.SetOwnership(stack, False)
    can._hists.append(stack)
    leg = topRightLegend(can, 0.225, 0.325)
    can._leg = leg
    leg.SetBorderSize(0)
    leg._reversedEntries = []
    for group, histo in sortedAs(histosBkg, stack_order) :
        histo.SetFillColor(getGroupColor(group))
        histo.SetLineWidth(2)
        histo.SetLineColor(r.kBlack)
        stack.Add(histo)
        can._hists.append(histo)
        leg._reversedEntries.append((histo, group, 'F'))
    for h, g, o in leg._reversedEntries[::-1] : leg.AddEntry(h, g, o) # stack goes b-t, legend goes t-b
    stack.Draw('hist same')
    histoData.SetMarkerStyle(r.kFullCircle)
    histoData.SetLineWidth(2)
    dataGraph = graphWithPoissonError(histoData)
    dataGraph.Draw('same p')
    if histoSignal :
        histoSignal.SetLineColor(getGroupColor('signal'))
        histoSignal.SetLineWidth(2)
        histoSignal.Draw('histo same')
        leg.AddEntry(histoSignal, '(m_{C1},m_{N1})=(130, 0)GeV', 'l')
    if statErrBand and drawStatErr :
        statErrBand.SetFillStyle(3006)
        statErrBand.Draw('E2 same')
        leg.AddEntry(statErrBand, 'stat', 'f')
    if systErrBand and drawSystErr :
        systErrBand.SetFillStyle(3007)
        systErrBand.Draw('E2 same')
        leg.AddEntry(systErrBand, 'syst', 'f')
    totErrBand = systUtils.addErrorBandsInQuadrature(statErrBand, systErrBand)
    if totErrBand :
        totErrBand.Draw('E2 same')
        totErrBand.SetFillStyle(3005)
        leg.AddEntry(totErrBand, 'stat+syst', 'f')
    leg.Draw('same')
    can.Update()
    tex = r.TLatex()
    tex.SetTextSize(0.5 * tex.GetTextSize())
    tex.SetNDC(True)
    label  = "%s tot bkg : "%(can.GetName())
    label += "%.3f #pm %.3f (stat)"%(integralAndError(histoTotBkg))
    if systErrBand :
        sysUp, sysDo = systUtils.totalUpDownVariation(systErrBand)
        label += "#pm #splitline{%.3f}{%.3f} (syst)"%(sysUp, sysDo)
    if drawYieldAndError :
        tex.DrawLatex(0.10, 0.95, label)
        can.SetTopMargin(2.0*can.GetTopMargin())
    drawAtlasLabel(can, xpos=0.125, align=13)
    yMin, yMax = getMinMax([histoData, dataGraph, histoTotBkg, histoSignal, totErrBand])
    padMaster.SetMinimum(0.0)
    padMaster.SetMaximum(1.1 * yMax)
    increaseAxisFont(padMaster.GetXaxis())
    increaseAxisFont(padMaster.GetYaxis())
    can.RedrawAxis()
    can.Update() # force stack to create padMaster
    for ext in ['png','eps'] : can.SaveAs(outdir+'/'+can.GetName()+'.'+ext)

def saveHistos(samplesPerGroup={}, histosPerGroup={}, outdir='./', verbose=False) :
    for groupname, histosThisGroup in histosPerGroup.iteritems() :
        group = samplesPerGroup[groupname].setHistosDir(outdir)
        outFilename = group.filenameHisto
        writeObjectsToFile(outFilename, histosThisGroup, verbose)

# this is duplicated with plot_fake_weight_correlation.py; put it in smth like tuple_utils
tlv = r.TLorentzVector
def FourMom2TLorentzVector(fm) :
    l = tlv()
    l.SetPxPyPzE(fm.px, fm.py, fm.pz, fm.E)
    return l
def addTlv(l) :
    if not hasattr(l, 'p4') : l.p4 = FourMom2TLorentzVector(l)
    return l


if __name__=='__main__' :
    main()
