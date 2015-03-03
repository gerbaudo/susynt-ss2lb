#!/bin/env python

# Generic utility functions for SusyTest0
#
# davide.gerbaudo@gmail.com
# 2013-07-25

import collections
import difflib
from functools import wraps
import glob
import json
import os
import re
import sys
import subprocess
import unittest

def getCommandOutput(command):
    "lifted from supy (https://github.com/elaird/supy/blob/master/utils/io.py)"
    p = subprocess.Popen(command, shell = True, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    stdout,stderr = p.communicate()
    return {"stdout":stdout, "stderr":stderr, "returncode":p.returncode}

def filterWithRegexp(stringList, regexp, func=lambda x : x) :
    return [d for d in stringList if re.search(regexp, func(d))]

def excludeWithRegexp(stringList, regexp, func=lambda x : x) :
    return [d for d in stringList if not re.search(regexp, func(d))]

def findLatestOneOrTwoRootFiles(dir) :
    files = filter(os.path.isfile, glob.glob(dir + "*.root"))
    files.sort(key=lambda x: os.path.getmtime(x))
    return files[-2:] if len(files)>=2 else files
def findLastRootFile(dir) : return findLatestOneOrTwoRootFiles(dir)[-1]
def guessMonthDayTag(name) :
    "extract a 'Xxx_yy' tag with a 3-char month and a day"
    match = re.search('(?P<tag>'
                      '(?:Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)_\d+)',
                      name)
    if match : return match.group('tag')
def guessMonthDayTagFromLastRootFile(dir, debug) :
    lastFile = findLastRootFile(dir)
    if lastFile : return guessMonthDayTag(lastFile)
    elif debug : print "guessMonthDayTagFromLastRootFile: no root files"
def guessLatestTagFromLatestRootFiles(dir, debug) :
    """
    Latest tag if there are at least 2 root files; otherwise empty string.
    The tag does not have to be Mmm_dd.
    """
    files = [f.replace('.root','').rstrip() for f in findLatestOneOrTwoRootFiles(dir)]
    def commonSuffix(ll) : return os.path.commonprefix([l[::-1] for l in ll])[::-1]
    suffix = commonSuffix(files)
    if not len(suffix) : return ''
    def cleanupSampleName(suff) :
        "if the suffix includes part of a sample name, it can be isolated by the '.'"
        lastDot = max(0, suffix.rfind('.'))
        return suffix[suffix.find('_', lastDot):]
    if debug : print "guessLatestTagFromLatestRootFiles: cleanup '%s'"%suffix
    return cleanupSampleName(suffix) if len(suffix) else suffix
def isMonthDayTag(tag) : return guessMonthDayTag(tag)
def commonPrefix(list) : return os.path.commonprefix(list)
def commonSuffix(list) : return os.path.commonprefix([l[::-1] for l in list])[::-1]
def longestCommonSubstring(s1, s2) :
    m = difflib.SequenceMatcher(None, s1, s2).find_longest_match(0, len(s1), 0, len(s2))
    return s1[m.a : m.a+m.size]

class Memoize :
    """A class to cache cpu-intensive functions.
    Arguments must be hashable.
    See for example
    http://stackoverflow.com/questions/1988804/what-is-memoization-and-how-can-i-use-it-in-python
    """
    def __init__(self, f) :
        self.f = f
        self.memo = {}
    def __call__(self, *args) :
        if not args in self.memo : self.memo[args] = self.f(*args)
        return self.memo[args]
def enumFromHeader(filename, enumName) :
    """
    Given a c header file, extract the enum as a dict of key:values.
    From:
    https://mail.python.org/pipermail/python-list/2009-August/548422.html
    Modified to also get the enum name.
    """
    verbose = False
    file_data = open(filename).read()
    # Remove comments and preprocessor directives
    file_data = ' '.join(line.split('//')[0].split('#')[0] for line in file_data.splitlines())
    file_data = ' '.join(re.split(r'\/\*.*?\*\/', file_data))
    # Look for enums: In the first { } block after the keyword "enum"
    enums = [(text.split('{')[0].replace('enum','').strip(), text.split('{')[1].split('}')[0])
             for text in re.split(r'\benum\b', file_data)[1:]]
    enum = dict()
    for enum_name, enum_keyvals in enums:
        last_value = -1
        for key_name in enum_keyvals.split(','):
            if '=' in key_name:
                key_name, key_value = key_name.split('=')
                key_value = int(str(eval(key_value)), 0) # int(str()) to catch shift << as well
            else:
                key_value = last_value + 1
            last_value = key_value
            key_name = key_name.strip()
            if enum_name == enumName : enum[key_name] = key_value
            if verbose : print '%s = %d' % (key_name, key_value)
        if verbose : print
    return enum
def dictKeysSortedByValue(aDict={}) :
    "Given a dict, return its keys sorted by their values"
    return [x[0] for x in sorted(aDict.iteritems(), key=operator.itemgetter(1))]
def dictSum(d0, d1) :
    "see http://stackoverflow.com/questions/6005066/adding-dictionaries-together-python"
    return dict(d0, **d1)
def first(listOrDict) :
    lod = listOrDict
    return lod.itervalues().next() if type(lod) is dict else lod[0] if lod else None
def json_write(obj, fname) :
    with open(fname, 'w') as out :
        json.dump(obj, out)
def json_read(fname) :
    with open(fname) as inp :
        return json.load(inp)
def rmIfExists(filename) :
    if os.path.exists(filename) : os.remove(filename)
def mkdirIfNeeded(dirname) :
    dest_dir = None
    if os.path.exists(dirname) and os.path.isdir(dirname) :
        dest_dir = dirname
    elif not os.path.exists(dirname) :
        os.makedirs(dirname)
        dest_dir = dirname
    return dest_dir
def verticalSlice(list2d) :
    "http://stackoverflow.com/questions/6253586/python-vertical-array-slicing"
    return zip(*list2d)
def linearTransform(values, targetRange=[0.0,1.0]) :
    xLoT, xHiT = targetRange[0], targetRange[1]
    xLoO, xHiO = min(values), max(values)
    oriRange, tarRange = (xHiO-xLoO), (xHiT-xLoT)
    return [(xLoT + (x-xLoO)*tarRange/oriRange) if oriRange else 0.0
            for x in values]
def cumsum(l, leftToRight=True) :
    #return numpy.cumsum(l) # not available ?
    return [sum(l[:i]) for i in range(1,len(l)+1)] if leftToRight else [sum(l[-i:]) for i in range(1,len(l)+1)][::-1]
def mergeOuter(bc, nOuter=2) : # add over/underflow in the first/last bin
    return [sum(bc[:nOuter])] + bc[nOuter:-nOuter] + [sum(bc[-nOuter:])]
def transposeDict(d) :
    "given a dict[key1][key2] return dict[key2][key1]"
    possible_k2s = [sorted(row.keys()) for row in d.values()]
    assert len(frozenset(possible_k2s[0]))==len(possible_k2s[0]),"ambigous keys, cannot transpose %s"%str(possible_k2s[0])
    assert len(frozenset([frozenset(ks) for ks in possible_k2s])),"rows with different keys, cannot transpose %s"%str(possible_k2s)
    k2s = first(possible_k2s)
    return dict([(k2, dict([(k1, d[k1][k2]) for k1 in d.keys()])) for k2 in k2s])
def renameDictKey(d, old, new) :
    d[new] = d.pop(old)
    return d
def sortedAs(d={}, sortedKeys=[]) :
    "take a dictionary and access its item with a specified order; unspecified keys go at the end"
    allKeys = d.keys()
    keys = [k for k in sortedKeys if k in allKeys] + [k for k in allKeys if k not in sortedKeys]
    #return collections.OrderedDict([(k, d[k]) for k in keys]) # OrderedDict not available in 2.6.5 ??
    return [(k, d[k]) for k in keys]
def remove_duplicates(seq=[]) :
    "see http://stackoverflow.com/questions/480214/how-do-you-remove-duplicates-from-a-list-in-python-whilst-preserving-order"
    seen = set()
    seen_add = seen.add
    return [ x for x in seq if x not in seen and not seen_add(x)]

def rootcoredir():
    return os.environ['ROOTCOREDIR']

def add_susynt_path():
    """
    provide access to the susynt modules
    """
    rootcoredir = os.environ['ROOTCOREBIN']
    susynt_path = os.path.realpath(os.path.join(rootcoredir, '../SusyNtuple'))
    if susynt_path not in sys.path:
        sys.path.append(susynt_path)
def import_susyntutils():
    add_susynt_path()
    import susyntuple.utils as susynt
    return susynt
def import_susyntcutflow():
    add_susynt_path()
    import susyntuple.cutflow as cutflow
    return cutflow

def print_running_conditions(parser, opts):
    print "working from {0}".format(os.getcwd())
    print "being called as : {0}".format(' '.join(os.sys.argv))
    allOptions = [x.dest for x in parser._get_all_options()[1:]]
    print "options parsed:\n"+'\n'.join("%s : %s"%(o, str(getattr(opts, o))) for o in allOptions)

#
# testing
#
class testGuessMdTag(unittest.TestCase) :
    def testKnownValues(self) :
        knownValues = [('out/foo/ttbar_Aug_23.root',     'Aug_23'),
                       ('foo/baz_Aug/ttbar_Aug_23.root', 'Aug_23'),
                       ('out/foo/ttbar_Aug_2.root',      'Aug_2'),
                       ('out/foo/ttbar_August_2.root',    None),
                       ]
        for s, tag in  knownValues :
            gTag = guessMonthDayTag(s)
            self.assertEqual(tag, gTag)

if __name__ == "__main__":
    unittest.main()
