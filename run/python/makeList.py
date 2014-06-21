#!/bin/env python

# Make lists of input root files for given lists of datasets
#
# davide.gerbaudo@gmail.com
# Jan 2013 -- first version (within SusyTest0)
# Jun 2014 -- rewrite (within SusyntHlfv)

import optparse
import os
import sys

import dataset
import utils

def main():
    options = parse_options()
    inputdf = options.input
    outdir  = options.output_dir
    regexp  = options.sample_regexp
    tag     = options.tag
    verbose = options.verbose
    debug   = options.debug

    if debug : dataset.Dataset.verbose_parsing = True
    datasets = (dataset.Dataset.parse_files_in_dir(inputdf) if os.path.isdir(inputdf) else
                dataset.Dataset.parse_datasets_from_file(inputdf))
    datasets = utils.filterWithRegexp(datasets, regexp, lambda _: _.name)
    counter = {'fail':0, 'pass':0}
    for d in datasets:
        outcome = 'pass' if  d.build_filelist(gpatlas_dir(d, tag), './filelist/', verbose) else 'fail'
        counter[outcome] += 1
    if verbose:
        print "created %d filelists (%d failures)" % (counter['pass'], counter['fail'])

def gpatlas_dir(dataset, tag):
    base_dir = '/gdata/atlas/ucintprod/SusyNt'
    # note to self: if we use susy samples, you need to add 'susy_<tag>'
    prefix = 'mc12_' if dataset.is_simulation else 'data12_' if dataset.is_data else None
    if not prefix : print "dataset must be either simulation or data '%s'"%dataset.name
    return base_dir+'/'+prefix+tag

def parse_options():
    usage = """usage: %prog [options]
    Example:
    %prog -t n0145 -s 'Z(ee|mumu|tautau)' -v
    """
    parser = optparse.OptionParser(usage=usage)
    parser.add_option('-i', '--input', default='samples/', help='input directory or file (default: ./samples/)')
    parser.add_option('-o', '--output-dir', default='filelist/', help='output directory')
    parser.add_option('-s', '--sample-regexp', default='.*', help="create filelists only for matching samples (default '.*')")
    parser.add_option('-t', '--tag', help='SusyNt production tag')
    parser.add_option('-v', '--verbose', action='store_true', default=False, help='print more details about what is going on')
    parser.add_option('-d', '--debug', action='store_true', default=False, help='print even more details, only useful to debug problems')
    (options, args) = parser.parse_args()
    if not options.tag : parser.error('tag is a required option')
    if not os.path.exists(options.input) : parser.error("invalid input '%s'"%parser.input)
    return options

if __name__=='__main__':
    main()
