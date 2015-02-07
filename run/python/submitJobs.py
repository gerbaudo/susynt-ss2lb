#!/bin/env python

# Submit the SusyntHlfv
#
#
# davide.gerbaudo@gmail.com
# Jan 2013 -- first version (within SusyTest0)
# Jun 2014 -- rewrite (within SusyntHlfv)

import optparse
import os


import dataset
import utils

def main():
    options = parse_options()
    exe     = options.executable
    inputdf = options.input
    include = options.include_regexp
    exclude = options.exclude_regexp
    tag     = options.tag
    verbose = options.verbose
    submit  = options.submit

    datasets = dataset.build_all_datasets_from_dir_or_file(inputdf)
    datasets = utils.filterWithRegexp (datasets, include, lambda _: _.name) if include else datasets
    datasets = utils.excludeWithRegexp(datasets, exclude, lambda _: _.name) if exclude else datasets

    for dset in datasets :
        if not get_filelist(dset.name):
            print "# skipping (missing filelist) {0}".format(dset.name)
            continue
        script = get_batch_script(dset, options)
        if not script:
            if verbose : print "skipping (do-not-overwrite) {0}".format(dset.name)
            continue
        cmd = "sbatch %s" % script
        if verbose : print cmd
        if submit :
            out = utils.getCommandOutput(cmd)
            if verbose : print out['stdout']
    if not submit : print "This was a dry run; use '--submit' to actually submit the jobs"



def valid_executables():
    return ['selection', 'seltuple', 'matrix_prediction'] # todo 'seltuple', 'faketupl'...]

def parse_options():
    parser = optparse.OptionParser()
    parser.add_option('--selection',  action='store_true', default=False, help='run Selector')
    parser.add_option('--seltuple',  action='store_true', default=False, help='run Selector to make ntuple')
    parser.add_option('--matrix-prediction',  action='store_true', default=False, help='run MatrixPrediction')
    parser.add_option('-i', '--input', default='samples/', help='input directory or file (default: ./samples/)')
    parser.add_option("-o", "--do-not-overwrite", action="store_true", default=False, help="do not overwrite existing batch scripts")
    parser.add_option("-O", "--other-opt", help="other options that will be passed on to the executable; double quotes if necessary")
    parser.add_option('-q', '--queue', default='atlas_all', help="batch queue, default atlas_all")
    parser.add_option('-s', '--include-regexp', help="select only matching samples (default '.*')")
    parser.add_option('-e', '--exclude-regexp', help="exclude matching samples")
    parser.add_option("-S", "--submit", action='store_true', default=False, help="submit jobs (default dry run)")
    parser.add_option("-t", "--tag", help="batch tag")
    parser.add_option("-C", "--use-cache", action='store_true', default=False, help="create TEventList cache (or use existing one)")
    parser.add_option("-v", "--verbose", action="store_true", default=False, help="print more details about what is going on")
    (options, args) = parser.parse_args()
    options.executable = get_executable(parser, options)
    if not options.executable : parser.error("specify one executable")
    if not os.path.exists(options.input) : parser.error("invalid input '%s'"%parser.input)
    if not options.tag : parser.error("specify a tag")
    return options

def get_executable(parser, options):
    """
    the executable is a label for each main executable; valid
    labels must be listed in 'valid_executables'. We then expect to
    have exe-specific cases for 'get_template_script' and 'get_subdir'
    """
    exe_name = None
    for exe in valid_executables():
        if getattr(options, exe)==True:
            exe_name = exe
    if not exe : parser.error("specify one executable")
    return exe_name

def get_batch_dir(exe) : return mk_dest_dir('batch/'+get_subdir(exe))
def get_cache_dir(exe) :
    mk_dest_dir('out/cache/'+get_subdir(exe)) # pre-make dir; need to then pass it with unaccessible env var
    return '${SLURM_SUBMIT_DIR}/out/cache/'+get_subdir(exe)
# log and out go to separate 'tag' subdirs to keep things tidy
def get_log_dir(exe, tag='') : return mk_dest_dir('/'.join(['log', get_subdir(exe), tag]))
def get_out_dir(exe, tag='') : return mk_dest_dir('/'.join(['out', get_subdir(exe), tag]))

def get_template_script(exe_name):
    template  = ''
    template += 'batch/templates/selection.sh' if exe_name=='selection' else ''
    template += 'batch/templates/seltuple.sh' if exe_name=='seltuple' else ''
    template += 'batch/templates/matrix_prediction.sh' if exe_name=='matrix_prediction' else ''
    # other cases here
    return template

def get_subdir(exe_name) :
    subdir = None
    if exe_name=='selection' : subdir = 'selection'
    if exe_name=='seltuple' : subdir = 'selection_tuple'
    if exe_name=='matrix_prediction' : subdir = 'matrix_prediction'
    # other cases here
    return subdir

def mk_dest_dir(dir_name) :
    if not os.path.isdir(dir_name)  : os.makedirs(dir_name)
    return dir_name

def get_filelist(dataset_name, filelist_dir='filelist/'):
    path = os.path.join(filelist_dir, dataset_name+'.txt')
    return path if os.path.exists(path) else None

def get_batch_script(dset, options):
    exe     = options.executable
    tag     = options.tag
    queue   = options.queue
    filelist = get_filelist(dset.name)
    outdir   = get_out_dir(exe, options.tag)
    logdir   = get_log_dir(exe, options.tag)
    batchdir = get_batch_dir(exe)
    cachedir = get_cache_dir(exe)

    dsname = dset.name
    jobname = dsname
    script_template = get_template_script(exe)
    batch_script = batchdir+'/'+dsname+'.sh'
    out_rootfile = outdir+'/'+dsname+'.root'
    out_logfile  = logdir+'/'+dsname+'.log'
    exe_options = ''
    exe_options += " --event-list %s"%(cachedir+'/'+dsname+'.root') if options.use_cache else ''

    if options.do_not_overwrite and os.path.exists(batch_script):
        if options.verbose : print "file exists : {0}".format(batch_script)
        return None
    out_file = open(batch_script, 'w')
    for line in open(script_template).readlines() :
        # note to self: could use string.format, but this can also handle special cases (e.g. output)
        line = line.replace('%(filelist)s', filelist)
        line = line.replace('%(jobname)s', jobname)
        line = line.replace('%(queue)s', queue)
        line = line.replace('%(logfile)s', out_logfile)
        line = line.replace('%(local_outfilename)s', os.path.basename(out_rootfile))
        line = line.replace('%(outfilename)s', out_rootfile) # will need special treatment for multiple output files
        line = line.replace('%(opt)s', exe_options)
        line = line.replace('%(samplename)s', dsname)
        out_file.write(line)
    out_file.close()
    return batch_script

if __name__=='__main__':
    main()
