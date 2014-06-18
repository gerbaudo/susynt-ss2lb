#!/bin/bash

# Template to execute `run_Selector` on the greenplanet queue
#
# The parameters specified as '%(parameter)s' should be filled in by
# the submission script (python/submitJobs.py).
# The bash variables '${VAR}' are either environmental variables (or
# they are built from the available job info).
# For more details on the slurm syntax, see the slurm documentation
# (http://slurm.schedmd.com/)
#
# davide.gerbaudo@gmail.com
# June 2014

#SBATCH -p atlas_all
#SBATCH --mem-per-cpu=1gb
#SBATCH --time=24:00
#SBATCH -o %(logfile)s
#SBATCH --job-name=%(jobname)s

SCRATCH=/scratch/${USER}/${SLURM_JOB_ID}

echo "Starting on `hostname`, `date`"
mkdir -p ${SCRATCH}
cd       ${SCRATCH}
echo "Working from ${PWD}"
run_Selector -f ${SLURM_SUBMIT_DIR}/%(filelist)s  -s %(samplename)s
echo "${PWD} contentents:"
ls -ltrh
# cp -p ${SCRATCH}/susyNt.root ${SLURM_SUBMIT_DIR}/%(outfilename)s # right now no output
echo "Done, `date`"
echo "Cleaning up ${SCRATCH}"
rm -rf $SCRATCH || exit $?
