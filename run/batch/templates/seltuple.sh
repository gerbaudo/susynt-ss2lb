#!/bin/bash

# Template to execute `run_Selector` and make nutples on the greenplanet queue
#
# For more details on the formatting, see 'selection.sh' and 'submitJobs.py'.
#
# davide.gerbaudo@gmail.com
# September 2014

#SBATCH -p %(queue)s
#SBATCH --mem-per-cpu=1gb
#SBATCH --time=04:00:00
#SBATCH -o %(logfile)s
#SBATCH --job-name=%(jobname)s

SCRATCH=/scratch/${USER}/${SLURM_JOB_ID}

echo "Starting on `hostname`, `date`"
mkdir -p ${SCRATCH}
cd       ${SCRATCH}
echo "Working from ${PWD}"

stdbuf --output=1M --error=1M \
 run_Selector \
 --sample %(samplename)s \
 --input ${SLURM_SUBMIT_DIR}/%(filelist)s \
 --tuple-out %(local_outfilename)s \
 %(opt)s


echo "${PWD} contentents:"
ls -ltrh
cp -p ${SCRATCH}/%(local_outfilename)s ${SLURM_SUBMIT_DIR}/%(outfilename)s
echo "Done, `date`"
echo "Cleaning up ${SCRATCH}"
rm -rf $SCRATCH || exit $?

