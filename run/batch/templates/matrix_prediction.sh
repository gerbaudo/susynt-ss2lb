#!/bin/bash

# Template to execute `run_MatrixPrediction` on the greenplanet queue
#
# For more details on the formatting, see 'selection.sh' and 'submitJobs.py'.
#
# davide.gerbaudo@gmail.com
# September 2014

#SBATCH -p atlas_all
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
run_MatrixPrediction \
 -i ${SLURM_SUBMIT_DIR}/%(filelist)s \
 -o %(local_outfilename)s \
 -s %(samplename)s \
 --matrix-file ${ROOTCOREDIR}/data/DileptonMatrixMethod/FakeMatrix_Jul_26.root \
 --etapt \
 %(opt)s

 # --matrix-file ${ROOTCOREDIR}/data/DileptonMatrixMethod/FinalFakeHist_May_20.root \


echo "${PWD} contentents:"
ls -ltrh
cp -p ${SCRATCH}/%(local_outfilename)s ${SLURM_SUBMIT_DIR}/%(outfilename)s
echo "Done, `date`"
echo "Cleaning up ${SCRATCH}"
rm -rf $SCRATCH || exit $?
