#!/bin/bash

#SBATCH -p %(queue)s
#SBATCH --mem-per-cpu=1gb
#SBATCH --time=04:00:00
#SBATCH -o %(logfile)s
#SBATCH --job-name=%(jobname)s

echo "Starting on `hostname`, `date`"
cd ${SLURM_SUBMIT_DIR}
./python/plot_emu.py %(opt)s
echo "Done, `date`"
