#!/bin/bash
#SBATCH -p investor
#SBATCH --qos pq_madlab
#SBATCH --account iacc_madlab
#SBATCH -o /scratch/madlab/Pruden_SEA/run_recon_out
#SBATCH -e /scratch/madlab/Pruden_SEA/run_recon_err
source ~/.bashrc;
sea_env;
module add glibc/2.14

python /home/data/madlab/Pruden_SEA/code/preproc/SEA_run_recon_EBC.py
