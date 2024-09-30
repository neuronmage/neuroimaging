#!/bin/bash
#SBATCH --qos pq_madlab
#SBATCH --account iacc_madlab
#SBATCH --partition 16C_128G
#SBATCH --tasks 1
#SBATCH --cpus-per-task 16
#SBATCH --job-name sing-3.5.3
#SBATCH --output /scratch/madlab/Pruden_SEA/mandy/dwi_preproc/sub-1017_sing-3.5.3_out
#SBATCH --error /scratch/madlab/Pruden_SEA/mandy/dwi_preproc/sub-1017_sing-3.5.3_err

#######################################################
# Set up environmental variables.
#######################################################

export NPROCS=`echo $SLURM_JOB_NODELIST | wc -w`
export OMP_NUM_THREADS=$NPROCS
#OR
#export OMP_NUM_THREADS=4

. $MODULESHOME/../global/profile.modules
module load singularity-3.5.3
#######################################################
#######################################################

bash ABCD_DWI_processing.bash \
-u=/scratch/madlab/Pruden_SEA/mandy/dwi_preproc/dset/sub-1017/dwi/sub-1017_dwi.nii.gz \
-ub=/scratch/madlab/Pruden_SEA/mandy/dwi_preproc/dset/sub-1017/dwi/sub-1017_dwi.bval \
-uv=/scratch/madlab/Pruden_SEA/mandy/dwi_preproc/dset/sub-1017/dwi/sub-1017_dwi.bvec \
-d=/scratch/madlab/Pruden_SEA/mandy/dwi_preproc/dset/sub-1017/dwi/sub-1017_dwi-PA.nii.gz \
-db=/scratch/madlab/Pruden_SEA/mandy/dwi_preproc/dset/sub-1017/dwi/sub-1017_dwi-PA.bval \
-dv=/scratch/madlab/Pruden_SEA/mandy/dwi_preproc/dset/sub-1017/dwi/sub-1017_dwi-PA.bvec \
-T1=/scratch/madlab/Pruden_SEA/mandy/dwi_preproc/dset/sub-1017/anat/sub-1017_T1w.nii.gz \
-T2=/scratch/madlab/Pruden_SEA/mandy/dwi_preproc/dset/sub-1017/anat/sub-1017_T2-tse.nii.gz \
--gw_coeffs=/scratch/madlab/Pruden_SEA/mandy/dwi_preproc/scripts/coeff.grad \
--output_res=1 \
-p=/scratch/madlab/Pruden_SEA/mandy/dwi_preproc/temp/sub-1017/sing-3.5.3 \
-o=/scratch/madlab/Pruden_SEA/mandy/dwi_preproc/output/sub-1017/sing-3.5.3
