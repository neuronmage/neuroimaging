#!/bin/bash
#SBATCH --qos pq_madlab
#SBATCH --account iacc_madlab
#SBATCH --partition IB_44C_512G
#SBATCH --tasks 1
#SBATCH --cpus-per-task 16
#SBATCH --job-name SEA-1016
#SBATCH --output /scratch/madlab/Pruden_SEA/mandy/dwi_preproc/sub-1016_out
#SBATCH --error /scratch/madlab/Pruden_SEA/mandy/dwi_preproc/sub-1016_err

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
sub=1016
workdir=/scratch/madlab/Pruden_SEA/mandy/dwi_preproc

bash ABCD_DWI_processing.bash \
-u=$workdir/dset/sub-$sub/dwi/sub-"$sub"_dwi.nii.gz \
-ub=$workdir/dset/sub-$sub/dwi/sub-"$sub"_dwi.bval \
-uv=$workdir/dset/sub-$sub/dwi/sub-"$sub"_dwi.bvec \
-d=$workdir/dset/sub-$sub/dwi/sub-"$sub"_dwi-PA.nii.gz \
-db=$workdir/dset/sub-$sub/dwi/sub-"$sub"_dwi-PA.bval \
-dv=$workdir/dset/sub-$sub/dwi/sub-"$sub"_dwi-PA.bvec \
-T1=$workdir/dset/sub-$sub/anat/sub-"$sub"_T1w.nii.gz \
-T2=$workdir/dset/sub-$sub/anat/sub-"$sub"_T2-tse.nii.gz \
--gw_coeffs=$workdir/scripts/coeff.grad \
-p=$workdir/temp/sub-"$sub" \
-o=$workdir/output/sub-"$sub"
