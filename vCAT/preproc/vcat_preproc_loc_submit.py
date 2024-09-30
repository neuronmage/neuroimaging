#!/usr/bin/env python

import os


sids = ['sub-014', 'sub-015', 'sub-016', 'sub-019', 'sub-020', 'sub-028']
	

workdir = '/scratch/madlab/crash/mandy/vcat/preproc_loc'
outdir = '/home/data/madlab/Mattfeld_vCAT/derivatives/preproc_loc'
surf_dir = '/home/data/madlab/Mattfeld_vCAT/derivatives/freesurfer'

for i, sid in enumerate(sids):
    convertcmd = ' '.join(['python', 'vcat_preproc_loc.py',
                           '-s', sid, #subject id
                           '-o', outdir, #output directory
                           '-w', workdir, #working directory - pipeline output (scratch/madlab/crash)
                           '-r 1 2', #run numbers
                           '-d', surf_dir, #FS surface directory
                           '-u 0.007 -l -1', #high and low-pass freq
                           '--do_slice_times=True', #toggle Nipy slice-timing correction
                           '--use_fsl_bp=False', #toggle FSL band-pass filtering
                           '-t 1.76 -k 5.0 -n 3', #TR, FWHM smoothing kernel, # of components to extract from ACompCor
                           '-v first']) #which reference volume
    outcmd = 'sbatch -J vCAT-{0} --partition IB_16C_96G --qos pq_madlab --account iacc_madlab \
             -e /scratch/madlab/crash/mandy/vcat/preproc_loc/err_{0} \
             -o /scratch/madlab/crash/mandy/vcat/preproc_loc/out_{0} --wrap="{1}"'.format(sid[4:], convertcmd)
    os.system(outcmd)
    continue

