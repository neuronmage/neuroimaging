#!/usr/bin/env python

import os

sids = ['sub-010', 'sub-022', 'sub-028', 'sub-029', 'sub-031']

sids1 = ['sub-005', 'sub-006', 'sub-007', 'sub-008', 'sub-010', 
        'sub-012', 'sub-013', 'sub-014', 'sub-015', 'sub-016',  
	'sub-018', 'sub-019', 'sub-020', 'sub-021', 'sub-022', 
        'sub-023', 'sub-024', 'sub-025', 'sub-026', 'sub-027',  
	'sub-028', 'sub-029', 'sub-030', 'sub-031', 'sub-032']


workdir = '/scratch/madlab/crash/mandy/vcat/preproc'
outdir = '/home/data/madlab/Mattfeld_vCAT/derivatives/preproc'
surf_dir = '/home/data/madlab/Mattfeld_vCAT/derivatives/freesurfer'

for i, sid in enumerate(sids):
    convertcmd = ' '.join(['python', 'vcat_preproc.py',
                           '-s', sid, #subject id
                           '-o', outdir, #output directory
                           '-w', workdir, #working directory - pipeline output (scratch/madlab/crash)
                           '-r 1 2 3 4', #run numbers
                           '-d', surf_dir, #FS surface directory
                           '-u 0.007 -l -1', #high and low-pass freq
                           '--do_slice_times=True', #toggle Nipy slice-timing correction
                           '--use_fsl_bp=False', #toggle FSL band-pass filtering
                           '-t 1.76 -k 5.0 -n 3', #TR, FWHM smoothing kernel, # of components to extract from ACompCor
                           '-v first']) #which reference volume
    outcmd = 'sbatch -J vCAT-{0} --partition IB_16C_96G --qos pq_madlab --account iacc_madlab \
             -e /scratch/madlab/crash/mandy/vcat/preproc/err_{0} \
             -o /scratch/madlab/crash/mandy/vcat/preproc/out_{0} --wrap="{1}"'.format(sid[4:], convertcmd)
    os.system(outcmd)
    continue

