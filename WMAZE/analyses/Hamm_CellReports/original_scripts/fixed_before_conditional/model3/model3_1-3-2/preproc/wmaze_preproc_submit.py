#!/usr/bin/env python

import os

#subjs = ['WMAZE_001']
subjs = ['WMAZE_002', 'WMAZE_003', 'WMAZE_004', 'WMAZE_005', 'WMAZE_006', 'WMAZE_007',
         'WMAZE_008', 'WMAZE_009', 'WMAZE_010', 'WMAZE_012', 'WMAZE_017', 'WMAZE_018',
         'WMAZE_019', 'WMAZE_020', 'WMAZE_021', 'WMAZE_022', 'WMAZE_023', 'WMAZE_024',
         'WMAZE_026', 'WMAZE_027']

workdir = '/scratch/madlab/wmaze/preproc'
outdir = '/home/data/madlab/data/mri/wmaze/preproc'
surf_dir = '/home/data/madlab/surfaces/wmaze'
for i, sid in enumerate(subjs):
    convertcmd = ' '.join(['python', 'wmaze_preproc.py',
                           '-s', sid,
                           '-o', outdir,
                           '-w', workdir,
                           '-r 1 2 3 4 5 6',
                           '-d', surf_dir,
                           '-u 0.007 -l -1',
                           '--do_slice_times=True',
                           '--use_fsl_bp=False',
                           '-t 2.0 -k 5.0 -n 3',
                           '-v middle'])
    script_file = 'preproc-%s.sh' % sid
    with open(script_file, 'wt') as fp:
        fp.writelines(['#!/bin/bash\n', convertcmd])
    outcmd = 'bsub -J atm-preproc-%s -q PQ_madlab -e /scratch/madlab/crash/wmaze_preproc_err_%s -o /scratch/madlab/crash/wmaze_preproc_out_%s < %s' % (sid, sid, sid, script_file)
    os.system(outcmd)
    continue

