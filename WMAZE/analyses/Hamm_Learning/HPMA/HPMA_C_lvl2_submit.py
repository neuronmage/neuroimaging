#!/usr/bin/env python
import os

subjs = ['WMAZE_001', 'WMAZE_002', 'WMAZE_004', 'WMAZE_005', 'WMAZE_006', 'WMAZE_007', 'WMAZE_008', 'WMAZE_009', 'WMAZE_010', 'WMAZE_012', 
         'WMAZE_017', 'WMAZE_018', 'WMAZE_019', 'WMAZE_020', 'WMAZE_021', 'WMAZE_022', 'WMAZE_023', 'WMAZE_024', 'WMAZE_026', 'WMAZE_027']

workdir = '/scratch/madlab/crash/mandy/learning/HPMA_C/lvl2'
outdir = '/home/data/madlab/data/mri/wmaze/scndlvl/learning/HPMA/cond/'

for i, sid in enumerate(subjs):
    convertcmd = ' '.join(['python', '/home/data/madlab/scripts/wmaze/mandy/learning/HPMA/HPMA_C_lvl2.py', '-s', sid, '-o', outdir, '-w', workdir])      
    outcmd = 'sbatch -J HPMA_lvl2-{0} -p investor --qos pq_madlab \
             -e /scratch/madlab/crash/mandy/learning/HPMA_C/lvl2/err_{0} \
             -o /scratch/madlab/crash/mandy/learning/HPMA_C/lvl2/out_{0} --wrap="{1}"'.format(sid, convertcmd)
    os.system(outcmd)
    continue
