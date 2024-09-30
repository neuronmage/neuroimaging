#!/usr/bin/env python
import os

subs = ['WMAZE_001', 'WMAZE_002', 'WMAZE_004', 'WMAZE_005', 'WMAZE_006', 'WMAZE_007', 'WMAZE_008', 'WMAZE_009', 'WMAZE_010', 'WMAZE_012', 
        'WMAZE_017', 'WMAZE_018', 'WMAZE_019', 'WMAZE_020', 'WMAZE_021', 'WMAZE_022', 'WMAZE_023', 'WMAZE_024', 'WMAZE_026', 'WMAZE_027']

workdir = '/scratch/madlab/crash/mandy/learning/HALVES/cond/lvl1'
outdir = '/home/data/madlab/Mattfeld_WMAZE/dset/analyses/learning/HALVES/lvl1/cond/1'

for i, sub in enumerate(subs):
    convertcmd = ' '.join(['python', '/home/data/madlab/Mattfeld_WMAZE/code/analyses/learning/HALVES/HALVES_C_lvl1.py', 
                           '-s', sub, '-o', outdir, '-w', workdir])
    outcmd = 'sbatch -J HALVES_C_lvl1-{0} --partition IB_16C_96G --qos pq_madlab --account iacc_madlab \
             -e /scratch/madlab/crash/mandy/learning/HALVES/cond/lvl1/1/err_{0} \
             -o /scratch/madlab/crash/mandy/learning/HALVES/cond/lvl1/1/out_{0} --wrap="{1}"'.format(sub, convertcmd)
    os.system(outcmd)
    continue
