#!/usr/bin/env python
import os

subs = ['WMAZE_001', 'WMAZE_002', 'WMAZE_004', 'WMAZE_005', 'WMAZE_006', 'WMAZE_007', 'WMAZE_008', 'WMAZE_009', 'WMAZE_010', 'WMAZE_012', 
        'WMAZE_017', 'WMAZE_018', 'WMAZE_019', 'WMAZE_020', 'WMAZE_021', 'WMAZE_022', 'WMAZE_023', 'WMAZE_024', 'WMAZE_026', 'WMAZE_027']

workdir = '/scratch/madlab/crash/mandy/learning/PMA/fixed/lvl2/status'
outdir = '/home/data/madlab/data/mri/wmaze/scndlvl/learning/PMA/fixed'

for i, sub in enumerate(subs):
    convertcmd = ' '.join(['python', '/home/data/madlab/scripts/wmaze/mandy/learning/PMA/PMA_fixed_lvl2.py', '-s', sub, '-o', outdir, '-w', workdir])
    outcmd = 'sbatch -J PMA_fixed_lvl2-{0} -p investor --qos pq_madlab \
             -e /scratch/madlab/crash/mandy/learning/PMA/fixed/lvl2/err_{0} \
             -o /scratch/madlab/crash/mandy/learning/PMA/fixed/lvl2/out_{0} --wrap="{1}"'.format(sub, convertcmd)
    os.system(outcmd)
    continue
