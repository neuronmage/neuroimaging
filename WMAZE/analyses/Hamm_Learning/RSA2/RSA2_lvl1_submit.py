#!/usr/bin/env python
import os

subs1 = ['WMAZE_001', 'WMAZE_002', 'WMAZE_004', 'WMAZE_005', 'WMAZE_006', 'WMAZE_007', 'WMAZE_008', 'WMAZE_009', 'WMAZE_010', 'WMAZE_012']
subs = ['WMAZE_017', 'WMAZE_018', 'WMAZE_019', 'WMAZE_020', 'WMAZE_021', 'WMAZE_022', 'WMAZE_023', 'WMAZE_024', 'WMAZE_026', 'WMAZE_027']

workdir = '/scratch/madlab/crash/mandy/Mattfeld_WMAZE/RSA2'
outdir = '/home/data/madlab/Mattfeld_WMAZE/Hamm_Learning/RSA2/lvl1'

for i, sub in enumerate(subs):
    convertcmd = ' '.join(['python', '/home/data/madlab/Mattfeld_WMAZE/code/analyses/Hamm_Learning/RSA2/RSA2_lvl1.py', 
                           '-s', sub, '-o', outdir, '-w', workdir])
    outcmd = 'sbatch -J RSA2_lvl1-{0} --partition centos7_16C_128G --qos pq_madlab --account iacc_madlab \
             -e /scratch/madlab/crash/mandy/Mattfeld_WMAZE/RSA2/err_{0} \
             -o /scratch/madlab/crash/mandy/Mattfeld_WMAZE/RSA2/out_{0} --wrap="{1}"'.format(sub, convertcmd)
    os.system(outcmd)
    continue
