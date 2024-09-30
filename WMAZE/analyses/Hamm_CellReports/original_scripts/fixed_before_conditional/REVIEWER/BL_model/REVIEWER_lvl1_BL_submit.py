#!/usr/bin/env python
import os

#subjs = ['WMAZE_001']

subjs = ['WMAZE_001', 'WMAZE_002', 'WMAZE_004', 'WMAZE_005', 'WMAZE_006',  
         'WMAZE_007', 'WMAZE_008', 'WMAZE_009', 'WMAZE_010', 'WMAZE_012', 
         'WMAZE_017', 'WMAZE_018', 'WMAZE_019', 'WMAZE_020', 'WMAZE_021', 
         'WMAZE_022', 'WMAZE_023', 'WMAZE_024', 'WMAZE_026', 'WMAZE_027']

workdir = '/scratch/madlab/crash/mandy_crash/REVIEWER/lvl1/'
outdir = '/home/data/madlab/data/mri/wmaze/frstlvl/wmaze_MRthesis/fixed_before_conditional/REVIEWER/'

for i, sid in enumerate(subjs):
    convertcmd = ' '.join(['python', '/home/data/madlab/scripts/wmaze/anal_MR_thesis/fixed_before_conditional/REVIEWER/BL_model/REVIEWER_lvl1_BL.py', 
                           '-s', sid, '-o', outdir, '-w', workdir])
      
    outcmd = 'sbatch -J atm-REVIEWER_lvl1_BL-{0} -p investor --qos pq_madlab \
             -e /scratch/madlab/crash/mandy_crash/REVIEWER/lvl1/err_{0} \
             -o /scratch/madlab/crash/mandy_crash/REVIEWER/lvl1/out_{0} --wrap="{1}"'.format(sid, convertcmd)
    os.system(outcmd)
    continue
