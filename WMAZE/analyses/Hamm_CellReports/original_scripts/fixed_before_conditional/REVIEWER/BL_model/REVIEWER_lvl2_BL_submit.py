#!/usr/bin/env python
import os

#subjs = ['WMAZE_002']

subjs = ['WMAZE_001', 'WMAZE_002', 'WMAZE_004', 'WMAZE_005', 'WMAZE_006',  
         'WMAZE_007', 'WMAZE_008', 'WMAZE_009', 'WMAZE_010', 'WMAZE_012', 
         'WMAZE_017', 'WMAZE_018', 'WMAZE_019', 'WMAZE_020', 'WMAZE_021', 
         'WMAZE_022', 'WMAZE_023', 'WMAZE_024', 'WMAZE_026', 'WMAZE_027']
    
workdir = '/scratch/madlab/crash/mandy_crash/REVIEWER/lvl2'
outdir = '/home/data/madlab/data/mri/wmaze/scndlvl/wmaze_MRthesis/fixed_before_conditional/REVIEWER'

for i, sid in enumerate(subjs):
    # Flexible command to execute the level 1 script with the kwarg flags and respective information using python
    convertcmd = ' '.join(['python', '/home/data/madlab/scripts/wmaze/anal_MR_thesis/fixed_before_conditional/REVIEWER/BL_model/REVIEWER_lvl2_BL.py', 
                           '-s', sid, '-o', outdir, '-w', workdir])
      
    # Submission statement of the shell file to the SLURM scheduler
    outcmd = 'sbatch -J atm-REVIEWER_lvl2_BL-{0} -p investor --qos pq_madlab \
             -e /scratch/madlab/crash/mandy_crash/REVIEWER/lvl2/err_{0} \
             -o /scratch/madlab/crash/mandy_crash/REVIEWER/lvl2/out_{0} --wrap="{1}"'.format(sid, convertcmd)
    os.system(outcmd)
    continue
