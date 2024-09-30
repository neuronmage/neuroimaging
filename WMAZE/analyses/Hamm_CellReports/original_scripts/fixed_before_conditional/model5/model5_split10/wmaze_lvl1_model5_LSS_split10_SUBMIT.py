#!/usr/bin/env python
import os

#subjs = ['WMAZE_001']

#subjs = ['WMAZE_001', 'WMAZE_002', 'WMAZE_004', 'WMAZE_005', 'WMAZE_006',  
#         'WMAZE_007', 'WMAZE_008', 'WMAZE_009', 'WMAZE_010', 'WMAZE_012', 
#         'WMAZE_017', 'WMAZE_018', 'WMAZE_019', 'WMAZE_020', 'WMAZE_021', 
#         'WMAZE_022', 'WMAZE_023', 'WMAZE_024', 'WMAZE_026', 'WMAZE_027']

subjs = ['WMAZE_007', 'WMAZE_010', 'WMAZE_012', 'WMAZE_021', 'WMAZE_024']

workdir = '/home/data/madlab/scripts/wmaze/anal_MR_thesis/fixed_before_conditional/model5/status/lvl1_LSS_split10'
outdir = '/home/data/madlab/data/mri/wmaze/frstlvl/wmaze_MRthesis/fixed_before_conditional/model5_LSS_split10'

for i, sid in enumerate(subjs):
    # Flexible command to execute the level 1 script with the kwarg flags and respective information using python
    convertcmd = ' '.join(['python', '/home/data/madlab/scripts/wmaze/anal_MR_thesis/fixed_before_conditional/model5/model5_split10/wmaze_lvl1_model5_LSS_split10.py', '-s', sid, '-o', outdir, '-w', workdir])
    
    # The shell file for each subject
    script_file = 'wmaze_lvl1_model5_LSS_split10-{0}.sh'.format(sid)
    
    # Creates and opens the shell script file
    with open(script_file, 'wt') as fp:
        # Writes the line to identify as bash
        fp.writelines(['#!/bin/bash\n', convertcmd])
        
    # Submission statement of the shell file to the scheduler
    outcmd = 'bsub -J atm-wmaze_lvl1_model5_LSS_split10-{0} -q PQ_madlab -e /scratch/madlab/crash/wmaze_MRthesis/model5_LSS_split10/lvl1/err_{1} -o /scratch/madlab/crash/wmaze_MRthesis/model5_LSS_split10/lvl1/out_{2} < {3}'.format(sid, sid, sid, script_file)
    
    # Execution of the submission command
    os.system(outcmd)
    continue
