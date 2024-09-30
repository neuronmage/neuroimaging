#!/usr/bin/env python
import os

# Array containing all subject id numbers

#sids = ['WMAZE_001']


sids = ['WMAZE_002', 'WMAZE_004', 'WMAZE_005', 'WMAZE_006',  
        'WMAZE_007', 'WMAZE_008', 'WMAZE_009', 'WMAZE_010', 'WMAZE_012', 
        'WMAZE_017', 'WMAZE_018', 'WMAZE_019', 'WMAZE_020', 'WMAZE_021', 
        'WMAZE_022', 'WMAZE_023', 'WMAZE_024', 'WMAZE_026', 'WMAZE_027']


# Establish working and output directories
workdir = '/home/data/madlab/scripts/wmaze/anal_MR_thesis/fixed_before_conditional/model3/model3_2-3-2/model3_2-3-2-5/status/merge_copes'
outdir = '/home/data/madlab/data/mri/wmaze/frstlvl/wmaze_MRthesis/fixed_before_conditional/model3_2-3-2-5/merge_copes'

# Iterate through all subjects
for sub in sids:
    # Flexible command to execute main script, providing the kwarg flags with their respective information (subject_id, output_dir, work_dir, run)
    convertcmd = ' '.join(['python', '/home/data/madlab/scripts/wmaze/anal_MR_thesis/fixed_before_conditional/model3/model3_2-3-2/model3_2-3-2-5/merge_copes_model3_2-3-2-5.py', '-s', sub, '-o', outdir, '-w', workdir])
    # The shell file for each subject
    script_file = 'merge_copes_model3_2-3-2-5-{0}.sh'.format(sub)
    # Creates and opens the shell script file
    with open(script_file, 'wt') as fp:
        # Writes the line to identify as bash
        fp.writelines(['#!/bin/bash\n', convertcmd])
    # Submission statement of the shell file to the scheduler
    outcmd = 'bsub -J atm-merge_copes_model3_2-3-2-5-{0} -q PQ_madlab -e /scratch/madlab/crash/wmaze_MRthesis/model3_2-3-2-5/merge_copes/err_{0} -o /scratch/madlab/crash/wmaze_MRthesis/model3_2-3-2-5/merge_copes/out_{0} < {1}'.format(sub, script_file)
    # Execution of the submission command
    os.system(outcmd)
    continue
