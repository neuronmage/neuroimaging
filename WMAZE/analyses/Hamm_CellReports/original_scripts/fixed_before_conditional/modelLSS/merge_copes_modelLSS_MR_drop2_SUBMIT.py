#!/usr/bin/env python
import os

# Array containing all subject id numbers

sids = ['WMAZE_001', 'WMAZE_002', 'WMAZE_004', 'WMAZE_005', 'WMAZE_006',  
        'WMAZE_007', 'WMAZE_008', 'WMAZE_009', 'WMAZE_010', 'WMAZE_012', 
        'WMAZE_017', 'WMAZE_018', 'WMAZE_019', 'WMAZE_020', 'WMAZE_021', 
        'WMAZE_022', 'WMAZE_023', 'WMAZE_024', 'WMAZE_026', 'WMAZE_027']

#sids = ['WMAZE_001']

# Establish working and output directories
workdir = '/scratch/madlab/crash/wmaze_MRthesis/modelLSS_MR_drop3/merge_copes'
outdir = '/home/data/madlab/data/mri/wmaze/frstlvl/wmaze_MRthesis/fixed_before_conditional/modelLSS_MR_drop3/merge_copes3'

# Iterate through all subjects
for sub in sids:
    # Flexible command to execute main script, providing the kwarg flags with their respective information (subject_id, output_dir, work_dir, run)
    convertcmd = ' '.join(['python', '/home/data/madlab/scripts/wmaze/anal_MR_thesis/fixed_before_conditional/modelLSS_MR/merge_copes_modelLSS_MR_drop3.py', '-s', sub, '-o', outdir, '-w', workdir])
    # The shell file for each subject
    script_file = 'merge_copes_modelLSS_MR_drop3-{0}.sh'.format(sub)
    # Creates and opens the shell script file
    with open(script_file, 'wt') as fp:
        # Writes the line to identify as bash
        fp.writelines(['#!/bin/bash\n', convertcmd])
    # Submission statement of the shell file to the scheduler
    outcmd = 'bsub -J merge_copes_modelLSS_MR_drop3{0} -q PQ_madlab -e {2}/merge_copes3_err_{0} -o {2}/merge_copes3_out_{0} < {1}'.format(sub, script_file, workdir)
    # Execution of the submission command
    os.system(outcmd)
    continue
