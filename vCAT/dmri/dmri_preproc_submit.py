#!/usr/bin/env python
import os

subs = ['105', '106', '107', '108', '110', '112', '113', '114', '115', '116',
        '118', '119', '120', '121', '122', '123', '124', '126', '126', '127',
        '128', '129', '130', '131', '132']

work_dir = '/scratch/madlab/crash/mandy/Mattfeld_vCAT/dmri'
out_dir = '/home/data/madlab/Mattfeld_vCAT/derivatives/preproc/dmri'

for i, sub in enumerate(subs):
    convertcmd = ' '.join(['python', '/home/data/madlab/Mattfeld_vCAT/code/dmri/dmri_preproc.py', 
                           '-s', sub, '-g', str(i % 2), '-o', out_dir, '-w', work_dir])
    outcmd = 'sbatch -J vCAT_dmri_preproc-{0} --partition centos7 --qos pq_madlab --account iacc_madlab \
                     -e /scratch/madlab/crash/mandy/Mattfeld_vCAT/dmri/dmripreproc_err_{0} \
                     -o /scratch/madlab/crash/mandy/Mattfeld_vCAT/dmri/dmripreproc_out_{0} --wrap="{1}"'.format(sub, convertcmd)
    os.system(outcmd)
    continue

