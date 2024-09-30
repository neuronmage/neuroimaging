#!/usr/bin/env python
import os

sids1 = ['sub-006', 'sub-007', 'sub-008', 'sub-010', 
        'sub-012', 'sub-013', 'sub-014', 'sub-015', 'sub-016',  
	'sub-018', 'sub-019', 'sub-020', 'sub-021', 'sub-022', 
        'sub-023', 'sub-024', 'sub-025', 'sub-026', 'sub-027',  
	'sub-028', 'sub-029', 'sub-030', 'sub-031', 'sub-032']
sids = ['sub-005', 'sub-030']
   
workdir = '/scratch/madlab/crash/mandy/vcat/MVPA-Loc/lvl1'
outdir = '/home/data/madlab/Mattfeld_vCAT/derivatives/MVPA-Loc/lvl1'

for sid in sids:
    #flexible command to execute script with kwarg flags and respective information using python
    convertcmd = ' '.join(['python', '/home/data/madlab/Mattfeld_vCAT/code/analyses/MVPA-Loc/MVPA-Loc_lvl1.py', '-s', sid, '-o', outdir, '-w', workdir])
    
    #submission statement of shell file to SLURM scheduler
    outcmd = 'sbatch -J lvl1-{0} -p investor --qos pq_madlab --account iacc_madlab\
             -e /scratch/madlab/crash/mandy/vcat/MVPA-Loc/lvl1/err_{0} \
             -o /scratch/madlab/crash/mandy/vcat/MVPA-Loc/lvl1/out_{0} --wrap="{1}"'.format(sid[4:], convertcmd)
    os.system(outcmd)
    continue
