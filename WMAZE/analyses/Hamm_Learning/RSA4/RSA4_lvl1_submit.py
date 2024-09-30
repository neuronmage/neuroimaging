#!/usr/bin/env python
import os

subsX = ['WMAZE_001', 'WMAZE_002', 'WMAZE_004', 'WMAZE_005', 'WMAZE_006', 'WMAZE_007', 'WMAZE_008', 'WMAZE_009', 'WMAZE_010', 'WMAZE_012', #all Ss
        'WMAZE_017', 'WMAZE_018', 'WMAZE_019', 'WMAZE_020', 'WMAZE_021', 'WMAZE_022', 'WMAZE_023', 'WMAZE_024', 'WMAZE_026', 'WMAZE_027']
subs = ['WMAZE_001'] #Ss who failed initial exec

workdir = '/scratch/madlab/crash/mandy/learning/RSA4' #dir for workflow
outdir = '/scratch/madlab/crash/mandy/learning/RSA4/lvl1' #dir for final data

for i, sub in enumerate(subs): #iterate through Ss
    convertcmd = ' '.join(['python', '/home/data/madlab/Mattfeld_WMAZE/code/analyses/Hamm_Learning/RSA4/RSA4_lvl1.py', 
                           '-s', sub, '-o', outdir, '-w', workdir]) #kwargs fed to lvl1 script via ArgumentParser
    #string to submit each Ss job, defining SLURM parameters and directing error and output files to desired directories 
    outcmd = 'sbatch -J RSA4_lvl1-{0} --partition investor --qos pq_madlab --account iacc_madlab \
             -e /scratch/madlab/crash/mandy/learning/RSA4/err_{0} \
             -o /scratch/madlab/crash/mandy/learning/RSA4/out_{0} --wrap="{1}"'.format(sub, convertcmd)
    os.system(outcmd) #exec of job command
    continue 
