#!/usr/bin/env python

import os
subjs = ['001', '002', '004', '005', '006',  
         '007', '008', '009', '010', '012',  
         '017', '018', '019', '020', '021', 
         '022', '023', '024', '026', '027']

workdir = '/home/data/madlab/scripts/wmaze/anal_MR_thesis/fixed_before_conditional/model3/model3_1-3-2-B/status/lvl1'
outdir = '/home/data/madlab/data/mri/wmaze/frstlvl/wmaze_MRthesis/fixed_before_conditional/model3_1-3-2-B'
for i, sid in enumerate(subjs):
    wmazefrstlvlcmd = ' '.join(['python', 'wmaze_lvl1_model3-1-3-2-B.py', '-s', sid, '-g', str(i % 2),
                              '-o', outdir, '-w', workdir])
    script_file = 'wmaze_lvl1_model3-1-3-2-B-{0}.sh'.format(sid)
    with open(script_file, 'wt') as fp:
        fp.writelines(['#!/bin/bash\n'])
        fp.writelines(['#SBATCH --job-name = wmaze_lvl1_model3-1-3-2-B_{0}\n'.format(sid)])
        fp.writelines(['#SBATCH --nodes 1\n'])
        fp.writelines(['#SBATCH --ntasks 1\n'])
        fp.writelines(['#SBATCH -p investor\n'])
        fp.writelines(['#SBATCH --qos pq_madlab\n'])
        fp.writelines(['#SBATCH -e /scratch/madlab/crash/wmaze_MRthesis/model3_1-3-2-B/lvl1/lvl1_model3_1-3-2-B_err_{0}\n'.format(sid)])
        fp.writelines(['#SBATCH -o /scratch/madlab/crash/wmaze_MRthesis/model3_1-3-2-B/lvl1/lvl1_model3_1-3-2-B_out_{0}\n'.format(sid)])
        fp.writelines([dwipreproccmd])
    outcmd = 'sbatch {0}'.format(script_file)
    os.system(outcmd)
    continue
