#!/bin/bash

#BSUB -J madlab_wmaze_preproc_rest
#BSUB -o /scratch/madlab/crash/wmaze_preproc_rest_out
#BSUB -e /scratch/madlab/crash/wmaze_preproc_rest_err

for subj in WMAZE_001; do
python resting_preproc.py -r 1 2 -s ${subj} \
      -d /home/data/madlab/surfaces/wmaze \
      -u 0.009 -l 0.08 \
      --do_slice_times=True \
      --use_fsl_bp=True \
      -t 2.0 -k 5.0 -n 3 \
      -v middle \
      -o /home/data/madlab/data/mri/wmaze_rest/preproc \
      -w /scratch/madlab/wmaze_rest/preproc \
      -p LSF --plugin_args "dict(bsub_args='-q PQ_madlab')"
done

