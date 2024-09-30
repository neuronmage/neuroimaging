#!/bin/bash

#BSUB -J madlab_wmaze_preproc
#BSUB -o /scratch/madlab/wmaze/report/wmaze_preproc_out
#BSUB -e /scratch/madlab/wmaze/report/wmaze_preproc_err

for subj in pilot_1; do
python wmaze_preproc.py -r 1 2 3 -s ${subj} \
      -d /home/data/madlab/surfaces/wmaze \
      -u 0.007 -l -1 \
      --do_slice_times=True \
      --use_fsl_bp=False \
      -t 2.0 -k 5.0 -n 3 \
      -v middle \
      -o /home/data/madlab/data/mri/wmaze/preproc \
      -w /scratch/madlab/wmaze/preproc \
      -p LSF --plugin_args "dict(bsub_args='-q PQ_madlab')"
done

