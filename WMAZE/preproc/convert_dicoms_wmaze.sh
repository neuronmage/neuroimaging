#!/bin/bash

#BSUB -J madlab_wmaze_dcm_convert
#BSUB -o /home/data/madlab/scripts/preproc_scripts/dcm_convert_out
#BSUB -e /home/data/madlab/scripts/preproc_scripts/dcm_convert_err

python dicomconvert2.py -d /home/data/madlab/dicoms/wmaze -o /home/data/madlab/data/mri/wmaze -f heuristic_wmaze.py -q PQ_madlab -c dcm2nii -s WMAZE_023

