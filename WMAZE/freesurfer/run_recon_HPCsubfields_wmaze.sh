#! /bin/bash

#BSUB -J madlab_recon_edits

#wmaze_subjlist= `ls /home/data/madlab/surfaces/wmaze`;

wmaze_subjlist='WMAZE_002 WMAZE_004 WMAZE_005 WMAZE_006 WMAZE_007 WMAZE_008 WMAZE_009 WMAZE_010 WMAZE_012 WMAZE_017 WMAZE_018 WMAZE_019 WMAZE_020 WMAZE_021 WMAZE_022 WMAZE_023
WMAZE_024 WMAZE_026 WMAZE_027'



for subj in $wmaze_subjlist; do

# re-running it after wm edits have been made
cmd="recon-all -subjid ${subj} -hippocampal-subfields-T1"

echo `echo ${cmd}` | bsub -q PQ_madlab -e /scratch/madlab/crash/wmaze_MRthesis/ROI/recon_HPCSF_${subj}_err -o /scratch/madlab/crash/wmaze_MRthesis/ROI/recon_HPCSF_${subj}_out

done
