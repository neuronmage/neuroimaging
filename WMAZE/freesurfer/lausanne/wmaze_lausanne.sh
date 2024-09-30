#! /bin/bash


# This will apply the atlas of about 500 ROIs extracted from The Connectome Mapper Toolkit (http://www.cmtk.org/) 
# and embed that atlas with a randomized colortable for ease of viewing.

# set SUBJECT_DIR
export SUBJECTS_DIR=/home/data/madlab/surfaces/wmaze

# change path to atlases location
export atlasdir=/data/madlab/scripts/freesurfer_scripts/lausanne/cmp_nipype.2.1.0/cmtklib/data/colortable_and_gcs


# Subjects list
wmaze_subs='WMAZE_001 WMAZE_002 WMAZE_004 WMAZE_005 WMAZE_006 WMAZE_007 WMAZE_008 WMAZE_009 WMAZE_010 WMAZE_012 WMAZE_017 WMAZE_018 WMAZE_019 WMAZE_020 WMAZE_021 WMAZE_022 WMAZE_023
WMAZE_024 WMAZE_026 WMAZE_027'

for subj in $wmaze_subs; do
	echo $subj
        # Specify atlases: other_atlases  "P1_16" "P17_28"  "p29_36"
        for atlas in 36 60 125 250; do

            # create aparc.annot
			mris_ca_label $subj lh lh.sphere.reg $atlasdir/my_atlas_gcs/myatlas_${atlas}_lh.gcs $SUBJECTS_DIR\/$subj\/label/lh.lausanne2008_$atlas\_aparc.annot
			mris_ca_label $subj rh rh.sphere.reg $atlasdir/my_atlas_gcs/myatlas_${atlas}_rh.gcs $SUBJECTS_DIR\/$subj\/label/rh.lausanne2008_$atlas\_aparc.annot
			
			# Create a new directory to place the extracted labels
			mkdir ${SUBJECTS_DIR}/${subj}/label/lausanne2008_${atlas}_label
			
			# Extract atlas label files 
			mri_annotation2label --annotation lausanne2008_${atlas}_aparc --hemi lh --subject $subj --outdir ${SUBJECTS_DIR}/${subj}/label/lausanne2008_${atlas}_label
			mri_annotation2label --annotation lausanne2008_${atlas}_aparc --hemi rh --subject $subj --outdir ${SUBJECTS_DIR}/${subj}/label/lausanne2008_${atlas}_label
			
			# Remove the previously created lausanne2008_aparc.annot.  
			# Another one will be created in the next step.
			rm ${SUBJECTS_DIR}/${subj}/label/lh.lausanne2008_${atlas}_aparc.annot
			rm ${SUBJECTS_DIR}/${subj}/label/rh.lausanne2008_${atlas}_aparc.annot
			
			# Create new lausanne2008_aparc.annot files and embed correct look up table
			mris_label2annot --ctab $atlasdir\/original_color_${atlas}_l.txt --subject $subj --ldir ${SUBJECTS_DIR}/${subj}/label/lausanne2008_${atlas}_label --no-unknown --annot lausanne2008_${atlas}_aparc --hemi lh
			mris_label2annot --ctab $atlasdir\/original_color_${atlas}_r.txt --subject $subj --ldir ${SUBJECTS_DIR}/${subj}/label/lausanne2008_${atlas}_label --no-unknown --annot lausanne2008_${atlas}_aparc --hemi rh
			
			# Writing output aseg to ${SUBJECTS_DIR}/${subj}/mri/lausanne2008_$atlas_aparc+aseg.mgz
			
			mri_aparc2aseg --s $subj --annot lausanne2008_${atlas}_aparc
			done
done
