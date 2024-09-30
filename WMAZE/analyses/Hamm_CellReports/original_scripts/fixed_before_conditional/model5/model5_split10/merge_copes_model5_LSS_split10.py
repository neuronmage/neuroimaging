#!/usr/bin/env python

"""
============================================================================================
wmaze_fMRI: Mandy Thesis -- Fixed before Condition - Merge COPES for beta series correlation
============================================================================================
Very simple Nipype script to merge beta series copes from wmaze_lvl1_model5_LSS_split10.py

This workflow makes use of:
- FSL
For example:
  python model3_merge_copes.py -s WMAZE_001
                                        -o /home/data/madlab/data/mri/wmaze/frstlvl/wmaze_MRthesis/fixed_before_conditional/model5_LSS_split10/merge_copes
                                        -w /scratch/madlab/crash/wmaze_MRthesis/model5_LSS_split10/merge_copes
"""

import os
from nipype.pipeline.engine import Workflow, Node, MapNode
from nipype.interfaces.utility import Function
from nipype.interfaces.io import DataGrabber
from nipype.interfaces.io import DataSink
from nipype.interfaces.fsl import Merge

# Function to merge copes for each participant
def cope_merge_wf(subject_id, 
	          sink_directory,
	          name = 'cope_merge_wf'): 
        # Name workflow   
	cope_merge_wf = Workflow(name = 'cope_merge_wf')
        # Dictionary for Datagrabber
    	info = dict(frst10_corr_cope = [['subject_id']],
		    frst10_incorr_cope = [['subject_id']],
                    last10_corr_cope = [['subject_id']],
		    last10_incorr_cope = [['subject_id']])

        # Node to grab corr and incorr cope files
	dsource = Node(DataGrabber(infields = ['subject_id'], 
                                   outfields = info.keys()), 
                       name = 'dsource')
	dsource.inputs.template = '*'
	dsource.inputs.subject_id = subject_id
	dsource.inputs.base_directory = os.path.abspath('/home/data/madlab/data/mri/wmaze/frstlvl/wmaze_MRthesis/fixed_before_conditional/model5_LSS_split10/')
	dsource.inputs.field_template = dict(frst10_corr_cope = '%s/modelfit/contrasts/_estimate_model*/cope*_frst10_FX_before_COND_corr_set*_trl*.nii.gz',
				             frst10_incorr_cope = '%s/modelfit/contrasts/_estimate_model*/cope*_frst10_FX_before_COND_incorr_set*_trl*.nii.gz',
                                             last10_corr_cope = '%s/modelfit/contrasts/_estimate_model*/cope*_last10_FX_before_COND_corr_set*_trl*.nii.gz',
				             last10_incorr_cope = '%s/modelfit/contrasts/_estimate_model*/cope*_last10_FX_before_COND_incorr_set*_trl*.nii.gz')
	dsource.inputs.template_args = info
	dsource.inputs.sort_filelist = True
	dsource.inputs.ignore_exception = False
	dsource.inputs.raise_on_empty = True


	# Node to merge frst10 CORRECT trials across all 3 sets
	merge_frst10_corr = Node(Merge(), 
                          name = 'merge_frst10_corr')
	merge_frst10_corr.inputs.dimension = 't'
	merge_frst10_corr.inputs.output_type = 'NIFTI_GZ'
        merge_frst10_corr.inputs.merged_file = 'frst10_cope_corr.nii.gz'
	merge_frst10_corr.inputs.tr = 2.00
	cope_merge_wf.connect(dsource, 'frst10_corr_cope', merge_frst10_corr, 'in_files')

        #### dsource (frst10_corr_cope) ----> merge_corr (in_files)


        # Node to merge frst10 INCORRECT trials across all 3 sets
	merge_frst10_incorr = Node(Merge(), 
                            name = 'merge_frst10_incorr')
	merge_frst10_incorr.inputs.dimension = 't'
	merge_frst10_incorr.inputs.output_type = 'NIFTI_GZ'
        merge_frst10_incorr.inputs.merged_file = 'frst10_cope_incorr.nii.gz'
	merge_frst10_incorr.inputs.tr = 2.00
	cope_merge_wf.connect(dsource, 'frst10_incorr_cope', merge_frst10_incorr, 'in_files')
 
         #### dsource (incorr_cope) ----> merge_incorr (in_files)

        # Node to merge last10 CORRECT trials across all 3 sets
	merge_last10_corr = Node(Merge(), 
                          name = 'merge_last10_corr')
	merge_last10_corr.inputs.dimension = 't'
	merge_last10_corr.inputs.output_type = 'NIFTI_GZ'
        merge_last10_corr.inputs.merged_file = 'last10_cope_corr.nii.gz'
	merge_last10_corr.inputs.tr = 2.00
	cope_merge_wf.connect(dsource, 'last10_corr_cope', merge_last10_corr, 'in_files')

        #### dsource (corr_cope) ----> merge_corr (in_files)


        # Node to merge last10 INCORRECT trials across all 3 sets
	merge_last10_incorr = Node(Merge(), 
                            name = 'merge_last10_incorr')
	merge_last10_incorr.inputs.dimension = 't'
	merge_last10_incorr.inputs.output_type = 'NIFTI_GZ'
        merge_last10_incorr.inputs.merged_file = 'last10_cope_incorr.nii.gz'
	merge_last10_incorr.inputs.tr = 2.00
	cope_merge_wf.connect(dsource, 'last10_incorr_cope', merge_last10_incorr, 'in_files')
 
         #### dsource (incorr_cope) ----> merge_incorr (in_files)

        # Node to output data
	dsink = Node(DataSink(), 
                     name = 'dsink')
	dsink.inputs.base_directory = sink_directory 
	dsink.inputs.container = subject_id
	cope_merge_wf.connect(merge_frst10_corr, 'merged_file', dsink, 'merged.@frst10_corr')
	cope_merge_wf.connect(merge_frst10_incorr, 'merged_file', dsink, 'merged.@frst10_incorr')
	cope_merge_wf.connect(merge_last10_corr, 'merged_file', dsink, 'merged.@last10_corr')
	cope_merge_wf.connect(merge_last10_incorr, 'merged_file', dsink, 'merged.@last10_incorr')
	
        #### merge_corr (merged_file) ----> dsink (corr)
        #### merge_incorr (merged_file) ----> dsink (incorr)

	return cope_merge_wf

	
"""
Creates the full workflow
"""


def create_workflow(args, 
                    name = 'model3_merge_copes'):
	kwargs = dict(subject_id = args.subject_id,
		      sink_directory = os.path.abspath(args.out_dir),
		      name = name)
    	cope_workflow = cope_merge_wf(**kwargs)
    	return cope_workflow


if __name__ == "__main__":
    	from argparse import ArgumentParser
	parser = ArgumentParser(description =__doc__)

    	parser.add_argument("-s", "--subject_id", 
                            dest = "subject_id",
                            help = "Current subject id", 
                            required = True)

    	parser.add_argument("-o", "--output_dir", 
                            dest = "out_dir",
                            help = "Output directory")

    	parser.add_argument("-w", "--work_dir", 
                            dest = "work_dir",
                            help = "Working directory")

    	args = parser.parse_args()

    	wf = create_workflow(args)

    	if args.work_dir:
        	work_dir = os.path.abspath(args.work_dir)
    	else:
        	work_dir = os.getcwd()

    	wf.config['execution']['crashdump_dir'] = '/scratch/madlab/crash/wmaze_MRthesis/model5_LSS_split10/merge_copes'
    	wf.base_dir = work_dir + '/' + args.subject_id
    	wf.run(plugin = 'LSF', plugin_args = {'bsub_args': '-q PQ_madlab'})
