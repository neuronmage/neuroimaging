#!/usr/bin/env python

"""
============================================================================================
wmaze_fMRI: Mandy Thesis -- Fixed before Condition - Merge COPES for beta series correlation
============================================================================================
Very simple Nipype script to merge beta series copes from wmaze_lvl1_modelLSS_atm12292017.py

This workflow makes use of:
- FSL
For example:
  python model3_merge_copes.py -s WMAZE_001
                                            -o /home/data/madlab/data/mri/wmaze/frstlvl/wmaze_MRthesis/fixed_before_conditional/model5_LSS/merge_copes
                                            -w /scratch/madlab/crash/wmaze_MRthesis/model5_LSS/merge_copes
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
    	info = dict(corr_cope = [['subject_id']],
		    incorr_cope = [['subject_id']])

        # Node to grab corr and incorr cope files
	dsource = Node(DataGrabber(infields = ['subject_id'], 
                                   outfields = info.keys()), 
                       name = 'dsource')
	dsource.inputs.template = '*'
	dsource.inputs.subject_id = subject_id
	dsource.inputs.base_directory = os.path.abspath('/home/data/madlab/data/mri/wmaze/frstlvl/wmaze_MRthesis/fixed_before_conditional/model5_LSS/')
	dsource.inputs.field_template = dict(corr_cope = '%s/modelfit/contrasts/_estimate_model*/cope*_FX_before_COND_corr_run*_trl*.nii.gz',
				             incorr_cope = '%s/modelfit/contrasts/_estimate_model*/cope*_FX_before_COND_incorr_run*_trl*.nii.gz')
	dsource.inputs.template_args = info
	dsource.inputs.sort_filelist = True
	dsource.inputs.ignore_exception = False
	dsource.inputs.raise_on_empty = True


	# Node to merge CORRECT trials across all 6 runs
	merge_corr = Node(Merge(), 
                          name = 'merge_corr')
	merge_corr.inputs.dimension = 't'
	merge_corr.inputs.output_type = 'NIFTI_GZ'
        merge_corr.inputs.merged_file = 'cope_corr.nii.gz'
	merge_corr.inputs.tr = 2.00
	cope_merge_wf.connect(dsource, 'corr_cope', merge_corr, 'in_files')

        #### dsource (corr_cope) ----> merge_corr (in_files)


        # Node to merge INCORRECT trials across all 6 runs
	merge_incorr = Node(Merge(), 
                            name = 'merge_incorr')
	merge_incorr.inputs.dimension = 't'
	merge_incorr.inputs.output_type = 'NIFTI_GZ'
        merge_incorr.inputs.merged_file = 'cope_incorr.nii.gz'
	merge_incorr.inputs.tr = 2.00
	cope_merge_wf.connect(dsource, 'incorr_cope', merge_incorr, 'in_files')
 
         #### dsource (incorr_cope) ----> merge_incorr (in_files)


        # Node to output data
	dsink = Node(DataSink(), 
                     name = 'dsink')
	dsink.inputs.base_directory = sink_directory 
	dsink.inputs.container = subject_id
	cope_merge_wf.connect(merge_corr, 'merged_file', dsink, 'merged.@corr')
	cope_merge_wf.connect(merge_incorr, 'merged_file', dsink, 'merged.@incorr')
	
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

    	wf.config['execution']['crashdump_dir'] = '/scratch/madlab/crash/wmaze_MRthesis/model5_LSS/merge_copes'
    	wf.base_dir = work_dir + '/' + args.subject_id
    	wf.run(plugin = 'LSF', plugin_args = {'bsub_args': '-q PQ_madlab'})
