#!/usr/bin/env python

"""
============================================================================================
wmaze_fMRI: Mandy Thesis -- Fixed before Condition - Merge COPES for beta series correlation
============================================================================================
Very simple Nipype script to merge beta series copes from wmaze_lvl1_model3_2-3-2-5.py

- WMAZE Model 3 Version 2.3.2.5
2  - LSS model -- Beta series correlation analysis
3  - Remove last three trials/volumes
2  - Use FSL ROI to recreate EPI data, removing last 3 volumes
     - Removed last 3 trials before EV creation
5  - Model 3 - Derivative learning


- python merge_copes_model3_2-3-2-5.py -s WMAZE_001
                                       -o /home/data/madlab/data/mri/wmaze/frstlvl/wmaze_MRthesis/fixed_before_conditional/model3_2-3-2-5/merge_copes
                                       -w /home/data/madlab/scripts/wmaze/anal_MR_thesis/fixed_before_conditional/model3/model3_2-3-2-5/status/merge_copes
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
    	info = dict(learning_cope = [['subject_id']],
		    nonlearning_cope = [['subject_id']])

        # Node to grab corr and incorr cope files
	dsource = Node(DataGrabber(infields = ['subject_id'], 
                                   outfields = info.keys()), 
                       name = 'dsource')
	dsource.inputs.template = '*'
	dsource.inputs.subject_id = subject_id
	dsource.inputs.base_directory = os.path.abspath('/home/data/madlab/data/mri/wmaze/frstlvl/wmaze_MRthesis/fixed_before_conditional/model3_2-3-2-5')
	dsource.inputs.field_template = dict(learning_cope = '%s/learning/*.nii.gz',
				             nonlearning_cope = '%s/nonlearning/*.nii.gz')
	dsource.inputs.template_args = info
	dsource.inputs.sort_filelist = True
	dsource.inputs.ignore_exception = False
	dsource.inputs.raise_on_empty = True


	# Node to merge learning trials across all 6 runs
	merge_learning = Node(Merge(), 
                              name = 'merge_learning')
	merge_learning.inputs.dimension = 't'
	merge_learning.inputs.output_type = 'NIFTI_GZ'
        merge_learning.inputs.merged_file = 'cope_learning.nii.gz'
	merge_learning.inputs.tr = 2.00
	cope_merge_wf.connect(dsource, 'learning_cope', merge_learning, 'in_files')

        #### dsource (learning_cope) ----> merge_learning (in_files)


        # Node to merge nonlearning trials across all 6 runs
	merge_nonlearning = Node(Merge(), 
                             name = 'merge_nonlearning')
	merge_nonlearning.inputs.dimension = 't'
	merge_nonlearning.inputs.output_type = 'NIFTI_GZ'
        merge_nonlearning.inputs.merged_file = 'cope_nonlearning.nii.gz'
	merge_nonlearning.inputs.tr = 2.00
	cope_merge_wf.connect(dsource, 'nonlearning_cope', merge_nonlearning, 'in_files')
 
         #### dsource (nonlearning_cope) ----> merge_nonlearning (in_files)


        # Node to output data
	dsink = Node(DataSink(), 
                     name = 'dsink')
	dsink.inputs.base_directory = sink_directory 
	dsink.inputs.container = subject_id
	cope_merge_wf.connect(merge_learning, 'merged_file', dsink, 'merged.@learning')
	cope_merge_wf.connect(merge_nonlearning, 'merged_file', dsink, 'merged.@nonlearning')
	
        #### merge_learning (merged_file) ----> dsink (learning)
        #### merge_nonlearning (merged_file) ----> dsink (nonlearning)

	return cope_merge_wf

	
"""
Creates the full workflow
"""


def create_workflow(args, 
                    name = 'model3_2-3-2-5_merge_copes'):
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

    	wf.config['execution']['crashdump_dir'] = '/scratch/madlab/crash/wmaze_MRthesis/model3_2-3-2-5/merge_copes'
    	wf.base_dir = work_dir + '/' + args.subject_id
    	wf.run(plugin = 'LSF', plugin_args = {'bsub_args': '-q PQ_madlab'})
