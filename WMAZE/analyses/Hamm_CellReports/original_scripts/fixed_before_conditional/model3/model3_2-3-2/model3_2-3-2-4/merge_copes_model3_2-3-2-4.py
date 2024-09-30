#!/usr/bin/env python

"""
============================================================================================
wmaze_fMRI: Mandy Thesis -- Fixed before Condition - Merge COPES for beta series correlation
============================================================================================
Very simple Nipype script to merge beta series copes from wmaze_lvl1_model3_2-3-2-4.py

- WMAZE Model 3 Version 2.3.2.4
2  - LSS model -- Beta series correlation analysis
3  - Remove last three trials/volumes
2  - Use FSL ROI to recreate EPI data, removing last 3 volumes
     - Removed last 3 trials before EV creation
4  - Model 3 - 11/11 window


- python merge_copes_model3_2-3-2.py -s WMAZE_001
                                     -o /home/data/madlab/data/mri/wmaze/frstlvl/wmaze_MRthesis/fixed_before_conditional/model3_2-3-2-4/merge_copes
                                     -w /home/data/madlab/scripts/wmaze/anal_MR_thesis/fixed_before_conditional/model3/model3_2-3-2-4/status/merge_copes
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
    	info = dict(during11_cope = [['subject_id']],
		    after11_cope = [['subject_id']])

        # Node to grab corr and incorr cope files
	dsource = Node(DataGrabber(infields = ['subject_id'], 
                                   outfields = info.keys()), 
                       name = 'dsource')
	dsource.inputs.template = '*'
	dsource.inputs.subject_id = subject_id
	dsource.inputs.base_directory = os.path.abspath('/home/data/madlab/data/mri/wmaze/frstlvl/wmaze_MRthesis/fixed_before_conditional/model3_2-3-2-4')
	dsource.inputs.field_template = dict(during11_cope = '11-11window/%s/learning_window/*.nii.gz',
				             after11_cope = '11-11window/%s/after_window/*.nii.gz')
	dsource.inputs.template_args = info
	dsource.inputs.sort_filelist = True
	dsource.inputs.ignore_exception = False
	dsource.inputs.raise_on_empty = True


	# Node to merge during11 trials across all 6 runs
	merge_during11 = Node(Merge(), 
                              name = 'merge_during11')
	merge_during11.inputs.dimension = 't'
	merge_during11.inputs.output_type = 'NIFTI_GZ'
        merge_during11.inputs.merged_file = 'cope_during11.nii.gz'
	merge_during11.inputs.tr = 2.00
	cope_merge_wf.connect(dsource, 'during11_cope', merge_during11, 'in_files')

        #### dsource (during11_cope) ----> merge_during11 (in_files)


        # Node to merge after11 trials across all 6 runs
	merge_after11 = Node(Merge(), 
                             name = 'merge_after11')
	merge_after11.inputs.dimension = 't'
	merge_after11.inputs.output_type = 'NIFTI_GZ'
        merge_after11.inputs.merged_file = 'cope_after11.nii.gz'
	merge_after11.inputs.tr = 2.00
	cope_merge_wf.connect(dsource, 'after11_cope', merge_after11, 'in_files')
 
         #### dsource (after11_cope) ----> merge_after11 (in_files)


        # Node to output data
	dsink = Node(DataSink(), 
                     name = 'dsink')
	dsink.inputs.base_directory = sink_directory 
	dsink.inputs.container = subject_id
	cope_merge_wf.connect(merge_during11, 'merged_file', dsink, 'merged.@during11')
	cope_merge_wf.connect(merge_after11, 'merged_file', dsink, 'merged.@after11')
	
        #### merge_during11 (merged_file) ----> dsink (during11)
        #### merge_after11 (merged_file) ----> dsink (after11)

	return cope_merge_wf

	
"""
Creates the full workflow
"""


def create_workflow(args, 
                    name = 'model3_2-3-2-4_merge_copes'):
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

    	wf.config['execution']['crashdump_dir'] = '/scratch/madlab/crash/wmaze_MRthesis/model3_2-3-2-4/merge_copes'
    	wf.base_dir = work_dir + '/' + args.subject_id
    	wf.run(plugin = 'LSF', plugin_args = {'bsub_args': '-q PQ_madlab'})
