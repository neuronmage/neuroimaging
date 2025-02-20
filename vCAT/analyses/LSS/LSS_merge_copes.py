#!/usr/bin/env python

"""
=============================================================================
LSS_fMRI-- Fixed before Conditional - Merge COPES for beta series correlation
=============================================================================
Very simple Nipype script to merge beta series copes from LSS_lvl1.py

Model LSS
2  - LSS model -- Beta series correlation analysis
3  - Remove last three trials/volumes
2  - Use FSL ROI to recreate EPI data, removing last 3 volumes
     - Removed last 3 trials before EV creation
5  - Model 3 - Derivative learning


- python LSS_merge_copes.py -s WMAZE_001
                            -o /home/data/madlab/data/mri/wmaze/frstlvl/model_LSS/merge_copes
                            -w /home/data/madlab/scripts/wmaze/model_LSS/status/merge_copes
"""

import os
from nipype.pipeline.engine import Workflow, Node, MapNode
from nipype.interfaces.utility import Function
from nipype.interfaces.io import DataGrabber
from nipype.interfaces.io import DataSink
from nipype.interfaces.fsl import Merge

# Function to merge copes for each participant
def cope_merge_wf(subject_id, sink_directory,
	          name = 'cope_merge_wf'):   
	cope_merge_wf = Workflow(name = 'cope_merge_wf')
        
    	info = dict(learning_cope = [['subject_id']], #dictionary for Datagrabber
		    nonlearning_cope = [['subject_id']])

        #node to grab corr and incorr cope files
	datasource = Node(DataGrabber(infields = ['subject_id'], outfields = info.keys()), 
                          name = 'datasource')
	datasource.inputs.template = '*'
	datasource.inputs.subject_id = subject_id
	datasource.inputs.base_directory = os.path.abspath('/home/data/madlab/data/mri/wmaze/frstlvl/model_LSS2')
	datasource.inputs.field_template = dict(learning_cope = '%s/deriv/learn/*.nii.gz',
				                nonlearning_cope = '%s/deriv/nonlearn/*.nii.gz')
	datasource.inputs.template_args = info
	datasource.inputs.sort_filelist = True
	datasource.inputs.ignore_exception = False
	datasource.inputs.raise_on_empty = True


	#node to merge learning trials across all 6 runs
	merge_learning = Node(Merge(), name = 'merge_learning')
	merge_learning.inputs.dimension = 't'
	merge_learning.inputs.output_type = 'NIFTI_GZ'
        merge_learning.inputs.merged_file = 'cope_learning.nii.gz'
	merge_learning.inputs.tr = 2.00
	cope_merge_wf.connect(datasource, 'learning_cope', merge_learning, 'in_files')


        #node to merge nonlearning trials across all 6 runs
	merge_nonlearning = Node(Merge(), name = 'merge_nonlearning')
	merge_nonlearning.inputs.dimension = 't'
	merge_nonlearning.inputs.output_type = 'NIFTI_GZ'
        merge_nonlearning.inputs.merged_file = 'cope_nonlearning.nii.gz'
	merge_nonlearning.inputs.tr = 2.00
	cope_merge_wf.connect(datasource, 'nonlearning_cope', merge_nonlearning, 'in_files')


        #node to output data
	dsink = Node(DataSink(), name = 'dsink')
	dsink.inputs.base_directory = sink_directory 
	dsink.inputs.container = subject_id
	cope_merge_wf.connect(merge_learning, 'merged_file', dsink, 'merged.@learning')
	cope_merge_wf.connect(merge_nonlearning, 'merged_file', dsink, 'merged.@nonlearning')
	
	return cope_merge_wf

	
###############################
## Creates the full workflow ##
###############################

def create_workflow(args, name = 'model_LSS_merge_copes'):
	kwargs = dict(subject_id = args.subject_id, sink_directory = os.path.abspath(args.out_dir), name = name)
    	cope_workflow = cope_merge_wf(**kwargs)
    	return cope_workflow

if __name__ == "__main__":
    	from argparse import ArgumentParser
	parser = ArgumentParser(description =__doc__)
    	parser.add_argument("-s", "--subject_id", dest = "subject_id", help = "Current subject id", required = True)
    	parser.add_argument("-o", "--output_dir", dest = "out_dir", help = "Output directory")
    	parser.add_argument("-w", "--work_dir", dest = "work_dir", help = "Working directory")
    	args = parser.parse_args()

    	wf = create_workflow(args)

    	if args.work_dir:
        	work_dir = os.path.abspath(args.work_dir)
    	else:
        	work_dir = os.getcwd()

    	wf.config['execution']['crashdump_dir'] = '/scratch/madlab/crash/mandy_crash/model_LSS2/merge_copes/'
   	wf.base_dir = work_dir + '/' + args.subject_id
	wf.run(plugin='SLURM', plugin_args={'sbatch_args': ('-p investor --qos pq_madlab -N 1 -n 1'), 'overwrite': True})

