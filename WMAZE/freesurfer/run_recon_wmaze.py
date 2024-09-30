#!/usr/bin/env python

#SBATCH -p investor
#SBATCH --qos pq_madlab
#SBATCH -o /scratch/madlab/surfaces/wmaze/run_recon_out
#SBATCH -e /scratch/madlab/surfaces/wmaze/run_recon_err

import os
from glob import glob

from nipype import Node, Function, Workflow, IdentityInterface
from nipype.interfaces.freesurfer import ReconAll
from nipype.interfaces.io import DataGrabber

# CURRENT PROJECT DATA DIRECTORY
data_dir = '/home/data/madlab/data/mri/wmaze/'

# CURRENT PROJECT SUBJECT IDS
sids = ['template']

'''sids = ['WMAZE_001', 'WMAZE_002', 'WMAZE_004', 'WMAZE_005', 'WMAZE_006', 'WMAZE_007',
        'WMAZE_008', 'WMAZE_009', 'WMAZE_010', 'WMAZE_012', 'WMAZE_017', 'WMAZE_018',
        'WMAZE_019', 'WMAZE_020', 'WMAZE_021', 'WMAZE_022', 'WMAZE_023', 'WMAZE_024',
        'WMAZE_026', 'WMAZE_027']'''

info = dict(T1 = [['subject_id']])

infosource = Node(IdentityInterface(fields = ['subject_id']), 
                  name = 'infosource')
infosource.iterables = ('subject_id', sids)

# Create a datasource node to get the T1 file
datasource = Node(DataGrabber(infields = ['subject_id'],
                              outfields = info.keys()),
                  name = 'datasource')
datasource.inputs.template = '%s/%s'
datasource.inputs.base_directory = os.path.abspath(data_dir)
datasource.inputs.field_template = dict(T1 = 'wmaze_T1_%s/T_wmaze_template.nii.gz')
datasource.inputs.template_args = info
datasource.inputs.sort_filelist = True

reconall_node = Node(ReconAll(), 
                     name = 'reconall_node')
reconall_node.inputs.openmp = 2
reconall_node.inputs.args = '-hippocampal-subfields-T1'
reconall_node.inputs.subjects_dir = '/home/data/madlab/surfaces/wmaze'
reconall_node.inputs.terminal_output = 'allatonce'
reconall_node.plugin_args={'sbatch_args': ('-p investor --qos pq_madlab -N 1 -n 2'), 'overwrite': True}

wf = Workflow(name = 'fsrecon')

wf.connect(infosource, 'subject_id', datasource, 'subject_id')
wf.connect(infosource, 'subject_id', reconall_node, 'subject_id')
wf.connect(datasource, 'T1', reconall_node, 'T1_files')

wf.base_dir = os.path.abspath('/scratch/madlab/wmaze/')
#wf.config['execution']['job_finished_timeout'] = 65

wf.run(plugin='SLURM', plugin_args={'sbatch_args': ('-p investor --qos pq_madlab -N 1 -n 1')})
