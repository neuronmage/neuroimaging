#!/usr/bin/env python
import os
from glob import glob
from nipype import Node, Function, Workflow, IdentityInterface
from nipype.interfaces.freesurfer import ReconAll
from nipype.interfaces.io import DataGrabber

#PROJECT DATA DIRECTORY
data_dir = '/home/data/madlab/Pruden_SEA/dset/'

for subject in glob('/home/data/madlab/Pruden_SEA/dset/sub-*'):
    if 'sub-' + subject[-4:] in os.listdir('/home/data/madlab/Pruden_SEA/derivatives/freesurfer_SEA/'):
        continue  
    sid = subject[-4:]

    info = dict(T1 = [['subject_id']])

    infosource = Node(IdentityInterface(fields = ['subject_id']), 
                      name = 'infosource')
    infosource.iterables = ('subject_id', sid)


    #datasource node to get the T1 files
    datasource = Node(DataGrabber(infields = ['subject_id'],
                                  outfields = list(info.keys())),
                      name = 'datasource')
    datasource.inputs.template = '%s/%s'
    datasource.inputs.base_directory = os.path.abspath(data_dir)
    datasource.inputs.field_template = dict(T1 = 'sub-%s/anat/sub-%s_T1w.nii.gz')
    datasource.inputs.template_args = info
    datasource.inputs.sort_filelist = True


    reconall_node = Node(ReconAll(), 
                         name = 'reconall_node')
    reconall_node.inputs.openmp = 4 #multiple processors
    reconall_node.inputs.hippocampal_subfields_T1 = True 
    reconall_node.inputs.subjects_dir = '/home/data/madlab/Pruden_SEA/derivatives/freesurfer/' #FS surface save location
    reconall_node.plugin_args={'sbatch_args': ('-p IB_44C_512G --account iacc_madlab --qos pq_madlab -n 4'), 'overwrite': True}


    wf = Workflow(name = 'fsrecon')
    wf.connect(infosource, 'subject_id', datasource, 'subject_id')
    wf.connect(infosource, 'subject_id', reconall_node, 'subject_id')
    wf.connect(datasource, 'T1', reconall_node, 'T1_files')

    wf.base_dir = os.path.abspath('/scratch/madlab/Pruden_SEA/freesurfer/')
    wf.run(plugin='SLURM', plugin_args={'sbatch_args': ('-p investor --account iacc_madlab --qos pq_madlab -n 1')})
