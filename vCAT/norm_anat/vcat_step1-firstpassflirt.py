#!/usr/bin/env python

import os
from nipype.pipeline.engine import Workflow, Node, MapNode, JoinNode
from nipype.interfaces.utility import IdentityInterface, Function
from nipype.interfaces import fsl
from nipype.interfaces import ants
import nipype.interfaces.utility as util
from glob import glob
from nipype.interfaces.io import DataSink
from mattfeld_utility_workflows.fs_skullstrip_util import create_freesurfer_skullstrip_workflow
import numpy as np

def pickfirst(func):
    if isinstance(func, list):
        return func[0]
    else:
        return func

def get_subs(subject_id):
    subs = []
    subs.append(('_subject_id_%s/' %subject_id, ''))
    return subs



fs_projdir = '/home/data/madlab/Mattfeld_vCAT/derivatives/freesurfer'
projdir = '/home/data/madlab/Mattfeld_vCAT'
workdir = '/scratch/madlab/crash/mandy/vcat/norm_anat'
sinkdir = '/home/data/madlab/Mattfeld_vCAT/derivatives'

workingdir = os.path.join(workdir,'flirt_fp') #working directory
if not os.path.exists(workingdir):
    os.makedirs(workingdir)


fs_skullstrip_wf = create_freesurfer_skullstrip_workflow()
fs_skullstrip_wf.inputs.inputspec.subjects_dir = fs_projdir


sids = ['vCAT_005', 'vCAT_006', 'vCAT_007', 'vCAT_008', 'vCAT_010', 'vCAT_012', 'vCAT_013', 'vCAT_014', 'vCAT_015', 'vCAT_016',  
	'vCAT_018', 'vCAT_019', 'vCAT_020', 'vCAT_021', 'vCAT_022', 'vCAT_023', 'vCAT_024', 'vCAT_025', 'vCAT_026', 'vCAT_027',  
	'vCAT_028', 'vCAT_029', 'vCAT_030', 'vCAT_031', 'vCAT_032']


#set up FreeSurfer skull-strip workflow
flirt_fp_wf = Workflow(name = 'flirt_fp_wf')
flirt_fp_wf.base_dir = workingdir


#create subject iterable node
subjID_infosource = Node(IdentityInterface(fields = ['subject_id','subjects_dir']), 
                         name = 'subjID_infosource')
subjID_infosource.iterables = ('subject_id', sids)


#connect subject id from subject iterable into skullstripping workflow
flirt_fp_wf.connect(subjID_infosource, 'subject_id', fs_skullstrip_wf, 'inputspec.subject_id')


#create function node to rename output files
getsubs = Node(util.Function(input_names = ['subject_id'],
                             output_names = ['subs'],
                             function = get_subs),
               name = 'getsubs')
getsubs.inputs.ignore_exception = False
flirt_fp_wf.connect(subjID_infosource, 'subject_id', getsubs, 'subject_id')


#create FLIRT node to rigid body transform (6 DOF) skullstripped brains to MNI template
firstpass_flirt = Node(fsl.FLIRT(),
                       name = 'firstpass_flirt')
firstpass_flirt.inputs.reference = '/home/data/madlab/data/mri/templates/MNI/OASIS-30_Atropos_template_in_MNI152.nii.gz'
firstpass_flirt.inputs.dof = 6
firstpass_flirt.inputs.coarse_search = 5
firstpass_flirt.inputs.fine_search = 40
flirt_fp_wf.connect(fs_skullstrip_wf, 'outputspec.skullstripped_file', firstpass_flirt, 'in_file')


#provide files with unique output name
name_unique = Node(util.Rename(format_string = 'vcat_skullstrip_%(subject_id)s'),
                   name = 'name_unique')
name_unique.inputs.keep_ext = True 
flirt_fp_wf.connect(subjID_infosource, 'subject_id', name_unique, 'subject_id')
flirt_fp_wf.connect(firstpass_flirt, 'out_file', name_unique, 'in_file')


#move results to designated results folder
datasink = Node(DataSink(), 
                name = "datasink")
datasink.inputs.base_directory = os.path.join(sinkdir, 'firstpass_flirt')
flirt_fp_wf.connect(getsubs, 'subs', datasink, 'substitutions')
flirt_fp_wf.connect(name_unique, 'out_file', datasink, 'skullstrip_flirt1')


#run workflow
flirt_fp_wf.run(plugin = 'SLURM', plugin_args={'sbatch_args': ('--partition IB_40C_512G --qos pq_madlab --account iacc_madlab -t 24:00:00 -N 1 -n 1'), 'overwrite':True})

