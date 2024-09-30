#!/usr/bin/env python

import os
from nipype.pipeline.engine import Workflow
from nipype.pipeline.engine import Node
from nipype.pipeline.engine import MapNode
from nipype.pipeline.engine import JoinNode
from nipype.interfaces.utility import IdentityInterface
from nipype.interfaces.utility import Function
from nipype.interfaces import fsl
from nipype.interfaces import ants
from glob import glob
from nipype.interfaces.io import DataSink
from nipype.interfaces.ants.legacy import buildtemplateparallel
from mattfeld_utility_workflows.fs_skullstrip_util import create_freesurfer_skullstrip_workflow

fs_projdir = '/home/data/madlab/surfaces/wmaze'
projdir = '/home/data/madlab/data/mri/wmaze'
workdir = '/scratch/madlab'
workingdir = os.path.join(workdir,'ants_template_workdir2') #working directory
if not os.path.exists(workingdir):
    os.makedirs(workingdir)

fs_skullstrip_wf = create_freesurfer_skullstrip_workflow()
fs_skullstrip_wf.inputs.inputspec.subjects_dir = fs_projdir

sids = ['WMAZE_001', 'WMAZE_002', 'WMAZE_003', 'WMAZE_004', 'WMAZE_005', 'WMAZE_006', 'WMAZE_007', 'WMAZE_008',
        'WMAZE_010', 'WMAZE_012', 'WMAZE_017', 'WMAZE_018', 'WMAZE_019', 'WMAZE_020', 'WMAZE_021', 'WMAZE_022',
        'WMAZE_023', 'WMAZE_024', 'WMAZE_026', 'WMAZE_027']

# Set up the FreeSurfer skull stripper work flow
ants_buildtemplate_wf = Workflow(name='ants_buildtemplate_wf')
ants_buildtemplate_wf.base_dir = workingdir

subjID_infosource = Node(IdentityInterface(fields=['subject_id','subjects_dir']), name = 'subjID_infosource')
#subjID_infosource.inputs.subject_ids = sids
subjID_infosource.iterables = ('subject_id', sids)

ants_buildtemplate_wf.connect(subjID_infosource, 'subject_id', fs_skullstrip_wf, 'inputspec.subject_id')

# Use a JoinNode to aggregrate all of the outputs from the fs_skullstrip_wf
skullstripped_images = JoinNode(IdentityInterface(fields=['brainonly_images']),
                             joinsource='subjID_infosource',
                             joinfield='brainonly_images', name='skullstripped_images')
ants_buildtemplate_wf.connect(fs_skullstrip_wf, 'outputspec.skullstripped_file', skullstripped_images, 'brainonly_images')

# Create a FLIRT node to rigid body transform (6 DOF) skullstripped brains to a MNI template
firstpass_flirt = MapNode(fsl.FLIRT(),
                          iterfield=['in_file'],
                          name='firstpass_flirt')
firstpass_flirt.inputs.reference = "/home/data/madlab/data/mri/templates/MNI/OASIS-30_Atropos_template_in_MNI152.nii.gz"
firstpass_flirt.inputs.dof = 6
ants_buildtemplate_wf.connect(skullstripped_images, 'brainonly_images', firstpass_flirt, 'in_file')

ants_template = Node(buildtemplateparallel(), name = 'ants_template')
ants_template.inputs.parallelization = 0
ants_buildtemplate_wf.connect(firstpass_flirt, 'out_file', ants_template, 'in_files')

# Move the results to a designated results folder
datasink = Node(DataSink(), name="datasink")
datasink.inputs.base_directory = os.path.join(projdir, "wmaze_T1_template2")

ants_buildtemplate_wf.connect(ants_template, 'final_template_file', datasink, 'PrimaryTemplate')

# Run the workflow
ants_buildtemplate_wf.run(plugin='LSF', plugin_args={'bsub_args' : ('-q PQ_madlab')})




