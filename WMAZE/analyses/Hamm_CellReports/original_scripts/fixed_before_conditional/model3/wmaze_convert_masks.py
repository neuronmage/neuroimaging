#!/usr/bin/env python

from nipype.pipeline.engine import Workflow, Node, MapNode
from nipype.interfaces.utility import IdentityInterface, Merge
from nipype.interfaces.io import DataGrabber, DataSink
from nipype.interfaces import ants
from mattfeld_utility_workflows.fs_skullstrip_util import create_freesurfer_skullstrip_workflow
from nipype.interfaces.c3 import C3dAffineTool
from nipype.interfaces.utility import Function

###############
## Functions ##
###############

get_len = lambda x: len(x)

# Project base directory
proj_dir = '/home/data/madlab/data/mri/wmaze'
# Where to find the surfaces
fs_projdir = '/home/data/madlab/surfaces/wmaze'
# Where to find the crash files
work_dir = '/scratch/madlab/crash/wmaze_MRthesis/fixed_before_conditional/model3/model3_1-3-2/volume_transform'
# Where the output will go
sink_dir = '/home/data/madlab/scripts/wmaze/anal_MR_thesis/fixed_before_conditional/model3/masks'

sids = ['WMAZE_020']


# Define project workflow
wf = Workflow("wf")
wf.base_dir = work_dir

# Node to iterate through the subjects
subj_iterable = Node(IdentityInterface(fields = ['subject_id'], 
                                       mandatory_inputs = True), 
                     name = 'subj_interable')
subj_iterable.iterables = ('subject_id', sids)


# Variable containing the dictionary keys for the datasource node
info = dict(mpfc_mask = [['subject_id']],
            dlpfc_mask = [['subject_id']],
            hpc_mask = [['subject_id']],
            caudate_mask = [['subject_id']],
            ants_warp = [['subject_id', 'subject_id', 'output']])

# Node to grab the various data files for each subject
datasource = Node(DataGrabber(infields = ['subject_id'],
                              outfields = info.keys()),
                  name = "datasource")
datasource.inputs.base_directory = proj_dir
datasource.inputs.field_template = dict(mpfc_mask = 'roi_analysis/mask/anat_masks/_subject_id_%s/_anatmask_xfm*/lh-mPFC_rac-cac_warped.nii.gz',
                                        dlpfc_mask = 'roi_analysis/mask/anat_masks/_subject_id_%s/_anatmask_xfm0/*h-dlPFC_lausanne_warped.nii.gz',
                                        hpc_mask = 'roi_analysis/mask/anat_masks/_subject_id_%s/_anatmask_xfm0/*h-hippocampus_warped.nii.gz',
                                        caudate_mask = 'roi_analysis/mask/anat_masks/_subject_id_%s/_anatmask_xfm0/*h_caudate_anat_mask_warped.nii.gz',
                                        ants_warp = 'norm_anat/%s/anat2targ_xfm/_subject_id_%s/%s*.h5')
datasource.inputs.ignore_exception = False
datasource.inputs.raise_on_empty = True
datasource.inputs.sort_filelist = True
datasource.inputs.subject_id = sids
datasource.inputs.template = '*'
datasource.inputs.template_args = info
wf.connect(subj_iterable, "subject_id", datasource, "subject_id")



mpfc2targ = Node(ants.ApplyTransforms(), 
                 name = 'mpfc2targ')
mpfc2targ.inputs.interpolation = 'NearestNeighbor'
mpfc2targ.inputs.invert_transform_flags = [False]
mpfc2targ.inputs.terminal_output = 'file'
mpfc2targ.inputs.args = '--float'
mpfc2targ.inputs.dimension = 3
# Specify the image whose space you are converting into
mpfc2targ.inputs.reference_image = '/home/data/madlab/data/mri/wmaze/wmaze_T1_template/T_wmaze_template.nii.gz'
wf.connect(datasource, 'mpfc_mask', mpfc2targ, 'input_image')
wf.connect(datasource, 'ants_warp', mpfc2targ, 'transforms')


# MapNode: Warp copes to target
dlpfc2targ = MapNode(ants.ApplyTransforms(), 
                     iterfield = ['input_image'], 
                     name = 'dlpfc2targ')
# Method of interpolation used in conversion
dlpfc2targ.inputs.interpolation = 'NearestNeighbor'
dlpfc2targ.inputs.invert_transform_flags = [False]
dlpfc2targ.inputs.terminal_output = 'file'
#dlpfc2targ.inputs.float = True
dlpfc2targ.inputs.args = '--float'
dlpfc2targ.inputs.dimension = 3
# Specify the image whose space you are converting into
dlpfc2targ.inputs.reference_image = '/home/data/madlab/data/mri/wmaze/wmaze_T1_template/T_wmaze_template.nii.gz'
wf.connect(datasource, 'dlpfc_mask', dlpfc2targ, 'input_image')
wf.connect(datasource, 'ants_warp', dlpfc2targ, 'transforms')


# MapNode: Warp copes to target
hpc2targ = MapNode(ants.ApplyTransforms(), 
                     iterfield = ['input_image'], 
                     name = 'hpc2targ')
# Method of interpolation used in conversion
hpc2targ.inputs.interpolation = 'NearestNeighbor'
hpc2targ.inputs.invert_transform_flags = [False]
hpc2targ.inputs.terminal_output = 'file'
hpc2targ.inputs.args = '--float'
hpc2targ.inputs.dimension = 3
# Specify the image whose space you are converting into
hpc2targ.inputs.reference_image = '/home/data/madlab/data/mri/wmaze/wmaze_T1_template/T_wmaze_template.nii.gz'
wf.connect(datasource, 'hpc_mask', hpc2targ, 'input_image')
wf.connect(datasource, 'ants_warp', hpc2targ, 'transforms')


# MapNode: Warp copes to target
caudate2targ = MapNode(ants.ApplyTransforms(), 
                     iterfield = ['input_image'], 
                     name = 'caudate2targ')
# Method of interpolation used in conversion
caudate2targ.inputs.interpolation = 'NearestNeighbor'
caudate2targ.inputs.invert_transform_flags = [False]
caudate2targ.inputs.terminal_output = 'file'
caudate2targ.inputs.args = '--float'
caudate2targ.inputs.dimension = 3
# Specify the image whose space you are converting into
caudate2targ.inputs.reference_image = '/home/data/madlab/data/mri/wmaze/wmaze_T1_template/T_wmaze_template.nii.gz'
wf.connect(datasource, 'caudate_mask', caudate2targ, 'input_image')
wf.connect(datasource, 'ants_warp', caudate2targ, 'transforms')


# Node: group.sinker
stats_sinker = Node(DataSink(infields = None), 
                         name = 'stats_sinker')
stats_sinker.inputs._outputs = {}
stats_sinker.inputs.base_directory = sink_dir
stats_sinker.inputs.ignore_exception = False
stats_sinker.inputs.parameterization = True
stats_sinker.inputs.remove_dest_dir = False
wf.connect(subj_iterable, 'subject_id', stats_sinker, 'container')
wf.connect(mpfc2targ, 'output_image', stats_sinker, 'mpfc_trans')
wf.connect(dlpfc2targ, 'output_image', stats_sinker, 'dlpfc_trans')
wf.connect(hpc2targ, 'output_image', stats_sinker, 'hpc_trans')
wf.connect(caudate2targ, 'output_image', stats_sinker, 'caudate_trans')


wf.config['execution']['crashdump_dir'] = work_dir
wf.run(plugin = 'LSF', plugin_args = {'bsub_args': '-q PQ_madlab'})


wf.config['execution']['crashdump_dir'] = work_dir
wf.run(plugin = 'LSF', plugin_args = {'bsub_args': '-q PQ_madlab'})
