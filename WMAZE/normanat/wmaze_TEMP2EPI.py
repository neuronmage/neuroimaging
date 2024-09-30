#!/usr/bin/env python

#SBATCH -p investor
#SBATCH --qos pq_madlab
#SBATCH -o /scratch/madlab/crash/temp2epi_out
#SBATCH -e /scratch/madlab/crash/temp2epi_err

"""
================================================================
wmaze_STRUCT: ANTs
================================================================


 
"""
from nipype.pipeline.engine import Workflow, Node, MapNode
from nipype.interfaces.utility import IdentityInterface, Merge
from nipype.interfaces.io import DataGrabber
from nipype.interfaces.io import DataSink
from nipype.interfaces import ants
from mattfeld_utility_workflows.fs_skullstrip_util import create_freesurfer_skullstrip_workflow
from nipype.interfaces.c3 import C3dAffineTool
from nipype.interfaces.utility import Function
import os
import nipype.interfaces.io as nio
import nipype.pipeline.engine as pe

# Functions
get_len = lambda x: len(x)

proj_dir = '/home/data/madlab/data/mri/wmaze'
fs_projdir = '/home/data/madlab/surfaces/wmaze'
work_dir = '/scratch/madlab/wmaze/temp2epi'
sink_dir = '/home/data/madlab/data/mri/wmaze/masks/temp2epi'

sids = ['WMAZE_023']
#sids = os.listdir(proj_dir + '/norm_stats/MRthesis_norm_stats/fixed_before_conditional/model3_1-3-2')[1:]

# Workflow
temp2epi_wf = Workflow("temp2epi_wf")
temp2epi_wf.base_dir = work_dir

# Node: subject_iterable
subj_iterable = pe.Node(IdentityInterface(fields=['subject_id'], mandatory_inputs=True), name='subj_interable')
subj_iterable.iterables = ('subject_id', sids)

# WORKFLOW: freesurfer skull stripping workflow
fs_skullstrip_wf = create_freesurfer_skullstrip_workflow()
fs_skullstrip_wf.inputs.inputspec.subjects_dir = fs_projdir
temp2epi_wf.connect(subj_iterable, 'subject_id', fs_skullstrip_wf, 'inputspec.subject_id')

# Node: datasource
datasource_norm = pe.Node(DataGrabber(infields=['subject_id'],
                                      outfields=['bbreg_xfm', 'ants_warp', 'mean_image']),
                          name="datasource_norm")
datasource_norm.inputs.base_directory = os.path.abspath(proj_dir)
datasource_norm.inputs.field_template = dict(bbreg_xfm='preproc/%s/bbreg/_fs_register0/%s*.mat',
                                             ants_warp='norm_anat/%s/targ2anat_xfm/_subject_id_%s/%s*.h5',
                                             mean_image='preproc/%s/ref/%s*.nii.gz')
datasource_norm.inputs.ignore_exception = False
datasource_norm.inputs.raise_on_empty = True
datasource_norm.inputs.sort_filelist = True
datasource_norm.inputs.template = '*'
datasource_norm.inputs.template_args = dict(bbreg_xfm=[['subject_id', 'wmaze']],
                                            ants_warp =[['subject_id', 'subject_id', 'output']],
                                            mean_image=[['subject_id', 'wmaze']])
temp2epi_wf.connect(subj_iterable, "subject_id", datasource_norm, "subject_id")

# Node: convert2itk INVERSE
convert2itk = pe.Node(C3dAffineTool(), name='convert2itk')
convert2itk.inputs.fsl2ras = True
convert2itk.inputs.itk_transform = True
convert2itk.inputs.args = '-inv'
temp2epi_wf.connect(datasource_norm, 'bbreg_xfm', convert2itk, 'transform_file')
temp2epi_wf.connect(datasource_norm, 'mean_image', convert2itk, 'source_file')
temp2epi_wf.connect(fs_skullstrip_wf, 'outputspec.skullstripped_file', convert2itk, 'reference_file')

# Node: concatenate the affine and ants transforms into a list
merge_xfm = pe.Node(Merge(2), iterfield = ['in2'], name='merge_xfm')
temp2epi_wf.connect(convert2itk, 'itk_transform', merge_xfm, 'in1')
temp2epi_wf.connect(datasource_norm, 'ants_warp', merge_xfm, 'in2')

# MapNode: Warp target to EPI
targ2EPI = pe.MapNode(ants.ApplyTransforms(), iterfield=['input_image'], name='targ2EPI')
targ2EPI.inputs.input_image_type = 3
targ2EPI.inputs.interpolation = 'NearestNeighbor'
targ2EPI.inputs.invert_transform_flags = [False, False]
targ2EPI.inputs.args = '--float'
targ2EPI.inputs.num_threads = 4
targ2EPI.plugin_args = {'sbatch_args': '--ntasks=1 --cpus-per-task=4'}
targ2EPI.inputs.input_image = proj_dir + '/grplvl/grplvl_MRthesis_randomise/model3_1-3-2/all_corr_minus_all_incorr/output/corrected/_contrast_all_corr_minus_all_incorr/mask_gt_0_99_bin.nii.gz'
temp2epi_wf.connect(datasource_norm, 'mean_image', targ2EPI, 'reference_image')
temp2epi_wf.connect(merge_xfm, 'out', targ2EPI, 'transforms')

# place to save outputs from transformation
temp2epi_sinker = pe.Node(nio.DataSink(), name="temp2epi_sinker") 
temp2epi_sinker.inputs.base_directory = sink_dir
temp2epi_wf.connect(subj_iterable, "subject_id", temp2epi_sinker, "container")
temp2epi_wf.connect(targ2EPI, "output_image", temp2epi_sinker, "temp2epi")

temp2epi_wf.config['execution']['crashdump_dir'] = '/scratch/madlab/crash/wmaze_temp2epi'
temp2epi_wf.run(plugin='SLURM', plugin_args={'sbatch_args': ('-p investor --qos pq_madlab -N 1 -n 1'), 'overwrite': True})

