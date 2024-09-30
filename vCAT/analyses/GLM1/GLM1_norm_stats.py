#!/usr/bin/env python

"""
==================================================================
GLM1_fMRI: Fixed before Conditional -- Model GLM1
==================================================================
Normal statistics workflow for UM GE 750 wmaze task data.

- vCAT GLM1 
  - Does not include fixed before conditional trials with intervening BLs
  - Use FSL ROI to recreate EPI data

- python GLM1_norm_stats.py -s sub-0??
                         -o /home/data/madlab/Mattfeld_vCAT/derivatives/GLM1/norm_stats/sub-0??
                         -w /scratch/madlab/crash/mandy/vcat/GLM1/norm_stats
"""

from nipype.pipeline.engine import Workflow, Node, MapNode
from nipype.interfaces.utility import IdentityInterface, Merge, Function
from nipype.interfaces.io import DataGrabber, DataSink
from nipype.interfaces import ants
from mattfeld_utility_workflows.fs_skullstrip_util import create_freesurfer_skullstrip_workflow
from nipype.interfaces.c3 import C3dAffineTool


###############
## Functions ##
###############

get_len = lambda x: len(x)


def get_contrasts(data_inputs):
    import os
    infiles = [os.path.split(d[0])[1] for d in data_inputs]
    contrasts = [inf[7:].split('.nii')[0] for inf in infiles]
    print(contrasts)
    return contrasts


def get_substitutions(subject_id, cons):
    subs = [('_subject_id_{0}'.format(subject_id),'')]
    for i, con in enumerate(cons):
        subs.append(('_cope2targ{0}/cope'.format(i), 'cope{0}'.format(con)))
        subs.append(('_varcope2targ{0}/varcope'.format(i), 'varcope{0}'.format(con)))
    return subs


proj_dir = '/home/data/madlab/Mattfeld_vCAT/derivatives'
fs_projdir = '/home/data/madlab/Mattfeld_vCAT/derivatives/freesurfer'
work_dir = '/scratch/madlab/crash/mandy/vcat/GLM1/norm_stats'
sink_dir = '/home/data/madlab/Mattfeld_vCAT/derivatives/GLM1/norm_stats'

sids = ['sub-005', 'sub-006', 'sub-007', 'sub-008', 'sub-010', 
        'sub-012', 'sub-013', 'sub-014', 'sub-015', 'sub-016',  
	'sub-018', 'sub-019', 'sub-020', 'sub-021', 'sub-022', 
        'sub-023', 'sub-024', 'sub-025', 'sub-026', 'sub-027',  
	'sub-028', 'sub-029', 'sub-030', 'sub-031', 'sub-032']

norm_stats_wf = Workflow('norm_stats_wf')
norm_stats_wf.base_dir = work_dir

subj_iterable = Node(IdentityInterface(fields = ['subject_id'], mandatory_inputs = True), 
                     name = 'subj_interable')
subj_iterable.iterables = ('subject_id', sids)


fs_skullstrip_wf = create_freesurfer_skullstrip_workflow() #skullstripping pipeline
fs_skullstrip_wf.inputs.inputspec.subjects_dir = fs_projdir
norm_stats_wf.connect(subj_iterable, 'subject_id', fs_skullstrip_wf, 'inputspec.subject_id')


info = dict(copes = [['subject_id']], #dict keys for datasource
            varcopes = [['subject_id']],
            bbreg_xfm = [['subject_id']],
            ants_warp = [['subject_id']],
            mean_image = [['subject_id']])


#acquire files
datasource = Node(DataGrabber(infields = ['subject_id'], outfields = list(info.keys())),
                       name = 'datasource')
datasource.inputs.base_directory = proj_dir                                           
datasource.inputs.field_template = dict(copes = 'GLM1/lvl2/%s/fixedfx/cope*.nii.gz', #copes and varcopes from 2nd level pipeline
                                        varcopes = 'GLM1/lvl2/%s/fixedfx/varcope*.nii.gz',                                            
                                        bbreg_xfm = 'preproc/%s/bbreg/_fs_register0/vcat*.mat', #BBReg trans matrix from preproc pipeline
                                        ants_warp = 'norm_anat/%s/anat2targ_xfm/output_Composite.h5',#ANTS trans matrix from antsreg pipeline                                          
                                        mean_image = 'preproc/%s/ref/vcat*.nii.gz') #mean reference img from preproc pipeline
datasource.inputs.sort_filelist = True
datasource.inputs.subject_id = sids
datasource.inputs.template = '*'
datasource.inputs.template_args = info
norm_stats_wf.connect(subj_iterable, 'subject_id', datasource, 'subject_id')


#convert FreeSurfer-style Affine registration into ANTS compatible itk format
convert2itk = Node(C3dAffineTool(), #ITK = Insight Tool Kit -- medical processing library
                   name = 'convert2itk')
convert2itk.inputs.fsl2ras = True #transform to ITK format 
convert2itk.inputs.itk_transform = True #export ITK transform
norm_stats_wf.connect(datasource, 'bbreg_xfm', convert2itk, 'transform_file') 
norm_stats_wf.connect(datasource, 'mean_image', convert2itk, 'source_file')
norm_stats_wf.connect(fs_skullstrip_wf, 'outputspec.skullstripped_file', convert2itk, 'reference_file')


#concatenate ITK affine and ANTS transforms into list
merge_xfm = Node(Merge(2), iterfield = ['in2'], 
                 name = 'merge_xfm')
norm_stats_wf.connect(convert2itk, 'itk_transform', merge_xfm, 'in2')
norm_stats_wf.connect(datasource, 'ants_warp', merge_xfm, 'in1')


#warp copes to target
cope2targ = MapNode(ants.ApplyTransforms(), iterfield = ['input_image'], 
                    name = 'cope2targ')
cope2targ.inputs.interpolation = 'LanczosWindowedSinc' #interpolation method used
cope2targ.inputs.invert_transform_flags = [False, False]
cope2targ.inputs.args = '--float'
cope2targ.inputs.num_threads = 4
cope2targ.inputs.dimension = 3
cope2targ.plugin_args = {'sbatch_args': '--ntasks=4'}
cope2targ.inputs.reference_image = '/home/data/madlab/Mattfeld_vCAT/derivatives/vcat_T1_grptemplate.nii.gz' #specify img whose space you are converting into
norm_stats_wf.connect(datasource, 'copes', cope2targ, 'input_image')
norm_stats_wf.connect(merge_xfm, 'out', cope2targ, 'transforms')


#warp varcopes to target
varcope2targ = MapNode(ants.ApplyTransforms(), iterfield = ['input_image'], 
                       name = 'varcope2targ')
varcope2targ.inputs.input_image_type = 3 #define input type as "timeseries"
varcope2targ.inputs.interpolation = 'LanczosWindowedSinc'
varcope2targ.inputs.invert_transform_flags = [False, False]
varcope2targ.inputs.args = '--float'
varcope2targ.inputs.num_threads = 4
varcope2targ.inputs.dimension = 3
varcope2targ.plugin_args = {'sbatch_args': '--ntasks=4 --tasks-per-node=4'}
varcope2targ.inputs.reference_image = '/home/data/madlab/Mattfeld_vCAT/derivatives/vcat_T1_grptemplate.nii.gz'
norm_stats_wf.connect(datasource, 'varcopes', varcope2targ, 'input_image')
norm_stats_wf.connect(merge_xfm, 'out', varcope2targ, 'transforms')


#define the contrasts from the names of the copes
getcontrasts = Node(Function(input_names = ['data_inputs'], output_names = ['contrasts'],
                             function = get_contrasts),
                    name = 'getcontrasts')
norm_stats_wf.connect(datasource, 'copes', getcontrasts, 'data_inputs')


#rename output files with something meaningful
getsubs = Node(Function(input_names = ['subject_id', 'cons'], output_names = ['subs'],
                        function = get_substitutions),
               name = 'getsubs')
norm_stats_wf.connect(subj_iterable, 'subject_id', getsubs, 'subject_id')
norm_stats_wf.connect(getcontrasts, 'contrasts', getsubs, 'cons')


#save output to sink_dir
sinkd = Node(DataSink(infields = None), 
                      name = 'sinkd')
sinkd.inputs.base_directory = sink_dir
sinkd.inputs.remove_dest_dir = False
norm_stats_wf.connect(subj_iterable, 'subject_id', sinkd, 'container')
norm_stats_wf.connect(cope2targ, 'output_image', sinkd, 'norm_copes')
norm_stats_wf.connect(varcope2targ, 'output_image', sinkd, 'norm_varcopes')
norm_stats_wf.connect(getsubs, 'subs', sinkd, 'substitutions')

norm_stats_wf.config['execution']['crashdump_dir'] = work_dir
norm_stats_wf.run(plugin='SLURM', plugin_args={'sbatch_args': ('-p investor --qos pq_madlab --account iacc_madlab -N 1 -n 1'), 'overwrite': True})
