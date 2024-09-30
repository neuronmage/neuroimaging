#!/usr/bin/env python

"""
=============================================
GLM2_fMRI: Fixed vs. Conditional -- vCAT GLM2
=============================================
Group level (Random Effects) workflow for Seimens vCAT task data

- vCAT GLM2 
  - Use FSL ROI to recreate EPI data
  - EV directory (GLM2) --- /home/data/madlab/Mattfeld_vCAT/behav/sub-0??/GLM2

- python GLM2_grplvl.py -s sub-0??
                        -o /home/data/madlab/Mattfeld_vCAT/derivatives/GLM2/grp_lvl
                        -w /scratch/madlab/crash/mandy/vcat/GLM2/grp_lvl
 
"""

from nipype.pipeline.engine import Workflow, Node, MapNode
from nipype.interfaces.utility import IdentityInterface, Merge
from nipype.interfaces.io import DataGrabber, DataSink
from nipype.interfaces import fsl
from nipype.interfaces.fsl.model import L2Model, FLAMEO, Randomise
from nipype.interfaces.fsl.utils import ImageMaths


###############
## Functions ##
###############


get_len = lambda x: len(x) #get length


def pickfirst(func): #pick first volume in a timeseries
    if isinstance(func, list):
        return func[0]
    else:
        return func


def get_substitutions(contrast):
    subs = [('_contrast{0}'.format(contrast),''), ('output','')]
    for i in range(0, 23):
        subs.append(('_z2pval{0}'.format(i),''))
        subs.append(('_cluster{0}'.format(i), ''))
        subs.append(('_fdr{0}'.format(i),''))
    return subs


contrasts = ['AllVsBase', 'fixed_corr', 'fixed_incorr', 'cond_corr', 'cond_incorr', 'remaining',
             'fixedCorr_minus_fixedIncorr', 'fixedIncorr_minus_fixedCorr', 'all_fixed',
             'condCorr_minus_condIncorr', 'condIncorr_minus_condCorr', 'all_cond',
             'fixedCorr_minus_condCorr', 'condCorr_minus_fixedCorr', 'all_corr',
             'fixedIncorr_minus_condIncorr', 'condIncorr_minus_fixedIncorr', 'all_incorr',
             'allCorr_minus_allIncorr', 'allIncorr_minus_allCorr',
             'allFixed_minus_allCond', 'allCond_minus_allFixed']


proj_dir = '/home/data/madlab/Mattfeld_vCAT/derivatives'
fs_projdir = '/home/data/madlab/Mattfeld_vCAT/derivatives/freesurfer'
work_dir = '/scratch/madlab/crash/mandy/vcat/GLM2/grp_lvl'
sink_dir = '/home/data/madlab/Mattfeld_vCAT/derivatives/GLM2/grp_lvl'

sids = ['sub-005', 'sub-006', 'sub-007', 'sub-008', 'sub-010', 
        'sub-012', 'sub-013', 'sub-014', 'sub-015', 'sub-016', 
        'sub-018', 'sub-019', 'sub-020', 'sub-021', 'sub-022', 
        'sub-023', 'sub-024', 'sub-025', 'sub-026', 'sub-027', 
        'sub-028', 'sub-029', 'sub-030', 'sub-031', 'sub-032']


group_wf = Workflow('group_wf')
group_wf.base_dir = work_dir


#iterate through contrasts
contrast_iterable = Node(IdentityInterface(fields = ['contrast'], mandatory_inputs = True), 
                         name = 'contrast_iterable')
contrast_iterable.iterables = ('contrast', contrasts)


#dictionary keys for datasource
info = dict(copes = [['subject_id', 'contrast']], 
            varcopes = [['subject_id', 'contrast']])

#grab cope and varcope data for each subject and contrast
datasource = Node(DataGrabber(infields = ['subject_id', 'contrast'], outfields = list(info.keys())),
                  name = 'datasource')
datasource.inputs.base_directory = proj_dir
datasource.inputs.field_template = dict(copes = 'GLM2/norm_stats/%s/norm_copes/cope_%s_trans.nii.gz',
                                        varcopes = 'GLM2/norm_stats/%s/norm_varcopes/varcope_%s_trans.nii.gz')
datasource.inputs.raise_on_empty = True
datasource.inputs.sort_filelist = True
datasource.inputs.template = '*'
datasource.inputs.template_args = info
datasource.inputs.subject_id = sids
group_wf.connect(contrast_iterable, 'contrast', datasource, 'contrast')


#hold and store info for input into future nodes
inputspec = Node(IdentityInterface(fields = ['copes', 'varcopes', 'brain_mask', 'run_mode'], mandatory_inputs = True),
                 name = 'inputspec')
inputspec.inputs.brain_mask = '/home/data/madlab/Mattfeld_vCAT/derivatives/vcat_T1_grptemplate.nii.gz' #anatomical group template mask
inputspec.inputs.run_mode = 'flame1'
group_wf.connect(datasource, 'copes', inputspec, 'copes')
group_wf.connect(datasource, 'varcopes', inputspec, 'varcopes')


#concatenate varcopes into single image across time
merge_varcopes = Node(fsl.utils.Merge(), 
                      name = 'merge_varcopes')
merge_varcopes.inputs.dimension = 't' #concatenate across time
merge_varcopes.inputs.environ = {'FSLOUTPUTTYPE': 'NIFTI_GZ'}
merge_varcopes.inputs.output_type = 'NIFTI_GZ'
group_wf.connect(inputspec, 'varcopes', merge_varcopes, 'in_files')


#group-specific level 2 model 
l2model = Node(L2Model(), 
               name = 'l2model')
group_wf.connect(inputspec, ('copes', get_len), l2model, 'num_copes')


#concatenate copes into single image across time
merge_copes = Node(fsl.utils.Merge(), 
                   name = 'merge_copes')
merge_copes.inputs.dimension = 't' #concatenate across time
merge_copes.inputs.environ = {'FSLOUTPUTTYPE': 'NIFTI_GZ'}
merge_copes.inputs.output_type = 'NIFTI_GZ'
group_wf.connect(inputspec, 'copes', merge_copes, 'in_files')


#execute Randomise
randomise = Node(fsl.Randomise(), 
                 name = 'randomise')
randomise.inputs.environ = {'FSLOUTPUTTYPE': 'NIFTI_GZ'}
randomise.inputs.tfce = True
randomise.inputs.base_name = 'oneSampT'
randomise.inputs.num_perm = 5000
randomise.inputs.output_type = 'NIFTI_GZ'
randomise.plugin_args = {'sbatch_args': '-p IB_40C_1.5T'}
group_wf.connect(l2model, 'design_mat', randomise, 'design_mat')
group_wf.connect(l2model, 'design_con', randomise, 'tcon')
group_wf.connect(merge_copes, 'merged_file', randomise, 'in_file')
group_wf.connect(inputspec, 'brain_mask', randomise, 'mask')


#save output to sink_dir
sinkd = Node(DataSink(infields = None), 
             name = 'sinkd')
sinkd.inputs.base_directory = sink_dir
group_wf.connect(randomise, 't_corrected_p_files', sinkd, 'output.corrected.@tcorr_p_files')
group_wf.connect(randomise, 'tstat_files', sinkd, 'output.@tstat_files')
group_wf.connect(contrast_iterable, 'contrast', sinkd, 'container')

group_wf.config['execution']['crashdump_dir'] = '/scratch/madlab/crash/mandy/vcat/GLM2/grp_lvl'
group_wf.base_dir = work_dir
group_wf.run(plugin='SLURM', plugin_args={'sbatch_args': ('-p investor --qos pq_madlab --account iacc_madlab -N 1 -n 1'), 'overwrite': True})
