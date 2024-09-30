#!/usr/bin/env python

"""
=====================================================================================
wmaze_fMRI: Mandy Thesis -- Fixed before Conditional -- Model 3 Version 1.3.2-B_split
=====================================================================================
Group level (Random Effects) workflow for UM GE 750 wmaze task data.

- WMAZE Model 3 Version 1.3.2-B_split 
  - Use FSL ROI to recreate EPI data, removing last 3 volumes
  - Removed last 3 trials before EV creation
  - EV directory (Model 3) --- /home/data/madlab/data/mri/wmaze/scanner_behav/WMAZE_001/MRthesis/model3_1-3-2-B_split


- python wmaze_grplvl.py -s WMAZE_001
                         -o /home/data/madlab/data/mri/wmaze/grp_lvl/MRthesis_norm_stats/fixed_before_conditional/model3_1-3-2-B_split/
                         -w /scratch/madlab/crash/wmaze_MRthesis/fixed_before_conditional/model3_1-3-2-B_split/grp_lvl
 
"""

from nipype.pipeline.engine import Workflow, Node, MapNode
from nipype.interfaces.utility import IdentityInterface, Merge
from nipype.interfaces.io import DataGrabber, DataSink
from nipype.interfaces import fsl
from nipype.interfaces.utility import Function
from nipype.interfaces.fsl.model import L2Model, Randomise
from nipype.interfaces.fsl.utils import ImageMaths
 

#################
### Functions ###
#################

get_len = lambda x: len(x)


def get_substitutions(contrast):
    subs = [('_contrast_{0}'.format(contrast), ''),
            ('output','')]
    for i in range(0, 11):
        subs.append(('_z2pval{0}'.format(i), ''))
        subs.append(('_cluster{0}'.format(i), ''))
        subs.append(('_fdr{0}'.format(i), ''))
    return subs


def get_contrasts():
    cont_1st_minus_2nd = ['FrstHalf_minus_ScndHalf', 'T', ['01_lrnperiod'], [1]]
    cont_2nd_minus_1st = ['ScndHalf_minus_FrstHalf', 'T', ['01_lrnperiod'], [-1]]
    grplvl_contrasts = [cont_1st_minus_2nd, cont_2nd_minus_1st]
    return grplvl_contrasts


def get_regressors(sids):
    import numpy as np
    import collections
    reg = {}
    sids_array = np.asarray(sids)
    for sid in sids:
        curr_reg = sid == sids_array
        curr_reg_int = curr_reg.astype(int)
        curr_reg_1sthalf = curr_reg_int.tolist()
        curr_reg_2ndhalf = curr_reg_int.tolist()
        curr_reg_all = curr_reg_1sthalf + curr_reg_2ndhalf
        reg[sid] = curr_reg_all
    reg['01_lrnperiod'] = [1]*len(sids_array) + [-1]*len(sids_array)
    reg_ord = collections.OrderedDict(sorted(reg.items()))
    group = [1]*(len(sids_array)*2)
    return group, reg_ord


################
### Pipeline ###
################

contrasts = ['corr_B', 'incorr_B', 'all_remaining',  
             'all_corr_minus_incorr', 'all_incorr_minus_corr']

proj_dir = '/home/data/madlab/data/mri/wmaze'
fs_projdir = '/home/data/madlab/surfaces/wmaze'
work_dir = '/scratch/madlab/crash/wmaze_MRthesis/model3_1-3-2-B_split/grp_lvl/grp_wf'

sids = ['WMAZE_001', 'WMAZE_002', 'WMAZE_004', 'WMAZE_005', 'WMAZE_006', 
        'WMAZE_007', 'WMAZE_008', 'WMAZE_009', 'WMAZE_010', 'WMAZE_012', 
        'WMAZE_017', 'WMAZE_018', 'WMAZE_019', 'WMAZE_020', 'WMAZE_021', 
        'WMAZE_022', 'WMAZE_023', 'WMAZE_024', 'WMAZE_026', 'WMAZE_027']


# Workflow
group_wf = Workflow('group_wf')
group_wf.base_dir = work_dir

# Contrast_iterable
contrast_iterable = Node(IdentityInterface(fields = ['contrast'],
                                           mandatory_inputs = True), 
                         name = 'contrast_iterable')
contrast_iterable.iterables = ('contrast', contrasts)


info = dict(copes1sthalf = [['subject_id', 'contrast']],
            varcopes1sthalf = [['subject_id', 'contrast']],
            copes2ndhalf = [['subject_id', 'contrast']],
            varcopes2ndhalf = [['subject_id', 'contrast']])


# Datasource
datasource = Node(DataGrabber(infields = ['subject_id', 'contrast'],
                              outfields = info.keys()),
                  name = 'datasource')
datasource.inputs.base_directory = proj_dir
datasource.inputs.field_template = dict(copes1sthalf = 'norm_stats/MRthesis_norm_stats/fixed_before_conditional/model3_1-3-2-B_split/%s/norm_copes/frsthalf/cope_%s_trans.nii.gz',
                                        varcopes1sthalf = 'norm_stats/MRthesis_norm_stats/fixed_before_conditional/model3_1-3-2-B_split/%s/norm_varcopes/frsthalf/varcope_%s_trans.nii.gz',
                                        copes2ndhalf = 'norm_stats/MRthesis_norm_stats/fixed_before_conditional/model3_1-3-2-B_split/%s/norm_copes/scndhalf/cope_%s_trans.nii.gz',
                                        varcopes2ndhalf = 'norm_stats/MRthesis_norm_stats/fixed_before_conditional/model3_1-3-2-B_split/%s/norm_varcopes/scndhalf/varcope_%s_trans.nii.gz')
datasource.inputs.ignore_exception = False
datasource.inputs.raise_on_empty = True
datasource.inputs.sort_filelist = True
datasource.inputs.template = '*'
datasource.inputs.template_args = info
datasource.inputs.subject_id = sids
group_wf.connect(contrast_iterable, 'contrast', datasource, 'contrast')


# Function node to define the contrasts for the experiment
getcontrasts = Node(Function(input_names = [],
                             output_names = ['grouplvl_contrasts'],
                             function = get_contrasts),
                    name = 'getcontrasts')
getcontrasts.inputs.ignore_exception = False


# Function node to define the regressors for the experiment
getregressors = Node(Function(input_names = ['sids'],
                              output_names = ['group', 'regressors'],
                              function = get_regressors),
                    name = 'getregressors')
getregressors.inputs.ignore_exception = False
getregressors.inputs.sids = sids

#group, regressors = get_regressors(sids)


# Identity interface node to hold and store information for input into future nodes
inputspec = Node(IdentityInterface(fields = ['copes1st', 'varcopes1st', 'copes2nd', 
                                             'varcopes2nd', 'brain_mask', 'run_mode',
                                             'group', 'regressors', 'group_contrasts'],
                                   mandatory_inputs = True),
                 name = 'inputspec')
# Provide the wmaze anatomical group template mask
inputspec.inputs.brain_mask = '/home/data/madlab/data/mri/wmaze/wmaze_T1_template/wmaze_grptemplate_mask.nii.gz'
inputspec.inputs.run_mode = 'flame1'
group_wf.connect(datasource, 'copes1sthalf', inputspec, 'copes1st')
group_wf.connect(datasource, 'varcopes1sthalf', inputspec, 'varcopes1st')
group_wf.connect(datasource, 'copes2ndhalf', inputspec, 'copes2nd')
group_wf.connect(datasource, 'varcopes2ndhalf', inputspec, 'varcopes2nd')
group_wf.connect(getcontrasts, 'grouplvl_contrasts', inputspec, 'group_contrasts')
group_wf.connect(getregressors, 'regressors', inputspec, 'regressors') 


# Create a merge node to create a single list of 1st and 2nd half copes for the fsl merge node below
merge_1st2ndhalf_copes = Node(Merge(2), 
                              name = 'merge_1st2ndhalf_copes')
group_wf.connect(inputspec, 'copes1st', merge_1st2ndhalf_copes, 'in1')
group_wf.connect(inputspec, 'copes2nd', merge_1st2ndhalf_copes, 'in2')


# Merge 1st and 2nd half copes into a single list for the fsl merge node below
merge_1st2ndhalf_varcopes = Node(Merge(2), 
                                 name = 'merge_1st2ndhalf_varcopes')
group_wf.connect(inputspec, 'varcopes1st', merge_1st2ndhalf_varcopes, 'in1')
group_wf.connect(inputspec, 'varcopes2nd', merge_1st2ndhalf_varcopes, 'in2')


# Concatenate varcopes into single image across time
merge_ALL_varcopes = Node(fsl.utils.Merge(), 
                          name = 'merge_1sthalf_varcopes')
merge_ALL_varcopes.inputs.dimension = 't'
merge_ALL_varcopes.inputs.environ = {'FSLOUTPUTTYPE': 'NIFTI_GZ'}
merge_ALL_varcopes.inputs.ignore_exception = False
merge_ALL_varcopes.inputs.output_type = 'NIFTI_GZ'
merge_ALL_varcopes.inputs.terminal_output = 'stream'
group_wf.connect(inputspec, 'varcopes1st', merge_ALL_varcopes, 'in_files')


# Concatenate copes into single image across time
merge_ALL_copes = Node(fsl.utils.Merge(), 
                       name = 'grp_merge_copesALL')
merge_ALL_copes.inputs.dimension = 't'
merge_ALL_copes.inputs.environ = {'FSLOUTPUTTYPE': 'NIFTI_GZ'}
merge_ALL_copes.inputs.ignore_exception = False
merge_ALL_copes.inputs.output_type = 'NIFTI_GZ'
merge_ALL_copes.inputs.terminal_output = 'stream'
group_wf.connect(merge_1st2ndhalf_copes, 'out', merge_ALL_copes, 'in_files')


# Group l2model
grp_l2model = Node(fsl.MultipleRegressDesign(), 
                   name = 'grp_l2model')
grp_l2model.inputs.ignore_exception = False
group_wf.connect(inputspec, 'group_contrasts', grp_l2model, 'contrasts')
group_wf.connect(inputspec, 'regressors', grp_l2model, 'regressors')
#group_wf.connect(inputspec, 'group', grp_l2model, 'groups')


# Randomise
grp_randomise = Node(fsl.Randomise(), 
                     name = 'grp_randomise')
grp_randomise.inputs.environ = {'FSLOUTPUTTYPE': 'NIFTI_GZ'}
grp_randomise.inputs.ignore_exception = False
grp_randomise.inputs.tfce = True
grp_randomise.inputs.base_name = 'oneSampT'
grp_randomise.inputs.num_perm = 5000
grp_randomise.inputs.output_type = 'NIFTI_GZ'
grp_randomise.inputs.terminal_output = 'stream'
grp_randomise.plugin_args = {'bsub_args': '-m IB_40C_1.5T_1'}
group_wf.connect(grp_l2model, 'design_mat', grp_randomise, 'design_mat')
group_wf.connect(grp_l2model, 'design_con', grp_randomise, 'tcon')
group_wf.connect(merge_ALL_copes, 'merged_file', grp_randomise, 'in_file')
group_wf.connect(inputspec, 'brain_mask', grp_randomise, 'mask')


# Randomise sinker
group_randomise_sinker = Node(DataSink(infields = None), 
                              name = 'group_randomise_sinker')
group_randomise_sinker.inputs.base_directory = '/home/data/madlab/data/mri/wmaze/grplvl/grplvl_randomise_exhaustive/MRthesis_z3.07/model3_1-3-2-B_split/all'
group_randomise_sinker.inputs.ignore_exception = False
group_randomise_sinker.inputs.parameterization = True
group_wf.connect(grp_randomise, 't_corrected_p_files', group_randomise_sinker, 'output.corrected.@tcorr_p_files')
group_wf.connect(grp_randomise, 'tstat_files', group_randomise_sinker, 'output.@tstat_files')
group_wf.connect(contrast_iterable, 'contrast', group_randomise_sinker, 'container')


group_wf.config['execution']['crashdump_dir'] = work_dir
group_wf.base_dir = work_dir
group_wf.run(plugin = 'LSF', plugin_args = {'bsub_args': '-q PQ_madlab'})
