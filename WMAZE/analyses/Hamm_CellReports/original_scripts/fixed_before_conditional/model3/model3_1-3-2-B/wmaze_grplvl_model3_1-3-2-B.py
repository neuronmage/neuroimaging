#!/usr/bin/env python

"""
===============================================================================
wmaze_fMRI: Mandy Thesis -- Fixed before Conditional -- Model 3 Version 1.3.2-B
===============================================================================
Group level (Random Effects) workflow for UM GE 750 wmaze task data.

- WMAZE Model 3 Version 1.3.2-B 
  - Use FSL ROI to recreate EPI data, removing last 3 volumes
  - Removed last 3 trials before EV creation
  - EV directory (Model 3) --- /home/data/madlab/data/mri/wmaze/scanner_behav/WMAZE_001/MRthesis/model3_1-3-2-B


- python wmaze_grplvl.py -s WMAZE_001
                         -o /home/data/madlab/data/mri/wmaze/grp_lvl/MRthesis_norm_stats/fixed_before_conditional/model3/model3_1-3-2-B/
                         -w /scratch/madlab/crash/wmaze_MRthesis/fixed_before_conditional/model3_1-3-2-B/grp_lvl
 
"""

from nipype.pipeline.engine import Workflow, Node, MapNode
from nipype.interfaces.utility import IdentityInterface, Merge
from nipype.interfaces.io import DataGrabber
from nipype.interfaces import fsl
from nipype.interfaces.fsl.model import L2Model, FLAMEO
from nipype.interfaces.fsl.utils import ImageMaths
from nipype.interfaces.fsl.model import SmoothEstimate, Cluster
from nipype.interfaces.io import DataSink

###############
## Functions ##
###############

# Get length
get_len = lambda x: len(x)

# Functions
def pickfirst(func):
    if isinstance(func, list):
        return func[0]
    else:
        return func


def get_substitutions(contrast):
    subs = [('_contrast{0}'.format(contrast),''),
            ('output','')]
    for i in range(0, 23):
        subs.append(('_z2pval{0}'.format(i),''))
        subs.append(('_cluster{0}'.format(i), ''))
        subs.append(('_fdr{0}'.format(i),''))
    return subs

# Array containing all 12 contrasts
contrasts = ['corr_B', 'incorr_B', 'all_remaining',  
             'all_corr_minus_incorr', 'all_incorr_minus_corr']


proj_dir = '/home/data/madlab/data/mri/wmaze'
fs_projdir = '/home/data/madlab/surfaces/wmaze'
work_dir = '/scratch/madlab/wmaze/grplvl/wmaze_MRthesis/fixed_before_conditional/model3/model3_1-3-2-B'

# Array containing subject IDs

sids = ['WMAZE_001', 'WMAZE_002', 'WMAZE_004', 'WMAZE_005', 'WMAZE_006',   
        'WMAZE_007', 'WMAZE_008', 'WMAZE_009', 'WMAZE_010', 'WMAZE_012',  
        'WMAZE_017', 'WMAZE_018', 'WMAZE_019', 'WMAZE_020', 'WMAZE_021', 
        'WMAZE_022', 'WMAZE_023', 'WMAZE_024', 'WMAZE_026', 'WMAZE_027']


# Define workflow
group_wf = Workflow('group_wf')
group_wf.base_dir = work_dir


# Node to iterate through all 13 contrasts
contrast_iterable = Node(IdentityInterface(fields = ['contrast'], 
                                           mandatory_inputs = True), 
                         name = 'contrast_iterable')
contrast_iterable.iterables = ('contrast', contrasts)

# Variable containing the dictionary keys for the datasource node
info = dict(copes = [['subject_id', 'contrast']],
            varcopes = [['subject_id', 'contrast']])

# Node to grab the cope and varcope data for each subject and contrast
datasource = Node(DataGrabber(infields = ['subject_id', 'contrast'],
                              outfields = info.keys()),
                  name = 'datasource')
# Creation of a base directory pathway
datasource.inputs.base_directory = proj_dir
# Dictionary containing the pathways necessary to obtain the files 
# %s serve as placeholders for the values in the info dictionary
datasource.inputs.field_template = dict(copes = 'norm_stats/MRthesis_norm_stats/fixed_before_conditional/model3_1-3-2-B/%s/norm_copes/cope_%s_trans.nii.gz',
                                        varcopes = 'norm_stats/MRthesis_norm_stats/fixed_before_conditional/model3_1-3-2-B/%s/norm_varcopes/varcope_%s_trans.nii.gz')
datasource.inputs.ignore_exception = False
datasource.inputs.raise_on_empty = True
# Forces DataGrabber to return data in sorted order when using wildcards
datasource.inputs.sort_filelist = True
# Indicates the string template to match (in this case, any that match the field template)
datasource.inputs.template = '*'
# Inputs from the infields argument ('subject_id', 'contrast') that satisfy the template
datasource.inputs.template_args = info
# Input for subject_id from sids
datasource.inputs.subject_id = sids
group_wf.connect(contrast_iterable, 'contrast', datasource, 'contrast')

#### contrast_iterable (contrast) ----> datasource (contrast)



# Identity interface node to hold and store information for input into future nodes
inputspec = Node(IdentityInterface(fields = ['copes', 'varcopes', 'brain_mask', 'run_mode'],
                                   mandatory_inputs = True),
                 name = 'inputspec')
# Provide the wmaze anatomical group template mask
inputspec.inputs.brain_mask = '/home/data/madlab/data/mri/wmaze/wmaze_T1_template/wmaze_grptemplate_mask.nii.gz'
# Defines "run_mode" as "flame1"
inputspec.inputs.run_mode = 'flame1'
group_wf.connect(datasource, 'copes', inputspec, 'copes')
group_wf.connect(datasource, 'varcopes', inputspec, 'varcopes')

#### datasource (copes) ----> inputspec (copes)
#### datasource (varcopes) ----> inputspec (varcopes)



# Node to concatenate the varcopes into single image across time
grp_merge_varcopes = Node(fsl.utils.Merge(), 
                          name = 'grp_merge_varcopes')
# Concatenate across time
grp_merge_varcopes.inputs.dimension = 't'
# Define dictionary for the FSL output type as NIFTI_GZ
grp_merge_varcopes.inputs.environ = {'FSLOUTPUTTYPE': 'NIFTI_GZ'}
grp_merge_varcopes.inputs.ignore_exception = False
# Declare output type as NIFTI_GZ
grp_merge_varcopes.inputs.output_type = 'NIFTI_GZ'
grp_merge_varcopes.inputs.terminal_output = 'stream'
group_wf.connect(inputspec, 'varcopes', grp_merge_varcopes, 'in_files')

#### inputspec (varcopes) ----> grp_merge_varcopes (in_files)



# Node to create group-specific level 2 model 
grp_l2model = Node(L2Model(), 
                   name = 'grp_l2model')
grp_l2model.inputs.ignore_exception = False
# get_len function provides the number of copes for the level 2 model node
group_wf.connect(inputspec, ('copes', get_len), grp_l2model, 'num_copes')

### inputspec (copes) ----> grp_l2model (num_copes)



# Node to concatenate the copes into single image across time
grp_merge_copes = Node(fsl.utils.Merge(), 
                       name = 'grp_merge_copes')
# Concatenate across time
grp_merge_copes.inputs.dimension = 't'
grp_merge_copes.inputs.environ = {'FSLOUTPUTTYPE': 'NIFTI_GZ'}
grp_merge_copes.inputs.ignore_exception = False
grp_merge_copes.inputs.output_type = 'NIFTI_GZ'
grp_merge_copes.inputs.terminal_output = 'stream'
group_wf.connect(inputspec, 'copes', grp_merge_copes, 'in_files')

#### inputspec (copes) ----> grp_merge_copes (in_files)


# Node: Randomise
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
group_wf.connect(grp_merge_copes, 'merged_file', grp_randomise, 'in_file')
group_wf.connect(inputspec, 'brain_mask', grp_randomise, 'mask')


# Node: group_randomise.sinker
group_randomise_sinker = Node(DataSink(infields = None), 
                              name = 'group_randomise_sinker')
group_randomise_sinker.inputs.base_directory = '/home/data/madlab/data/mri/wmaze/grplvl/grplvl_randomise_exhaustive/fixed_before_cond/MRthesis_z3.07/model3_1-3-2-B'
group_randomise_sinker.inputs.ignore_exception = False
group_randomise_sinker.inputs.parameterization = True
group_wf.connect(grp_randomise, 't_corrected_p_files', group_randomise_sinker, 'output.corrected.@tcorr_p_files')
group_wf.connect(grp_randomise, 'tstat_files', group_randomise_sinker, 'output.@tstat_files')
group_wf.connect(contrast_iterable, 'contrast', group_randomise_sinker, 'container')


# Run things and writes the inevitable crash files
group_wf.config['execution']['crashdump_dir'] = work_dir
# Defines the workflow directory as the working directory
group_wf.base_dir = work_dir
# Executes the script and submits to appropriate queue
group_wf.run(plugin = 'LSF', plugin_args = {'bsub_args': '-q PQ_madlab'})
