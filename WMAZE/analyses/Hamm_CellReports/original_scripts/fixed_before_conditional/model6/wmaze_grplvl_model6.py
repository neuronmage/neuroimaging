#!/usr/bin/env python

"""
=============================================================================
wmaze_fMRI: Mandy Thesis -- Fixed before Conditional -- Model 6 Version 1.3.2
=============================================================================
Group level (Random Effects) workflow for UM GE 750 wmaze task data.

- WMAZE Model 6 Version 1.3.2 
  - Use FSL ROI to recreate EPI data, removing last 3 volumes
  - Removed last 3 trials before EV creation
  - EV directory (Model 6) --- /home/data/madlab/data/mri/wmaze/scanner_behav/WMAZE_001/MRthesis/model6


- python wmaze_grplvl.py -s WMAZE_001
                         -o /home/data/madlab/data/mri/wmaze/grp_lvl/MRthesis_norm_stats/fixed_before_conditional/model6
                         -w /scratch/madlab/crash/wmaze_MRthesis/fixed_before_conditional/model6/grp_lvl
 
"""

from nipype.pipeline.engine import Workflow, Node, MapNode
from nipype.interfaces.utility import IdentityInterface, Merge
from nipype.interfaces.io import DataGrabber
from nipype.interfaces import fsl
from nipype.interfaces.fsl.model import L2Model
from nipype.interfaces.fsl.model import FLAMEO
from nipype.interfaces.fsl.utils import ImageMaths
from nipype.interfaces.fsl.model import SmoothEstimate
from nipype.interfaces.fsl.model import Cluster
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

# The contrast names 
contrasts = ['AllVsBase', 'fixed_corr_before_cond_corr', 'fixed_corr_before_cond_incorr', 'all_remaining',
                 'fixed_b4_condCorr_minus_condIncorr', 'fixed_b4_condIncorr_minus_condCorr']


proj_dir = '/home/data/madlab/data/mri/wmaze'
fs_projdir = '/home/data/madlab/surfaces/wmaze'
work_dir = '/scratch/madlab/wmaze/grplvl/wmaze_MRthesis/fixed_before_conditional/model6'

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
datasource.inputs.field_template = dict(copes = 'norm_stats/MRthesis_norm_stats/fixed_before_conditional/model6/%s/norm_copes/cope_%s_trans.nii.gz',
                                        varcopes = 'norm_stats/MRthesis_norm_stats/fixed_before_conditional/model6/%s/norm_varcopes/varcope_%s_trans.nii.gz')
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


'''
# Node: flameo
grp_flameo = Node(FLAMEO(), 
                  name = 'grp_flameo')
grp_flameo.inputs.environ = {'FSLOUTPUTTYPE': 'NIFTI_GZ'}
grp_flameo.inputs.ignore_exception = False
grp_flameo.inputs.log_dir = 'stats'
grp_flameo.inputs.output_type = 'NIFTI_GZ'
grp_flameo.inputs.terminal_output = 'stream'
# Design matrix file
group_wf.connect(grp_l2model, 'design_mat', grp_flameo, 'design_file')
# ASCII matrix specifying t-contrasts
group_wf.connect(grp_l2model, 'design_con', grp_flameo, 't_con_file')
# ASCII matrix specifying the groups the covariance is split into
group_wf.connect(grp_l2model, 'design_grp', grp_flameo, 'cov_split_file')
# Cope regressor data file
group_wf.connect(grp_merge_copes, 'merged_file', grp_flameo, 'cope_file')
# Varcope weightings data file
group_wf.connect(grp_merge_varcopes, 'merged_file', grp_flameo, 'var_cope_file')
# Inference to perform
group_wf.connect(inputspec, 'run_mode', grp_flameo, 'run_mode')
# Brain mask file
group_wf.connect(inputspec, 'brain_mask', grp_flameo, 'mask_file')

#### grp_l2model (design_mat) ----> grp_flameo (design_file)
#### grp_l2model (design_con) ----> grp_flameo (t_con_file)
#### grp_l2model (design_grp) ----> grp_flameo (cov_split_file)
#### grp_merge_copes (merged_file) ----> grp_flameo (cope_file)
#### grp_merge_varcopes (merged_file) ----> grp_flameo (var_cope_file)
#### inputspec (run_mode) ----> grp_flameo (run_mode)
#### inputspec (brain_mask) ----> grp_flameo (mask_file)



# Node to convert group z-statistic from Flameo into uncorrected p-value 
grp_z2pval = MapNode(ImageMaths(), 
                     iterfield = ['in_file'], 
                     name = 'grp_z2pval')
grp_z2pval.inputs.environ = {'FSLOUTPUTTYPE': 'NIFTI_GZ'}
grp_z2pval.inputs.ignore_exception = False
# String defining the operation
grp_z2pval.inputs.op_string = '-ztop'
grp_z2pval.inputs.output_type = 'NIFTI_GZ'
# Suffix for the output file
grp_z2pval.inputs.suffix = '_pval'
grp_z2pval.inputs.terminal_output = 'stream'
# Z-statistic from FLAMEO
group_wf.connect(grp_flameo, 'zstats', grp_z2pval, 'in_file')

#### grp_flameo (zstats) ----> grp_z2pval (in_file)



# IdentityInterface node to collect and hold output from the z2pval & FLAMEO nodes 
outputspec = Node(IdentityInterface(fields = ['zstat', 'tstat', 'cope', 'varcope', 'mrefvars', 
                                              'pes', 'res4d', 'mask', 'tdof', 'weights', 'pstat'],
                                    mandatory_inputs = True),
                  name = 'outputspec')
# P-value converted from group z-stat
group_wf.connect(grp_z2pval, 'out_file', outputspec, 'pstat')
# Contrast estimates for each contrast
group_wf.connect(grp_flameo, 'copes', outputspec, 'cope')
# Variance estimates for each contrast
group_wf.connect(grp_flameo, 'var_copes', outputspec, 'varcope')
# f-stat file for each contrast
group_wf.connect(grp_flameo, 'mrefvars', outputspec, 'mrefvars')
# Parameter estimates for each column of the design matrix for each voxel
group_wf.connect(grp_flameo, 'pes', outputspec, 'pes')
# Model fit residual mean-squared error for each time point
group_wf.connect(grp_flameo, 'res4d', outputspec, 'res4d')
# Weights file for each contrast
group_wf.connect(grp_flameo, 'weights', outputspec, 'weights')
# z-stat file for each contrast
group_wf.connect(grp_flameo, 'zstats', outputspec, 'zstat')
# t-stat file for each contrast
group_wf.connect(grp_flameo, 'tstats', outputspec, 'tstat')
# Temporal degrees of freedom (dof) file for each contrast
group_wf.connect(grp_flameo, 'tdof', outputspec, 'tdof')

#### grp_z2pval (out_file) ----> outputspec (pstat)
#### grp_flameo (copes) ----> outputspec (cope)
#### grp_flameo (var_copes) ----> outputspec (varcope)
#### grp_flameo (mrefvars) ----> outputspec (mrefvars)
#### grp_flameo (pes) ----> outputspec (pes)
#### grp_flameo (res4d) ----> outputspec (res4d)
#### grp_flameo (weights) ----> outputspec (weights)
#### grp_flameo (zstats) ----> outputspec (zstat)
#### grp_flameo (tstats) ----> outputspec (tstat)
#### grp_flameo (tdof) ----> outputspec (tdof)



# IdentityInterface node for holding input for the smooth_estimate and cluster nodes
cluster_inputspec = Node(IdentityInterface(fields = ['zstat', 'mask', 'zthreshold', 
                                                     'pthreshold', 'connectivity'],
                                           mandatory_inputs = True),
                         name = 'cluster_inputspec')
# Define connectivity as 26
cluster_inputspec.inputs.connectivity = 26
# Group anatomical template mask
cluster_inputspec.inputs.mask = '/home/data/madlab/data/mri/wmaze/wmaze_T1_template/wmaze_grptemplate_mask.nii.gz'
# Define p-value threshold as 0.05
cluster_inputspec.inputs.pthreshold = 0.05
# Define z-stat threshold as 2.3 or 3.07 (typical thresholds)
cluster_inputspec.inputs.zthreshold = 3.07
group_wf.connect(outputspec, 'zstat', cluster_inputspec, 'zstat')

#### outputspec (zstat) ----> cluster_inputspec (zstat)



# Mapnode to estimate the smoothness of each mask/zstat file
smooth_estimate = MapNode(SmoothEstimate(), iterfield = ['zstat_file'],
                          name = 'smooth_estimate')
smooth_estimate.inputs.environ = {'FSLOUTPUTTYPE': 'NIFTI_GZ'}
smooth_estimate.inputs.ignore_exception = False
smooth_estimate.inputs.output_type = 'NIFTI_GZ'
smooth_estimate.inputs.terminal_output = 'stream'
group_wf.connect(cluster_inputspec, 'zstat', smooth_estimate, 'zstat_file')
group_wf.connect(cluster_inputspec, 'mask', smooth_estimate, 'mask_file')

#### cluster_inputspec (zstat) ----> smooth_estimate (zstat_file)
#### cluster_inputspec (mask) ----> smooth_estimate (mask_file)



# Mapnode to perform clustering on statistical output
cluster = MapNode(Cluster(), iterfield = ['in_file', 'dlh', 'volume'], 
                  name = 'cluster')
cluster.inputs.environ = {'FSLOUTPUTTYPE': 'NIFTI_GZ'}
cluster.inputs.ignore_exception = False
# Output of cluster index (in size order -- Boolean)
cluster.inputs.out_index_file = True
# Local maxima text file (Boolean)
cluster.inputs.out_localmax_txt_file = True
# Output of local maxima volume (Boolean)
cluster.inputs.out_localmax_vol_file = True
# Output of local maxima volume (Boolean)
cluster.inputs.out_pval_file = True
# Thresholded image (Boolean)
cluster.inputs.out_threshold_file = True
# FSL output type
cluster.inputs.output_type = 'NIFTI_GZ'
cluster.inputs.terminal_output = 'stream'
# Threshold for input volume
group_wf.connect(cluster_inputspec, 'zthreshold', cluster, 'threshold')
# p-threshold for clusters
group_wf.connect(cluster_inputspec, 'pthreshold', cluster, 'pthreshold')
# The connectivity of voxels (default 26)
group_wf.connect(cluster_inputspec, 'connectivity', cluster, 'connectivity')
# Input volume
group_wf.connect(cluster_inputspec, 'zstat', cluster, 'in_file')
# Smoothness estimate = sqrt(det(Lambda))
group_wf.connect(smooth_estimate, 'dlh', cluster, 'dlh')
# Number of voxels in the mask
group_wf.connect(smooth_estimate, 'volume', cluster, 'volume')

#### cluster_inputspec (zthreshold) ----> cluster (threshold)
#### cluster_inputspec (pthreshold) ----> cluster (pthreshold)
#### cluster_inputspec (connectivity) ----> cluster (connectivity)
#### cluster_inputspec (zstat) ----> cluster (in_file)
#### smooth_estimate (dlh) ----> cluster (dlh)
#### smooth_estimate (volume) ----> cluster (volume)



# Stores the output of cluster node
cluster_outputspec = Node(IdentityInterface(fields = ['corrected_z', 'localmax_txt', 'index_file', 'localmax_vol'],
                                            mandatory_inputs = True),
                          name = 'cluster_outputspec')
# Thresholded image as the corrected z-stat file
group_wf.connect(cluster, 'threshold_file', cluster_outputspec, 'corrected_z')
# Output of cluster index (in size order)
group_wf.connect(cluster, 'index_file', cluster_outputspec, 'index_file')
# Output of local maxima volume
group_wf.connect(cluster, 'localmax_vol_file', cluster_outputspec, 'localmax_vol')
# Output of local maxima text file
group_wf.connect(cluster, 'localmax_txt_file', cluster_outputspec, 'localmax_txt')

#### cluster (threshold_file) ----> cluster_outputspec (corrected_z)
#### cluster (index_file) ----> cluster_outputspec (index_file)
#### cluster (localmax_vol_file) ----> cluster_outputspec (localmax_vol)
#### cluster (localmax_txt_file) ----> cluster_outputspec (localmax_txt)
'''

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
group_randomise_sinker.inputs.base_directory = '/home/data/madlab/data/mri/wmaze/grplvl/grplvl_MRthesis_randomise/model6'
group_randomise_sinker.inputs.ignore_exception = False
group_randomise_sinker.inputs.parameterization = True
group_wf.connect(grp_randomise, 't_corrected_p_files', group_randomise_sinker, 'output.corrected.@tcorr_p_files')
group_wf.connect(grp_randomise, 'tstat_files', group_randomise_sinker, 'output.@tstat_files')
group_wf.connect(contrast_iterable, 'contrast', group_randomise_sinker, 'container')

'''
# Save the relevant data into an output directory
group_sinker = Node(DataSink(infields = None), 
                    name = 'group_sinker')
# A dictionary with keys which are a unicode string and with values which are any value, nipype default value: {}
group_sinker.inputs._outputs = {}
# Sets the base directory for all datasink output folders
group_sinker.inputs.base_directory = '/home/data/madlab/data/mri/wmaze/grplvl/grplvl_MRthesis_randomise/model6'
group_sinker.inputs.ignore_exception = False
# Stores output in parametrized structure
group_sinker.inputs.parameterization = True
# Remove dest directory when copying dirs
group_sinker.inputs.remove_dest_dir = False
group_wf.connect(cluster_outputspec, 'corrected_z', group_sinker, 'output.corrected.@zthresh')
group_wf.connect(cluster_outputspec, 'localmax_txt', group_sinker, 'output.corrected.@localmax_txt')
group_wf.connect(cluster_outputspec, 'index_file', group_sinker, 'output.corrected.@index')
group_wf.connect(cluster_outputspec, 'localmax_vol', group_sinker, 'output.corrected.@localmax_vol')
group_wf.connect(contrast_iterable, ('contrast', get_substitutions), group_sinker, "substitutions")
group_wf.connect(contrast_iterable, 'contrast', group_sinker, 'container')
group_wf.connect(outputspec, 'cope', group_sinker, 'output.@cope')
group_wf.connect(outputspec, 'varcope', group_sinker, 'output.@varcope')
group_wf.connect(outputspec, 'mrefvars', group_sinker, 'output.@mrefvars')
group_wf.connect(outputspec, 'pes', group_sinker, 'output.@pes')
group_wf.connect(outputspec, 'res4d', group_sinker, 'output.@res4d')
group_wf.connect(outputspec, 'weights', group_sinker, 'output.@weights')
group_wf.connect(outputspec, 'zstat', group_sinker, 'output.@zstat')
group_wf.connect(outputspec, 'tstat', group_sinker, 'output.@tstat')
group_wf.connect(outputspec, 'pstat', group_sinker, 'output.@pstat')
group_wf.connect(outputspec, 'tdof', group_sinker, 'output.@tdof')
'''


#### cluster_outputspec (corrected_z) ----> group_sinker (output.corrected.@zthresh)
#### cluster_outputspec (localmax_txt) ----> group_sinker (output.corrected.@localmax_txt)
#### cluster_outputspec (index_file) ----> group_sinker (output.corrected.@index)
#### cluster_outputspec (localmax_vol) ----> group_sinker (output.corrected.@localmax_vol)
#### contrast_iterable (contrast, get_substitution) ----> group_sinker (substitutions)
#### contrast_iterable (contrast) ----> group_sinker (container)
#### outputspec (cope) ----> group_sinker (output.@cope)
#### outputspec (varcope) ----> group_sinker (output.@varcope)
#### outputspec (mrefvars) ----> group_sinker (output.@mrefvars)
#### outputspec (pes) ----> group_sinker (output.@pes)
#### outputspec (res4d) ----> group_sinker (output.@res4d)
#### outputspec (weights) ----> group_sinker (output.@weights)
#### outputspec (zstat) ----> group_sinker (output.@zstat)
#### outputspec (pstat) ----> group_sinker (output.@pstat)
#### outputspec (tdof) ----> group_sinker (output.@tdof)


# Run things and writes the inevitable crash files
group_wf.config['execution']['crashdump_dir'] = work_dir
# Defines the workflow directory as the working directory
group_wf.base_dir = work_dir
# Executes the script and submits to appropriate queue
group_wf.run(plugin = 'LSF', plugin_args = {'bsub_args': '-q PQ_madlab'})


