#!/usr/bin/env python

#SBATCH --partition IB_44C_512G
#SBATCH --account iacc_madlab
#SBATCH --qos pq_madlab
#SBATCH -o /scratch/madlab/crash/mandy/vcat/norm_anat/antsreg/antsnorm_out
#SBATCH -e /scratch/madlab/crash/mandy/vcat/norm_anat/antsreg/antsnorm_err

'''
ANTS Registration registers a moving image to a fixed image, using a predefined cost function(s) and transformation operations.
The cost function is defined using one or more 'metrics', specifically local cross-correlation (CC), mean square (MeansSquare), Demons (Demons),
global correlation (GC), or mutual information (Mattes or MI).
ANTS can use both linear (Translation, Rigid, Affine, CompositeAffine, or Translation) and nonlinear (BSpline, GaussianDisplacementField, 
TimeVaryingBSplineVelocityField, SyN, BSplineSyN, Exponential, or BSplineExponential) transformations. Usually, registration is performed in multiple stages.
'''

import os
from glob import glob
from nipype.interfaces import fsl
from nipype.interfaces import ants
from nipype.pipeline.engine import Workflow, Node, MapNode, JoinNode
from nipype.interfaces.utility import IdentityInterface, Function
from nipype.interfaces.io import DataSink
from mattfeld_utility_workflows.fs_skullstrip_util import create_freesurfer_skullstrip_workflow


fs_projdir = '/home/data/madlab/Mattfeld_vCAT/derivatives/freesurfer'
projdir = '/home/data/madlab/Mattfeld_vCAT'
crashdir = '/scratch/madlab/crash/mandy/vcat/norm_anat'
sinkdir = '/home/data/madlab/Mattfeld_vCAT/derivatives'
workdir = os.path.join(crashdir,'antsreg') #workflow output directory
if not os.path.exists(workdir):
    os.makedirs(workdir)


fs_skullstrip_wf = create_freesurfer_skullstrip_workflow()
fs_skullstrip_wf.inputs.inputspec.subjects_dir = fs_projdir


sids = ['vCAT_005', 'vCAT_006', 'vCAT_007', 'vCAT_008', 'vCAT_010', 'vCAT_012', 'vCAT_013', 'vCAT_014', 'vCAT_015', 'vCAT_016',  
	'vCAT_018', 'vCAT_019', 'vCAT_020', 'vCAT_021', 'vCAT_022', 'vCAT_023', 'vCAT_024', 'vCAT_025', 'vCAT_026', 'vCAT_027',  
	'vCAT_028', 'vCAT_029', 'vCAT_030', 'vCAT_031', 'vCAT_032']


#set up overall workflow
antsreg_wf = Workflow(name = 'antsreg_wf')
antsreg_wf.base_dir = workdir


subjID_infosource = Node(IdentityInterface(fields = ['subject_id','subjects_dir']), 
                         name = 'subjID_infosource')
subjID_infosource.iterables = ('subject_id', sids)
antsreg_wf.connect(subjID_infosource, 'subject_id', fs_skullstrip_wf, 'inputspec.subject_id')


#aggregrate all outputs from fs_skullstrip_wf
reg = Node(ants.Registration(), 
           name = 'antsRegister')
reg.inputs.fixed_image = os.path.join(projdir, 'derivatives/T1_grptemplate/vcat_T1_grptemplate.nii.gz')
reg.inputs.metric = ['Mattes'] * 2 + [['Mattes', 'CC']] #metric(s) for each stage
reg.inputs.metric_weight = [1] * 2 + [[0.5, 0.5]] #metric weights for each stage
reg.inputs.shrink_factors = [[3, 2, 1]]*2 + [[4, 2, 1]]
reg.inputs.smoothing_sigmas = [[4, 2, 1]] * 2 + [[1, 0.5, 0]]
reg.inputs.transforms = ['Rigid', 'Affine', 'SyN'] #list of transforms to be performed, in order
reg.inputs.collapse_output_transforms = True #combines all adjacent linear transforms and composes all adjacent displacement field transforms before writing results
reg.inputs.convergence_threshold = [1.e-8] * 2 + [-0.01] #requires number_of_iterations
reg.inputs.convergence_window_size = [20] * 2 + [5] #requires convergence_threshold
reg.inputs.dimension = 3 #image dimension [2 or 3]
reg.inputs.float = True #use float instead of double for computations
reg.inputs.initial_moving_transform_com = True #align moving_image and fixed_image before registration using geometric center of images [0], image intensities [1], or image origin [2]
#reg.inputs.num_threads = 4 #number of ITK threads to use
reg.inputs.number_of_iterations = [[10000, 11110, 11110]] * 2 + [[100, 100, 50]]
reg.inputs.output_transform_prefix = 'output_'
reg.inputs.radius_or_number_of_bins = [32] * 2 + [[32, 4]] #number of bins in each stage for MI and Mattes metric, radius for other metrics
reg.inputs.sampling_strategy = ['Regular'] * 2 + [[None, None]] #metric sampling strategies for each stage
reg.inputs.sampling_percentage = [0.3] * 2 + [[None, None]] #metric sampling percentage(s) for each stage
reg.inputs.sigma_units = ['vox'] * 3 #units for smoothing sigmas
reg.inputs.transform_parameters = [(0.1,), (0.1,), (0.2, 3.0, 0.0)] #list of tuples
reg.inputs.use_estimate_learning_rate_once = [True] * 3
reg.inputs.use_histogram_matching = [False] * 2 + [True] #histogram match images before registration
reg.inputs.winsorize_lower_quantile = 0.005 #lower quartile to clip image ranges
reg.inputs.winsorize_upper_quantile = 0.995 #upper quartile to clip image ranges
reg.inputs.write_composite_transform = True
reg.inputs.output_warped_image = 'output_warped_image.nii.gz'
antsreg_wf.connect(fs_skullstrip_wf, 'outputspec.skullstripped_file', reg, 'moving_image')


#move results to designated folder
datasink = Node(DataSink(), 
                name = 'datasink')
datasink.inputs.base_directory = os.path.join(sinkdir, 'proj_tmpl')
antsreg_wf.connect(subjID_infosource, 'subject_id', datasink, 'container')
antsreg_wf.connect(reg, 'composite_transform', datasink, 'anat2targ_xfm')
antsreg_wf.connect(reg, 'inverse_composite_transform', datasink, 'targ2anat_xfm')
antsreg_wf.connect(reg, 'warped_image', datasink, 'warped_image')


#execute workflow
antsreg_wf.run(plugin='SLURM', plugin_args={'sbatch_args': ('--partition IB_44C_512G --account iacc_madlab --qos pq_madlab -N 1 -n 1'), 'overwrite': True})

