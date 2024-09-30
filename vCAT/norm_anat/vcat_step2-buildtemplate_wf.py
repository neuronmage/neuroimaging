#!/usr/bin/env python

#SBATCH --partition IB_44C_512G
#SBATCH --account iacc_madlab
#SBATCH --qos pq_madlab
#SBATCH -o /scratch/madlab/crash/mandy/vcat/norm_anat/ants_template/template_out
#SBATCH -e /scratch/madlab/crash/mandy/vcat/norm_anat/ants_template/template_err

import os
from glob import glob
from nipype.interfaces import ants
from nipype.interfaces.ants.legacy import buildtemplateparallel
from nipype.pipeline.engine import Workflow, Node, MapNode, JoinNode
from nipype.interfaces.utility import IdentityInterface, Function
from nipype.interfaces.io import DataGrabber, DataSink



fs_projdir = '/home/data/madlab/Mattfeld_vCAT/derivatives/freesurfer'
projdir = '/home/data/madlab/Mattfeld_vCAT'
crashdir = '/scratch/madlab/crash/mandy/vcat/norm_anat'
sinkdir = '/home/data/madlab/Mattfeld_vCAT/derivatives' #final output directory
workdir = os.path.join(crashdir,'ants_template') #workflow output directory
if not os.path.exists(workdir): #create workdir if nonexistent
    os.makedirs(workdir)


#define ANTs buildtemplate workflow and set base directory
buildtemplate_wf = Workflow(name = 'buildtemplate_wf')
buildtemplate_wf.base_dir = workdir


#node to execute buildtemplateparallel in Nipype
ants_template = Node(buildtemplateparallel(), 
                     name = 'ants_template')
ants_template.inputs.dimension = 3 # -d / image dimension
ants_template.inputs.bias_field_correction = True # -n / applies bias field correction to moving img
ants_template.inputs.gradient_step_size = 0.20 # -g / smaller magnitude means more cautious steps
ants_template.inputs.iteration_limit = 4 # -i / iterations of template construction
ants_template.inputs.max_iterations = [100, 70, 50, 20] # -m / maximum number of iterations, beginning with coarsest resolution
ants_template.inputs.out_prefix = 'vcat_antsTMPL_' # -o / prefix for all output files
ants_template.inputs.parallelization = 5 # -c / control for parallel processing [0 = serial, 1 = PBS, 2 = PEXEC, 3 = Apple XGrid, 4 = ??, 5 = SLURM]
ants_template.inputs.num_cores = 20 #-j / sets # of cpu cores
ants_template.inputs.rigid_body_registration = True # -r / registers input before creating template
ants_template.inputs.similarity_metric = 'CC' # -s / registration [CC = cross corr, MI = mutual info, PR = prob mapping, MSQ = mean square diff]
ants_template.inputs.transformation_model = 'GR' # -t / registration [EL = elastic trans, SY = SyN w/ time, S2 = 2 time, GR = greedy SyN, EX = expon, DD = diff demons]
ants_template.inputs.use_first_as_target = False # use first volume as target of all inputs - if False, unbiased average image is used to start
ants_template.inputs.in_files = glob(os.path.join(projdir, 'derivatives/vcat_firstpass_flirt/vcat_skullstrip_sub-0*.nii.gz')) #files used to generate template


#results to designated folder
datasink = Node(DataSink(), 
                name = 'datasink')
datasink.inputs.base_directory = os.path.join(sinkdir, 'vcat_T1_template')
buildtemplate_wf.connect(ants_template, 'final_template_file', datasink, 'PrimaryTemplate')


#execute workflow
buildtemplate_wf.run(plugin='SLURM', plugin_args={'sbatch_args': ('--partition IB_40C_512G --qos pq_madlab --account iacc_madlab -N 2 -n 20'), 'overwrite':True})

