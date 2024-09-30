#!/usr/bin/env python

"""
=========================================================
WMAZE_fMRI: Fixed Before Conditional Trials -- Model RSA4
=========================================================
First level workflow for UM GE 750 RSA1 task data

- Model RSA4
  - Use FSL ROI to recreate EPI data, removing last 3 volumes
  - Removed last 3 trials before EV creation
  - Pattern similarity analysis using realigned data
  - Parametric modulation analysis
  - EV directory (Model RSA4) --- /home/data/madlab/Mattfeld_WMAZE/sourcedata/behav/WMAZE_001/model_RSA4
- python RSA4_lvl1.py -s WMAZE_001
                       	   -o /home/data/madlab/Mattfeld_WMAZE/dset/analyses/model_RSA4/lvl1/cond
                      	   -w /scratch/data/crash/mandy/learning/RSA4/cond/lvl1
Note: DOF file is writing out in numpy hexadecimal format
      Example: 0x1.64p+7 
      print((1 + 6./0x10 + 4./0x100) * 2**7) = 178
"""
import os
from nipype.pipeline.engine import Workflow, Node, MapNode
from nipype.interfaces.utility import IdentityInterface, Function, Merge
from nipype.interfaces.io import DataGrabber, DataSink
from nipype.algorithms.modelgen import SpecifyModel
from nipype.interfaces.fsl.model import Level1Design, FEATModel, FILMGLS, ContrastMgr 
from nipype.interfaces.fsl.utils import ImageMaths, ExtractROI


###################
#### Functions ####
###################

#grab the first dimension of an array/matrix
pop_lambda = lambda x : x[0]

def subjectinfo(subject_id): #function to gather contrast info into Bunch
    import os
    from nipype.interfaces.base import Bunch
    from copy import deepcopy
    import numpy as np
    base_proj_dir = '/home/data/madlab/Mattfeld_WMAZE/behav'    
    output = [] #contains Bunch for each run (index 1-6)

    
    for curr_run in range(1,7): #iterate through runs
        names = []
        onsets = []
        durations = []
        amplitudes = []
 
        data_learn_event = np.genfromtxt(base_proj_dir + '/{0}/model_RSA4/run{1}_learn_event.txt'.format(subject_id,curr_run),dtype = str) #load event data
        data_nonlearn_event = np.genfromtxt(base_proj_dir + '/{0}/model_RSA4/run{1}_nonlearn_event.txt'.format(subject_id,curr_run),dtype = str)
        data_learn_mod = np.genfromtxt(base_proj_dir + '/{0}/model_RSA4/run{1}_learn_mod.txt'.format(subject_id,curr_run),dtype = str) #load pmod data
        data_nonlearn_mod = np.genfromtxt(base_proj_dir + '/{0}/model_RSA4/run{1}_nonlearn_mod.txt'.format(subject_id,curr_run),dtype = str)
        data_remaining = np.genfromtxt(base_proj_dir + '/{0}/model_RSA4/run{1}_remaining.txt'.format(subject_id,curr_run),dtype = str) 

        sequence = ['event', 'mod']
        for curr_type in sequence:
            learn_array_name = eval('data_learn_{0}'.format(curr_type))
            nonlearn_array_name = eval('data_nonlearn_{0}'.format(curr_type))
            if nonlearn_array_name.size > 0: #more than one nonlearn trial
                curr_names = ['learn_{0}'.format(curr_type), 'nonlearn_{0}'.format(curr_type)]
                curr_learn_onsets = map(float, learn_array_name[:,0]) #map funct converts all onset values from 1st column of EV file to float
                curr_learn_durations = map(float, learn_array_name[:,1]) #map funct converts all duration values from 2nd column
                curr_learn_amplitudes = map(float, learn_array_name[:,2]) #map funct converts all amplitude values from 3rd column
                if nonlearn_array_name.size == 3: #only one nonlearn trial  
                    curr_nonlearn_onsets = [float(nonlearn_array_name[0])] #converts single trial onset value from 1st column of EV file to float
                    curr_nonlearn_durations = [float(nonlearn_array_name[1])] #converts single trial duration value from 2nd column
                    curr_nonlearn_amplitudes = [float(nonlearn_array_name[2])] #converts single trial amplitude value from 3rd column
                else: #more than one nonlearning trial
                    curr_nonlearn_onsets = map(float, nonlearn_array_name[:,0])
                    curr_nonlearn_durations = map(float, nonlearn_array_name[:,1])
                    curr_nonlearn_amplitudes = map(float, nonlearn_array_name[:,2])
                curr_onsets = [curr_learn_onsets, curr_nonlearn_onsets] #add learn and nonlearn onset arrays into larger list
                curr_durations = [curr_learn_durations, curr_nonlearn_durations] #add learn and nonlearn duration arrays into larger list
                curr_amplitudes = [curr_learn_amplitudes, curr_nonlearn_amplitudes] #add learn and nonlearn amplitude arrays into larger list
            else: #no nonlearn trials
                curr_names = ['learn_{0}'.format(curr_type)]
                curr_learn_onsets = map(float, learn_array_name[:,0])
                curr_learn_durations = map(float, learn_array_name[:,1])
                curr_learn_amplitudes = map(float, learn_array_name[:,2])
                curr_onsets = [curr_learn_onsets]
                curr_durations = [curr_learn_durations]
                curr_amplitudes = [curr_learn_amplitudes]            
            names.append(curr_names) #append current run contrast names to current Bunch
            onsets.append(curr_onsets) #append current run contrast onsets to current Bunch
            durations.append(curr_durations) #append run contrast durations to current Bunch
            amplitudes.append(curr_amplitudes) #append run contrast amplitudes to current Bunch
  
        #remaining trials
        curr_names = ['remaining']
        curr_onsets = map(float, data_remaining[:,0])
        curr_durations = map(float, data_remaining[:,1])
        curr_amplitudes = map(float, data_remaining[:,2])
        curr_onsets = [curr_onsets]
        curr_durations = [curr_durations]
        curr_amplitudes = [curr_amplitudes]          
        names.append(curr_names)  
        onsets.append(curr_onsets)
        durations.append(curr_durations)
        amplitudes.append(curr_amplitudes)        
                   
        
        if any(isinstance(el, list) for el in names): #if any element in names is a list instead of a single value, for those elements       
            names = [el for sublist in names for el in sublist] #unpacks subarrays into one mega array! 
        if any(isinstance(el, list) for el in onsets):
            onsets = [el_o for sublist_o in onsets for el_o in sublist_o]
        if any(isinstance(el, list) for el in durations):
            durations = [el_d for sublist_d in durations for el_d in sublist_d]
        if any(isinstance(el, list) for el in amplitudes):
            amplitudes = [el_a for sublist_a in amplitudes for el_a in sublist_a]
      
        output.insert(curr_run,  #insert contents of each run at index of curr_run (1-6) 
                      Bunch(conditions = names,
                            onsets = deepcopy(onsets),
                            durations = deepcopy(durations),
                            amplitudes = deepcopy(amplitudes),
                            tmod = None, pmod = None,
                            regressor_names = None, regressors = None))
    return output

#obtain and create contrasts *flexibly* if there are not enough incorrect trials
def get_contrasts(subject_id, info): #receives Bunch from subject_info 
    contrasts = []
    for i, j in enumerate(info): #iterate through the six run Bunches
        curr_run_contrasts = []
        cont_all = ['AllVsBase', 'T', j.conditions, [1. / len(j.conditions)] * len(j.conditions)] #all types vs. baseline
        curr_run_contrasts.append(cont_all) #append to current run contrast info
        for curr_cond in j.conditions: #iterate through contrast types
            curr_cont = [curr_cond, 'T', [curr_cond], [1]] #iterable contrast definition list
            curr_run_contrasts.append(curr_cont) #append to current run   
        if 'learn' in j.conditions and 'nonlearn' in j.conditions: #if a subject has both learn and nonlearn trial types
            cont_learn_vs_nonlearn = ['mod_learn_minus_nonlearn', 'T', ['learn_mod','nonlearn_mod'], [1, -1]] #contrast info
            cont_nonlearn_vs_learn = ['mod_nonlearn_minus_learn', 'T', ['learn_mod','nonlearn_mod'], [-1, 1]]
            curr_run_contrasts.append(cont_learn_vs_nonlearn) #append to current run
            curr_run_contrasts.append(cont_nonlearn_vs_learn)
        contrasts.append(curr_run_contrasts) #append to all subject contrasts
    return contrasts
                                                                            

#set naming convention for output types
def get_subs(cons):
    subs = []
    for run_cons in cons: #iterate through runs
        run_subs = []
        for i, con in enumerate(run_cons): #iterate through contrasts in run
            run_subs.append(('cope%d.'%(i + 1), 'cope%02d_%s.'%(i + 1, con[0]))) #set cope number and name
            run_subs.append(('varcope%d.'%(i + 1), 'varcope%02d_%s.'%(i + 1, con[0])))
            run_subs.append(('zstat%d.'%(i + 1), 'zstat%02d_%s.'%(i + 1, con[0])))
            run_subs.append(('tstat%d.'%(i + 1), 'tstat%02d_%s.'%(i + 1, con[0])))
        subs.append(run_subs)        
    return subs


#extracting motion parameters from noise files
def motion_noise(subjinfo, files): #received from subject_info and motion noise files 
    import numpy as np
    motion_noise_params = []
    motion_noi_par_names = []
    if not isinstance(files, list): #converts to list if not so already
        files = [files]
    if not isinstance(subjinfo, list):
        subjinfo = [subjinfo]
    for j, i in enumerate(files): #iterates through motion noise files
        curr_mot_noi_par_names = ['Pitch (rad)', 'Roll (rad)', 'Yaw (rad)', 'Tx (mm)', 'Ty (mm)', 'Tz (mm)', #establishes names for motion regressors
                                  'Pitch_1d', 'Roll_1d', 'Yaw_1d', 'Tx_1d', 'Ty_1d', 'Tz_1d',
                                  'Norm (mm)', 'LG_1stOrd', 'LG_2ndOrd', 'LG_3rdOrd', 'LG_4thOrd']
        a = np.genfromtxt(i) #open current motion noise file
        motion_noise_params.append([[]] * a.shape[1]) #append number of empty sublists to match number of columns in motion noise file
        if a.shape[1] > 17: #if number of columns exceed established names
            for num_out in range(a.shape[1] - 17): #for excess motion noise columns
                out_name = 'out_{0}'.format(num_out + 1) #name the additional column
                curr_mot_noi_par_names.append(out_name) #append additional name to established names
        for z in range(a.shape[1]): #for each motion noise column
            motion_noise_params[j][z] = a[:, z].tolist() #set column for current noise parameter in current file
        motion_noi_par_names.append(curr_mot_noi_par_names) 
    for j, i in enumerate(subjinfo): #iterate through Bunches
        if i.regressor_names == None: #if current Bunch regressor is emptyy 
            i.regressor_names = [] #make empty list
        if i.regressors == None: 
            i.regressors = []
        for j2, i2 in enumerate(motion_noise_params[j]): #for each motion noise parameter
            i.regressor_names.append(motion_noi_par_names[j][j2]) #append current name
            i.regressors.append(i2) #append current regressor           
    return subjinfo


#####################################
## Function defining lvl1 analysis ##
#####################################


def firstlevel_wf(subject_id,sink_directory, name = 'wmaze_frstlvl_wf'):
    frstlvl_wf = Workflow(name = 'frstlvl_wf')
  

    #function with name, onset, duration, and amplitude info 
    subject_info = Node(Function(input_names = ['subject_id'], output_names = ['output'],
                                 function = subjectinfo),
                        name = 'subject_info')
    subject_info.inputs.ignore_exception = False
    subject_info.inputs.subject_id = subject_id


    #function to define contrasts
    getcontrasts = Node(Function(input_names = ['subject_id', 'info'], output_names = ['contrasts'],
                                 function = get_contrasts),
                        name = 'getcontrasts')
    getcontrasts.inputs.ignore_exception = False
    getcontrasts.inputs.subject_id = subject_id
    frstlvl_wf.connect(subject_info, 'output', getcontrasts, 'info')

    
    #function to substitute names of output folders and files
    getsubs = Node(Function(input_names = ['cons'], output_names = ['subs'],
                            function = get_subs),
                   name = 'getsubs')
    getsubs.inputs.ignore_exception = False
    getsubs.inputs.subject_id = subject_id
    frstlvl_wf.connect(subject_info, 'output', getsubs, 'info')
    frstlvl_wf.connect(getcontrasts, 'contrasts', getsubs, 'cons')


    #dictionary holding wildcards used in datasource
    info = dict(task_mri_files = [['subject_id']], 
                motion_noise_files = [['subject_id']])
    

    #get task_mri and motion-noise files
    datasource = Node(DataGrabber(infields = ['subject_id'], outfields = info.keys()), 
                      name = 'datasource')
    datasource.inputs.template = '*'
    datasource.inputs.subject_id = subject_id
    datasource.inputs.base_directory = os.path.abspath('/home/data/madlab/Mattfeld_WMAZE/derivatives/preproc/')
    datasource.inputs.field_template = dict(task_mri_files = '%s/func/realigned/*wmaze*.nii.gz',
                                            motion_noise_files = '%s/noise/filter_regressor??.txt')
    datasource.inputs.template_args = info
    datasource.inputs.sort_filelist = True
    datasource.inputs.ignore_exception = False
    datasource.inputs.raise_on_empty = True


    #function to remove last three volumes from func data                                
    fslroi_epi = MapNode(ExtractROI(t_min = 0, t_size = 197), #start from the first volume and end on the -3 volume
                         iterfield = ['in_file'],
                         name = 'fslroi_epi')
    fslroi_epi.output_type = 'NIFTI_GZ'
    fslroi_epi.terminal_output = 'stream'
    frstlvl_wf.connect(datasource, 'task_mri_files', fslroi_epi, 'in_file')


    #function to modify the motion and noise files to single regressors
    motionnoise = Node(Function(input_names = ['subjinfo', 'files'], output_names = ['subjinfo'],
                                function = motion_noise),
                       name = 'motionnoise')
    motionnoise.inputs.ignore_exception = False
    frstlvl_wf.connect(subject_info, 'output', motionnoise, 'subjinfo')
    frstlvl_wf.connect(datasource, 'motion_noise_files', motionnoise, 'files')


    #makes model specification compatible with FSL (requires subject info be in Bunch format)
    specify_model = Node(SpecifyModel(), 
                         name = 'specify_model')
    specify_model.inputs.high_pass_filter_cutoff = -1.0
    specify_model.inputs.ignore_exception = False
    specify_model.inputs.input_units = 'secs'
    specify_model.inputs.time_repetition = 2.0
    frstlvl_wf.connect(fslroi_epi, 'roi_file', specify_model, 'functional_runs') 
    frstlvl_wf.connect(motionnoise, 'subjinfo', specify_model, 'subject_info')

    
    modelfit_inputspec = Node(IdentityInterface(fields = ['session_info', 'interscan_interval', 'contrasts', 'film_threshold', 
                                                          'functional_data', 'bases', 'model_serial_correlations'], 
                                                mandatory_inputs = True),
                              name = 'modelfit_inputspec')
    modelfit_inputspec.inputs.bases = {'dgamma':{'derivs': False}}
    modelfit_inputspec.inputs.film_threshold = 0.0
    modelfit_inputspec.inputs.interscan_interval = 2.0
    modelfit_inputspec.inputs.model_serial_correlations = True
    frstlvl_wf.connect(fslroi_epi, 'roi_file', modelfit_inputspec, 'functional_data')
    frstlvl_wf.connect(getcontrasts, 'contrasts', modelfit_inputspec, 'contrasts')
    frstlvl_wf.connect(specify_model, 'session_info', modelfit_inputspec, 'session_info')
   

    #first level SPM design matrix to demonstrate contrasts and motion/noise regressors
    level1_design = MapNode(Level1Design(), iterfield = ['contrasts', 'session_info'],
                            name = 'level1_design')
    level1_design.inputs.ignore_exception = False
    frstlvl_wf.connect(modelfit_inputspec, 'interscan_interval', level1_design, 'interscan_interval')
    frstlvl_wf.connect(modelfit_inputspec, 'session_info', level1_design, 'session_info')
    frstlvl_wf.connect(modelfit_inputspec, 'contrasts', level1_design, 'contrasts')
    frstlvl_wf.connect(modelfit_inputspec, 'bases', level1_design, 'bases')
    frstlvl_wf.connect(modelfit_inputspec, 'model_serial_correlations', level1_design, 'model_serial_correlations')
    

    #generate design.mat file for each run
    generate_model = MapNode(FEATModel(), iterfield = ['fsf_file', 'ev_files'],
                             name = 'generate_model') 
    generate_model.inputs.environ = {'FSLOUTPUTTYPE': 'NIFTI_GZ'}
    generate_model.inputs.ignore_exception = False
    generate_model.inputs.output_type = 'NIFTI_GZ'
    generate_model.inputs.terminal_output = 'stream'
    frstlvl_wf.connect(level1_design, 'fsf_files', generate_model, 'fsf_file')
    frstlvl_wf.connect(level1_design, 'ev_files', generate_model, 'ev_files')
  

    #estimate model using FILMGLS -- fits design matrix to voxel timeseries
    estimate_model = MapNode(FILMGLS(), iterfield = ['design_file', 'in_file', 'tcon_file'],
                             name = 'estimate_model')
    estimate_model.inputs.environ = {'FSLOUTPUTTYPE': 'NIFTI_GZ'}
    estimate_model.inputs.ignore_exception = False
    estimate_model.inputs.mask_size = 5
    estimate_model.inputs.output_type = 'NIFTI_GZ'
    estimate_model.inputs.results_dir = 'results'
    estimate_model.inputs.smooth_autocorr = True
    estimate_model.inputs.terminal_output = 'stream'
    frstlvl_wf.connect(modelfit_inputspec, 'film_threshold', estimate_model, 'threshold')
    frstlvl_wf.connect(modelfit_inputspec, 'functional_data', estimate_model, 'in_file')
    frstlvl_wf.connect(generate_model, 'design_file', estimate_model, 'design_file')
    frstlvl_wf.connect(generate_model, 'con_file', estimate_model, 'tcon_file')



    #merge contrasts (necessary for fsl 5.0.7 and greater)
    merge_contrasts = MapNode(Merge(2), iterfield = ['in1'], 
                              name = 'merge_contrasts')
    frstlvl_wf.connect(estimate_model, 'zstats', merge_contrasts, 'in1')



    #transform z2pval
    z2pval = MapNode(ImageMaths(), iterfield = ['in_file'], 
                     name='z2pval')
    z2pval.inputs.environ = {'FSLOUTPUTTYPE': 'NIFTI_GZ'}
    z2pval.inputs.ignore_exception = False
    z2pval.inputs.op_string = '-ztop'
    z2pval.inputs.output_type = 'NIFTI_GZ'
    z2pval.inputs.suffix = '_pval'
    z2pval.inputs.terminal_output = 'stream'
    frstlvl_wf.connect(merge_contrasts, ('out', pop_lambda), z2pval, 'in_file')



    #IdentityInterface() to receive info from estimate_model, merge_contrasts, z2pval, and generate_model
    modelfit_outputspec = Node(IdentityInterface(fields = ['copes', 'varcopes', 'dof_file', 'pfiles', 'parameter_estimates', 'zstats', 
                                                           'design_image', 'design_file', 'design_cov', 'sigmasquareds'],
                                                 mandatory_inputs = True),
                               name = 'modelfit_outputspec')
    frstlvl_wf.connect(estimate_model, 'copes', modelfit_outputspec, 'copes')
    frstlvl_wf.connect(estimate_model, 'varcopes', modelfit_outputspec, 'varcopes')
    frstlvl_wf.connect(merge_contrasts, 'out', modelfit_outputspec, 'zstats')
    frstlvl_wf.connect(z2pval, 'out_file', modelfit_outputspec, 'pfiles')
    frstlvl_wf.connect(generate_model, 'design_image', modelfit_outputspec, 'design_image')
    frstlvl_wf.connect(generate_model, 'design_file', modelfit_outputspec, 'design_file')
    frstlvl_wf.connect(generate_model, 'design_cov', modelfit_outputspec, 'design_cov')
    frstlvl_wf.connect(estimate_model, 'param_estimates', modelfit_outputspec, 'parameter_estimates')
    frstlvl_wf.connect(estimate_model, 'dof_file', modelfit_outputspec, 'dof_file')
    frstlvl_wf.connect(estimate_model, 'sigmasquareds', modelfit_outputspec, 'sigmasquareds')


    #save output from multiple points in pipeline
    sinkd = MapNode(DataSink(), iterfield = ['substitutions', 'modelfit.contrasts.@copes', 'modelfit.contrasts.@varcopes',
                                             'modelfit.estimates', 'modelfit.contrasts.@zstats'],
                    name = 'sinkd')
    sinkd.inputs.base_directory = sink_directory 
    sinkd.inputs.container = subject_id
    frstlvl_wf.connect(getsubs, 'subs', sinkd, 'substitutions')
    frstlvl_wf.connect(modelfit_outputspec, 'parameter_estimates', sinkd, 'modelfit.estimates')
    frstlvl_wf.connect(modelfit_outputspec, 'sigmasquareds', sinkd, 'modelfit.estimates.@sigsq')
    frstlvl_wf.connect(modelfit_outputspec, 'dof_file', sinkd, 'modelfit.dofs')
    frstlvl_wf.connect(modelfit_outputspec, 'copes', sinkd, 'modelfit.contrasts.@copes')
    frstlvl_wf.connect(modelfit_outputspec, 'varcopes', sinkd, 'modelfit.contrasts.@varcopes')
    frstlvl_wf.connect(modelfit_outputspec, 'zstats', sinkd, 'modelfit.contrasts.@zstats')
    frstlvl_wf.connect(modelfit_outputspec, 'design_image', sinkd, 'modelfit.design')
    frstlvl_wf.connect(modelfit_outputspec, 'design_cov', sinkd, 'modelfit.design.@cov')
    frstlvl_wf.connect(modelfit_outputspec, 'design_file', sinkd, 'modelfit.design.@matrix')
    frstlvl_wf.connect(modelfit_outputspec, 'pfiles', sinkd, 'modelfit.contrasts.@pstats')

    return frstlvl_wf


###############################
## Creates the full workflow ##
###############################

def create_frstlvl_workflow(args, name = 'wmaze_MR_frstlvl'): #func to execute workflow
    kwargs = dict(subject_id = args.subject_id, sink_directory = os.path.abspath(args.out_dir), name = name)
    frstlvl_workflow = firstlevel_wf(**kwargs) #first level workflow with kwargs as parameters
    return frstlvl_workflow

if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser(description = __doc__)
    parser.add_argument("-s", "--subject_id", dest = "subject_id", help = "Current subject id", required = True)
    parser.add_argument("-o", "--output_dir", dest = "out_dir", help = "Output directory base")
    parser.add_argument("-w", "--work_dir", dest = "work_dir", help = "Working directory base")
    args = parser.parse_args()

    wf = create_frstlvl_workflow(args)
    if args.work_dir:
        work_dir = os.path.abspath(args.work_dir)    
    else:
        work_dir = os.getcwd()

    wf.config['execution']['crashdump_dir'] = '/scratch/madlab/crash/mandy/Mattfeld_WMAZE/RSA4/' #define location to write crashfiles
    wf.base_dir = work_dir + '/' + args.subject_id #dir where workflow will be written
    wf.run(plugin='SLURM', plugin_args={'sbatch_args': ('-p centos7_16C_128G --account iacc_madlab --qos pq_madlab -N 1 -n 1'), 'overwrite': True})
