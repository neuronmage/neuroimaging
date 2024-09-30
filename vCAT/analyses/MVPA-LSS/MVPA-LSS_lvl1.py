#!/usr/bin/env python

"""
======================================
MVPA_fMRI: LSS Task  -- Model MVPA-LSS
======================================
First level workflow for Seimens vCAT MVPA task data

- vCAT MVPA-LSS
  - Use FSL ROI to recreate EPI data
  - EV directory --- /home/data/madlab/Mattfeld_vCAT/behav/vCAT_0??/MVPA-LSS

- python MVPA-LSS_lvl1.py  -s vCAT_0??
                           -o /home/data/madlab/Mattfeld_vCAT/derivatives/MVPA/lvl1
                           -w /scratch/madlab/crash/mandy/vcat/MVPA/lvl1

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

# Grab the first dimension of an array/matrix
pop_lambda = lambda x : x[0]

def subjectinfo(subject_id):
    import os
    from nipype.interfaces.base import Bunch
    from copy import deepcopy
    import numpy as np
    proj_dir = '/home/data/madlab/Mattfeld_vCAT/behav/'
    output = [] 

    
    model_counter = 0 #MODEL COUNTER WILL UPDATE FOR EACH TRIAL OF A GIVEN CONDITION
    for curr_run in range(1,5):
        data_FX_BL_face = np.genfromtxt(proj_dir + '/{0}/MVPA-LSS/run{1}_fixed_BL_face.txt'.format(subject_id,curr_run),dtype = str)
        data_FX_BL_scene = np.genfromtxt(proj_dir + '/{0}/MVPA-LSS/run{1}_fixed_BL_scene.txt'.format(subject_id,curr_run),dtype = str)
        data_all_remaining = np.genfromtxt(proj_dir + '/{0}/MVPA-LSS/run{1}_remaining.txt'.format(subject_id,curr_run),dtype = str)
        
        # CONSOLIDATE DATA INTO A DICTIONARY FOR ITERATION
        orig_FX_BL_data = {'FX_BL_face': data_FX_BL_face,
                           'FX_BL_scene': data_FX_BL_scene}

        # ITERATE OVER KEYS OF THE DICTIONARY TO ISOLATE CONDITIONS OF INTEREST
        for curr_key in orig_FX_BL_data.keys():
            # ESTABLISH TRIAL COUNTER FOR NAMING OF REGRESSORS
            trial_counter = 1
            # ISOLATE CURRENT CONDITION DATA USING POP FUNCTION
            # DICTIONARY WILL NO LONGER HAVE THAT KEY
            # I USE THAT FUNCTIONALITY TO ESTABLISH PENDING KEYS (NOT YET ITERATED OVER)
            copy_FX_BL_data = dict(orig_FX_BL_data)
            curr_condition_data = copy_FX_BL_data.pop(curr_key)

            if curr_condition_data.size == 3: # ONLY ONE EVENT OF THIS CONDITION DURING THIS RUN
                names = [curr_key + '_run%d_onset%0.2f_trl%d' %(curr_run, float(curr_condition_data[0]), trial_counter)]
                onsets = [[float(curr_condition_data[0])]]
                durations = [[float(curr_condition_data[1])]]
                amplitudes = [[float(curr_condition_data[2])]]
                # DEAL WITH REMAINING DATA THAT HASN'T BEEN ITERATED THROUGH YET (AKA PENDING)
                for pending_key in copy_FX_BL_data.keys():
                    names.append(pending_key)
                    pending_data = copy_FX_BL_data[pending_key]
                    if pending_data.size == 3: #ONLY ONE EVENT OF THIS CONDITION
                        onsets.append([float(pending_data[0])])
                        durations.append([float(pending_data[1])])
                        amplitudes.append([float(pending_data[2])])
                    else:
                        onsets.append(list(map(float,pending_data[:,0])))
                        durations.append(list(map(float,pending_data[:,1])))
                        amplitudes.append(list(map(float,pending_data[:,2])))
                # INSERT ALL REMAINING EV INTO MODEL
                names.append('all_remaining')
                onsets.append(list(map(float, data_all_remaining[:,0])))
                durations.append(list(map(float, data_all_remaining[:,1])))
                amplitudes.append(list(map(float, data_all_remaining[:,2])))

                # UPDATE TRIAL COUNTER
                trial_counter = trial_counter + 1

                # Insert contents of each run at index of model_counter
                output.insert(model_counter,
                              Bunch(conditions = names,
                                    onsets = deepcopy(onsets),
                                    durations = deepcopy(durations),
                                    amplitudes = deepcopy(amplitudes),
                                    tmod = None,
                                    pmod = None,
                                    regressor_names = None,
                                    regressors = None))

                # UPDATE MODEL COUNTER
                model_counter = model_counter + 1
            else: # THERE IS MORE THAN ONE EVENT OF THIS CONDITION DURING THIS RUN
                # ITERATE OVER NUMBER OF TRIALS WITHIN THAT CONDITION
                for curr_cond_trl in range(len(curr_condition_data)):
                    # ESTABLISH LISTS FOR NAMES, ONSETS, DURATIONS, AND AMPLITUDES FOR ALL MODELS
                    # WE WILL HAVE AS MANY MODELS AS TRIALS ACROSS RUNS FOR DIFFERENT CONDITIONS
                    names = []
                    onsets = []
                    durations = []
                    amplitudes = []
                    curr_cond_trl_name = curr_key + '_run%d_onset%0.2f_trl%d' %(curr_run, float(curr_condition_data[curr_cond_trl][0]), trial_counter)
                    curr_cond_trl_onset = [float(curr_condition_data[curr_cond_trl][0])]
                    curr_cond_trl_dur = [float(curr_condition_data[curr_cond_trl][1])]
                    curr_cond_trl_amp = [float(curr_condition_data[curr_cond_trl][2])]
                    
                    names.append(curr_cond_trl_name)
                    onsets.append(curr_cond_trl_onset)
                    durations.append(curr_cond_trl_dur)
                    amplitudes.append(curr_cond_trl_amp)
                
                    # ISOLATE REMAINING TRIALS FOR CURRENT CONDITION USING NUMPY DELETE FUNCTION
                    # THIS FUNCTION WILL NOT MODIFY ORIGINAL VARIABLE LIKE POP DOES ABOVE
                    curr_cond_remaining_data = np.delete(curr_condition_data, curr_cond_trl, 0)
                    curr_cond_remaining_name = curr_key + '_allbut_run%d_trl%d' %(curr_run, trial_counter)
                    curr_cond_remaining_onsets = list(map(float, curr_cond_remaining_data[:,0]))
                    curr_cond_remaining_durs = list(map(float, curr_cond_remaining_data[:,1]))
                    curr_cond_remaining_amps = list(map(float, curr_cond_remaining_data[:,2]))
                    
                    names.append(curr_cond_remaining_name)
                    onsets.append(curr_cond_remaining_onsets)
                    durations.append(curr_cond_remaining_durs)
                    amplitudes.append(curr_cond_remaining_amps)
                
                    # DEAL WITH PENDING DATA THAT HASN'T BEEN ITERATED THROUGH YET
                    # THIS IS WHERE THAT POP FUNCTION ABOVE CAME IN HANDY
                    for pending_key in copy_FX_BL_data.keys():
                        names.append(pending_key)
                        pending_data = copy_FX_BL_data[pending_key]
                        if pending_data.size == 3: #ONLY ONE EVENT OF THIS CONDITION
                            onsets.append([float(pending_data[0])])
                            durations.append([float(pending_data[1])])
                            amplitudes.append([float(pending_data[2])])
                        else:
                            onsets.append(list(map(float,pending_data[:,0])))
                            durations.append(list(map(float,pending_data[:,1])))
                            amplitudes.append(list(map(float,pending_data[:,2])))
               
                    # INSERT ALL REAMINING EV INTO MODEL
                    names.append('all_remaining')
                    onsets.append(list(map(float, data_all_remaining[:,0])))
                    durations.append(list(map(float, data_all_remaining[:,1])))
                    amplitudes.append(list(map(float, data_all_remaining[:,2])))

                    # UPDATE TRIAL COUNTER
                    trial_counter = trial_counter + 1

                    # Insert contents of each run at index of model_counter
                    output.insert(model_counter,
                                  Bunch(conditions = names,
                                        onsets = deepcopy(onsets),
                                        durations = deepcopy(durations),
                                        amplitudes = deepcopy(amplitudes),
                                        tmod = None,
                                        pmod = None,
                                        regressor_names = None,
                                        regressors = None))

                    # UPDATE MODEL COUNTER
                    model_counter = model_counter + 1
    return output



# Function to obtain and create contrasts *flexibly* in case there are not enough incorrect trials
def get_contrasts(subject_id, info):
    contrasts = []
    for i, j in enumerate(info):
        curr_run_contrasts = []
        for curr_cond in j.conditions:
            curr_cont = [curr_cond, 'T', [curr_cond], [1]]
            curr_run_contrasts.append(curr_cont)            
        contrasts.append(curr_run_contrasts)
    return contrasts


# Function to name output types
def get_subs(cons):
    subs = []
    for run_cons in cons:
        run_subs = []
        for i, con in enumerate(run_cons):
            run_subs.append(('cope%d.'%(i + 1), 'cope%02d_%s.'%(i + 1, con[0])))
            run_subs.append(('varcope%d.'%(i + 1), 'varcope%02d_%s.'%(i + 1, con[0])))
            run_subs.append(('zstat%d.'%(i + 1), 'zstat%02d_%s.'%(i + 1, con[0])))
            run_subs.append(('tstat%d.'%(i + 1), 'tstat%02d_%s.'%(i + 1, con[0])))
        subs.append(run_subs)       
    return subs


# Function to extract motion parameters from noise files
def motion_noise(subjinfo, files):
    import numpy as np
    motion_noise_params = []
    motion_noi_par_names = []
    if not isinstance(files, list):
        files = [files]
    if not isinstance(subjinfo, list):
        subjinfo = [subjinfo]
    for j,i in enumerate(files):
        curr_mot_noi_par_names = ['Pitch (rad)', 'Roll (rad)', 'Yaw (rad)', 'Tx (mm)', 'Ty (mm)', 'Tz (mm)',
                                  'Pitch_1d', 'Roll_1d', 'Yaw_1d', 'Tx_1d', 'Ty_1d', 'Tz_1d',
                                  'Norm (mm)', 'LG_1stOrd', 'LG_2ndOrd', 'LG_3rdOrd', 'LG_4thOrd']
        a = np.genfromtxt(i)
        motion_noise_params.append([[]] * a.shape[1])
        if a.shape[1] > 17:
            for num_out in range(a.shape[1] - 17):
                out_name = 'out_{0}'.format(num_out + 1)
                curr_mot_noi_par_names.append(out_name)
        for z in range(a.shape[1]):
            motion_noise_params[j][z] = a[:, z].tolist()
        motion_noi_par_names.append(curr_mot_noi_par_names)   
    for j,i in enumerate(subjinfo):
        if i.regressor_names == None: 
            i.regressor_names = []
        if i.regressors == None: 
            i.regressors = []
        if 'run1' in i.conditions[0]:
            curr_run = 0
        elif 'run2' in i.conditions[0]:
            curr_run = 1
        elif 'run3' in i.conditions[0]:
            curr_run = 2
        elif 'run4' in i.conditions[0]:
            curr_run = 3
        for j3, i3 in enumerate(motion_noise_params[curr_run]):
            i.regressor_names.append(motion_noi_par_names[curr_run][j3])
            i.regressors.append(i3)           
    return subjinfo


def expand_files(subjinfo, in_files):
    if not isinstance(in_files, list):
        in_files = [in_files]
    files_expanded = []
    for j,i in enumerate(subjinfo):
        if 'run1' in i.conditions[0]:
            files_expanded.append(in_files[0])
        elif 'run2' in i.conditions[0]:
            files_expanded.append(in_files[1])
        elif 'run3' in i.conditions[0]:
            files_expanded.append(in_files[2])
        elif 'run4' in i.conditions[0]:
            files_expanded.append(in_files[3])
    return files_expanded


###################################
## Function for 1st lvl analysis ##
###################################
def firstlevel_wf(subject_id,
                  sink_directory,
                  name = 'wmaze_frstlvl_wf'):
    frstlvl_wf = Workflow(name = 'frstlvl_wf')
    
    
    #wildcard used in datasource
    info = dict(task_mri_files = [['subject_id', 'vcat']],
                motion_noise_files = [['subject_id']])


    #calls subjectinfo function with name, onset, duration, and amplitude info 
    subject_info = Node(Function(input_names = ['subject_id'],
                                 output_names = ['output'],
                                 function = subjectinfo),
                        name = 'subject_info')
    subject_info.inputs.ignore_exception = False
    subject_info.inputs.subject_id = subject_id


    #define contrasts for experiment
    getcontrasts = Node(Function(input_names = ['subject_id', 'info'],
                                 output_names = ['contrasts'],
                                 function = get_contrasts),
                        name = 'getcontrasts')
    getcontrasts.inputs.ignore_exception = False
    getcontrasts.inputs.subject_id = subject_id
    frstlvl_wf.connect(subject_info, 'output', getcontrasts, 'info')


    #substitute names of folders and files created during pipeline
    getsubs = Node(Function(input_names = ['cons'],
                            output_names = ['subs'],
                            function = get_subs),
                   name = 'getsubs')
    getsubs.inputs.ignore_exception = False
    getsubs.inputs.subject_id = subject_id
    frstlvl_wf.connect(subject_info, 'output', getsubs, 'info')
    frstlvl_wf.connect(getcontrasts, 'contrasts', getsubs, 'cons')
  

    #get task_mri and motion-noise files
    datasource = Node(DataGrabber(infields = ['subject_id'], 
                                  outfields = list(info.keys())),
                      name = 'datasource')
    datasource.inputs.template = '*'
    datasource.inputs.subject_id = subject_id
    datasource.inputs.base_directory = os.path.abspath('/home/data/madlab/Mattfeld_vCAT/derivatives/preproc/')
    datasource.inputs.field_template = dict(task_mri_files = '%s/func/smoothed_fullspectrum/_maskfunc2*/*%s*.nii.gz', #func files
                                            motion_noise_files = '%s/noise/filter_regressor??.txt') #filter regressor noise files
    datasource.inputs.template_args = info
    datasource.inputs.sort_filelist = True
    datasource.inputs.ignore_exception = False
    datasource.inputs.raise_on_empty = True


    #modify motion and noise files to be single regressors
    motionnoise = Node(Function(input_names = ['subjinfo', 'files'],
                                output_names = ['subjinfo'],
                                function = motion_noise),
                       name = 'motionnoise')
    motionnoise.inputs.ignore_exception = False
    frstlvl_wf.connect(subject_info, 'output', motionnoise, 'subjinfo')
    frstlvl_wf.connect(datasource, 'motion_noise_files', motionnoise, 'files') 
 

    #expand task functional data
    expand_epi_files = Node(Function(input_names = ['subjinfo', 'in_files'],
                                     output_names = ['files_expanded'],
                                     function = expand_files),
                            name = 'expand_epi_files')
    expand_epi_files.inputs.ignore_exception = False
    frstlvl_wf.connect(motionnoise, 'subjinfo', expand_epi_files, 'subjinfo')
    frstlvl_wf.connect(datasource, 'task_mri_files', expand_epi_files, 'in_files') 


    #model specification compatible with spm/fsl designers
    #requires subjectinfo be received as a Bunch of a list of Bunch
    specify_model = Node(SpecifyModel(), 
                         name = 'specify_model')
    specify_model.inputs.high_pass_filter_cutoff = -1.0
    specify_model.inputs.input_units = 'secs'
    specify_model.inputs.time_repetition = 1.76
    frstlvl_wf.connect(expand_epi_files, 'files_expanded', specify_model, 'functional_runs')
    frstlvl_wf.connect(motionnoise, 'subjinfo', specify_model, 'subject_info')
  

    #basic interface class generates identity mappings
    modelfit_inputspec = Node(IdentityInterface(fields = ['session_info', 'interscan_interval', 'contrasts',
                                                          'film_threshold', 'functional_data', 'bases',
                                                          'model_serial_correlations'], 
                                                mandatory_inputs = True),
                              name = 'modelfit_inputspec')
    modelfit_inputspec.inputs.bases = {'dgamma':{'derivs': False}}
    modelfit_inputspec.inputs.film_threshold = 0.0
    modelfit_inputspec.inputs.interscan_interval = 1.76
    modelfit_inputspec.inputs.model_serial_correlations = True
    frstlvl_wf.connect(expand_epi_files, 'files_expanded', modelfit_inputspec, 'functional_data')
    frstlvl_wf.connect(getcontrasts, 'contrasts', modelfit_inputspec, 'contrasts')
    frstlvl_wf.connect(specify_model, 'session_info', modelfit_inputspec, 'session_info')
 
    
    #level1 design node
    level1_design = MapNode(Level1Design(),
                            iterfield = ['contrasts', 'session_info'],
                            name = 'level1_design')
    frstlvl_wf.connect(modelfit_inputspec, 'interscan_interval', level1_design, 'interscan_interval')
    frstlvl_wf.connect(modelfit_inputspec, 'session_info', level1_design, 'session_info')
    frstlvl_wf.connect(modelfit_inputspec, 'contrasts', level1_design, 'contrasts')
    frstlvl_wf.connect(modelfit_inputspec, 'bases', level1_design, 'bases')
    frstlvl_wf.connect(modelfit_inputspec, 'model_serial_correlations', level1_design, 'model_serial_correlations')
  

    #generate a model for each run as design.mat files
    generate_model = MapNode(FEATModel(),
                             iterfield = ['fsf_file', 'ev_files'],
                             name = 'generate_model') 
    generate_model.inputs.environ = {'FSLOUTPUTTYPE': 'NIFTI_GZ'}
    generate_model.inputs.output_type = 'NIFTI_GZ'
    frstlvl_wf.connect(level1_design, 'fsf_files', generate_model, 'fsf_file')
    frstlvl_wf.connect(level1_design, 'ev_files', generate_model, 'ev_files')
   

    #estimate model using FILMGLS -- fit a design matrix to a voxel timeseries
    estimate_model = MapNode(FILMGLS(),
                             iterfield = ['design_file', 'in_file', 'tcon_file'],
                             name = 'estimate_model')
    estimate_model.inputs.environ = {'FSLOUTPUTTYPE': 'NIFTI_GZ'}
    estimate_model.inputs.mask_size = 5
    estimate_model.inputs.output_type = 'NIFTI_GZ'
    estimate_model.inputs.results_dir = 'results'
    estimate_model.inputs.smooth_autocorr = True
    frstlvl_wf.connect(modelfit_inputspec, 'film_threshold', estimate_model, 'threshold')
    frstlvl_wf.connect(modelfit_inputspec, 'functional_data', estimate_model, 'in_file')
    frstlvl_wf.connect(generate_model, 'design_file', estimate_model, 'design_file') #ascii matrix for design    
    frstlvl_wf.connect(generate_model, 'con_file', estimate_model, 'tcon_file') #contrast vectors


    #merge contrasts - necessary for fsl 5.0.7 and greater
    merge_contrasts = MapNode(Merge(2), 
                              iterfield = ['in1'], 
                              name = 'merge_contrasts')
    frstlvl_wf.connect(estimate_model, 'zstats', merge_contrasts, 'in1')


    #transform z2pval
    z2pval = MapNode(ImageMaths(), 
                     iterfield = ['in_file'], 
                     name='z2pval')
    z2pval.inputs.environ = {'FSLOUTPUTTYPE': 'NIFTI_GZ'}
    z2pval.inputs.op_string = '-ztop'
    z2pval.inputs.output_type = 'NIFTI_GZ'
    z2pval.inputs.suffix = '_pval'
    frstlvl_wf.connect(merge_contrasts, ('out', pop_lambda), z2pval, 'in_file')


    #receive information from estimate_model, merge_contrasts, z2pval, generate_model, and estimate_model
    modelfit_outputspec = Node(IdentityInterface(fields = ['copes', 'varcopes', 'dof_file', 'pfiles', 'parameter_estimates', 
                                                           'zstats', 'design_image', 'design_file', 'design_cov', 'sigmasquareds'],
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
    
    
    #datasink node
    sinkd = MapNode(DataSink(),
                    iterfield = ['substitutions', 'modelfit.contrasts.@copes', 'modelfit.contrasts.@varcopes',
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

def create_frstlvl_workflow(args, name = 'wmaze_MR_frstlvl'):
    kwargs = dict(subject_id = args.subject_id, sink_directory = os.path.abspath(args.out_dir), name = name)
    frstlvl_workflow = firstlevel_wf(**kwargs)
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

    wf.config['execution']['crashdump_dir'] = '/scratch/madlab/crash/mandy/vcat/MVPA-LSS/lvl1'
    wf.base_dir = work_dir + '/' + args.subject_id
    wf.run(plugin='SLURM', plugin_args={'sbatch_args': ('-p IB_40C_512G --qos pq_madlab --account iacc_madlab -N 1 -n 1'), 'overwrite': True})

