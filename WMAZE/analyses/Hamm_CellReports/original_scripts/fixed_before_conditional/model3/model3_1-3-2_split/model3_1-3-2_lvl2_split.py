#!/usr/bin/env python

"""
=====================================================================================
wmaze_fMRI: Mandy Thesis -- Fixed before Conditional -- Model 3 Version 1.3.2_split
=====================================================================================
Second level (Fixed Effects) workflow for UM GE 750 wmaze task data

- WMAZE Model 3 Version 1.3.2_split 
  - Use FSL ROI to recreate EPI data, removing last 3 volumes
  - Removed last 3 trials before EV creation
  - EV directory (Model 3) --- /home/data/madlab/data/mri/wmaze/scanner_behav/WMAZE_001/MRthesis/model3_1-3-2_split
  - Splits set into individual runs to compare early and late learning for fixed before conditional trials


- python wmaze_lvl2.py -s WMAZE_001
                       -o /home/data/madlab/data/mri/wmaze/scndlvl/wmaze_MRthesis/fixed_before_conditional/model3_1-3-2_split
                       -w /home/data/madlab/scripts/wmaze/anal_MR_thesis/fixed_before_conditional/model3_1-3-2_split/status/lvl2 
"""
import os
from glob import glob
from nipype.pipeline.engine import Workflow, Node, MapNode
from nipype.interfaces.utility import IdentityInterface, Function
from nipype.utils.misc import getsource
from nipype.interfaces.fsl.model import L2Model, FLAMEO
from nipype.interfaces.io import DataGrabber, DataSink
from nipype.interfaces.fsl.utils import Merge

# Functions
def num_copes(files):
    if type(files[0]) is list:
        return len(files[0])
    else:
        return len(files)

def num_copes_range(files):
    if type(files[0]) is list:
        return range(0,len(files[0]))
    else:
        return range(0,len(files))

def doublelist(x):
    if isinstance(x[0],list):
        return x
    else:
        return [x]


def get_contrasts(data_inputs):
    import os
    infiles = [os.path.split(d[0])[1] for d in data_inputs]
    contrasts = [inf[7:].split('.nii')[0] for inf in infiles]
    print contrasts
    return contrasts


def get_subs(subject_id, cons):
    subs = []
    for i, con in enumerate(cons):
        subs.append(('_contrast_{0}/_flameo_fe1st0/cope1'.format(con), 'cope_{0}'.format(con)))
        subs.append(('_contrast_{0}/_flameo_fe1st0/varcope1'.format(con), 'varcope_{0}'.format(con)))
        subs.append(('_contrast_{0}/_flameo_fe1st0/tstat1'.format(con), 'tstat_{0}'.format(con)))
        subs.append(('_contrast_{0}/_flameo_fe1st0/zstat1'.format(con), 'zstat_{0}'.format(con)))
        subs.append(('_contrast_{0}/_flameo_fe1st0/res4d'.format(con), 'res4d_{0}'.format(con)))
        subs.append(('_contrast_{0}/_flameo_fe2nd0/cope1'.format(con), 'cope_{0}'.format(con)))
        subs.append(('_contrast_{0}/_flameo_fe2nd0/varcope1'.format(con), 'varcope_{0}'.format(con)))
        subs.append(('_contrast_{0}/_flameo_fe2nd0/tstat1'.format(con), 'tstat_{0}'.format(con)))
        subs.append(('_contrast_{0}/_flameo_fe2nd0/zstat1'.format(con), 'zstat_{0}'.format(con)))
        subs.append(('_contrast_{0}/_flameo_fe2nd0/res4d'.format(con), 'res4d_{0}'.format(con)))
    return subs

'''
def get_dofvolumes(dof_files, cope_files):
    import os
    import nibabel as nb
    import numpy as np
    filenames = []
    for curr_cope_file_num, curr_cope_file in enumerate(cope_files):
        curr_cope_file_name = os.path.basename(curr_cope_file)
        initial_name, _ = curr_cope_file_name.split('_merged')
        curr_cope_name = initial_name[7:]
        img = nb.load(curr_cope_file)
        curr_out_data = np.zeros(img.get_shape())
        if curr_out_data.shape[-1] != 56:
            for i in range(curr_out_data.shape[-1]):
                curr_out_data[:, :, :, i] = np.loadtxt(dof_files[curr_cope_file_num][i])
        else:
            curr_out_data = np.expand_dims(curr_out_data, 3)
            curr_out_data[:, :, :, 0] = np.loadtxt(dof_files[curr_cope_file_num][-1])
        curr_filename = os.path.join(os.getcwd(), 'dof_file_{0}.nii.gz'.format(curr_cope_name))
        curr_newimg = nb.Nifti1Image(curr_out_data, None, img.get_header())
        curr_newimg.to_filename(curr_filename)
        filenames.append(curr_filename)
    return filenames
'''

def get_dofvolumes(dof_files, cope_files):
    import os
    import nibabel as nb
    import numpy as np
    img = nb.load(cope_files[0])
    out_data = np.zeros(img.get_shape())
    for i in range(out_data.shape[-1]):
        dof = np.loadtxt(dof_files[i])
        out_data[:, :, :, i] = dof
    filename = os.path.join(os.getcwd(), 'dof_file.nii.gz')
    newimg = nb.Nifti1Image(out_data, None, img.get_header())
    newimg.to_filename(filename)
    return filename


def get_first_runs(subject_id, contrast): 
    import os
    from glob import glob
    proj_dir = os.path.abspath('/home/data/madlab/data/mri/wmaze/')  
    first_runs = []
    contrast_files = glob('{0}/frstlvl/wmaze_MRthesis/fix*/model3_1-3-2/{1}/mod*/cont*/_est*_model*/cope??_{2}.nii.gz'.format(proj_dir, 
                                                                                                                              subject_id, 
                                                                                                                              contrast))
    for curr_file in contrast_files: #for each of the contrast file lists
        if curr_file.split('/')[-2][-1] in ['0', '2', '4']: #if Run 1
            first_runs.append(curr_file.split('/')[-2][-1]) #append to first runs
    return first_runs    


def get_second_runs(subject_id, contrast): 
    import os
    from glob import glob 
    proj_dir = os.path.abspath('/home/data/madlab/data/mri/wmaze/') 
    second_runs = []
    contrast_files = glob('{0}/frstlvl/wmaze_MRthesis/fix*/model3_1-3-2/{1}/mod*/cont*/_est*_model*/cope??_{2}.nii.gz'.format(proj_dir, 
                                                                                                                              subject_id, 
                                                                                                                              contrast))
    for curr_file in contrast_files: #for each of the contrast file lists
        if curr_file.split('/')[-2][-1] in ['1', '3', '5']: #if Run 2
            second_runs.append(curr_file.split('/')[-2][-1]) #append to second runs
    return second_runs 


################
### Pipeline ###
################

def secondlevel_wf(subject_id,
                   sink_directory,
                   name = 'model3_1-3-2_split_scndlvl_wf'):
    
    scndlvl_wf = Workflow(name = 'scndlvl_wf')
    proj_dir = os.path.abspath('/home/data/madlab/data/mri/wmaze/')
   
    contrasts = ['all_before_B_corr', 'all_before_B_incorr', 'all_remaining',  
                 'all_corr_minus_all_incorr', 'all_incorr_minus_all_corr']

    contrasts = [f for f in contrasts if len(glob('{0}/frstlvl/wmaze_MRthesis/fix*/model3_1-3-2/{1}/mod*/cont*/_est*_model*/cope??_{2}.nii.gz'.format(proj_dir, subject_id, f))) >= 1]

    contrast_files = []
    #for each contrast, grab its copes from all runs (array of arrays -- 5)
    for curr_cont in contrasts:
        contrast_files.append(glob('{0}/frstlvl/wmaze_MRthesis/fix*/model3_1-3-2/{1}/mod*/cont*/_est*_model*/cope??_{2}.nii.gz'.format(proj_dir, 
                                                                                                                                        subject_id, 
                                                                                                                                        curr_cont)))
    '''
    first_runs = []
    second_runs = [] 
    for _ in contrasts:
        first_runs.append([])
        second_runs.append([])
    for i, curr_file_list in enumerate(contrast_files): #for each of the contrast file lists
        if not isinstance(curr_file_list, list): #if it is a single file, turn it into a list
            curr_file_list = [curr_file_list]
        for curr_file in curr_file_list: #for each file in the contrast list
            if curr_file.split('/')[-2][-1] in ['0', '2', '4']: #if Run1
                first_runs[i].append(curr_file.split('/')[-2][-1]) #append to first runs
            else:
                second_runs[i].append(curr_file.split('/')[-2][-1])#if Run2, append to second runs
    
   

    dof_runs = []
    for _ in contrasts:
        dof_runs.append([])
    for i, curr_file_list in enumerate(contrast_files):
        if not isinstance(curr_file_list, list):
            curr_file_list = [curr_file_list]
        for curr_file in curr_file_list:
            dof_runs[i].append(curr_file.split('/')[-2][-1])
    with open('/home/arenf001/stupidshit.txt', 'w') as f:
        f.write(str(dof_runs)) 
        f.write('\n')
        f.write(str(len(dof_runs)))
    '''

    contrast_iterable = Node(IdentityInterface(fields = ['contrast'], 
                                               mandatory_inputs = True), 
                             name = 'contrast_iterable')
    contrast_iterable.iterables = ('contrast', contrasts)

    # Get existing first runs for each contrast
    get_firstruns = Node(Function(input_names = ['subject_id', 'contrast'],
                                    output_names = ['first_runs'],
                                    function = get_first_runs),
                           name = 'get_firstruns')
    get_firstruns.inputs.ignore_exception = False
    get_firstruns.inputs.subject_id = subject_id
    scndlvl_wf.connect(contrast_iterable, 'contrast', get_firstruns, 'contrast')

    # Get existing second runs for each contrast
    get_secondruns = Node(Function(input_names = ['subject_id', 'contrast'],
                                    output_names = ['second_runs'],
                                    function = get_second_runs),
                           name = 'get_secondruns')
    get_secondruns.inputs.ignore_exception = False
    get_secondruns.inputs.subject_id = subject_id
    scndlvl_wf.connect(contrast_iterable, 'contrast', get_secondruns, 'contrast')


    ##ISSUE: 
    info = dict(copes1st = [['subject_id', 'first_runs', 'contrast']], 
                copes2nd = [['subject_id', 'second_runs', 'contrast']],
                varcopes1st = [['subject_id', 'first_runs', 'contrast']],
                varcopes2nd = [['subject_id', 'second_runs', 'contrast']],
                mask_file = [['subject_id', 'aparc+aseg_thresh']],
                dof_files1st = [['subject_id', 'first_runs', 'dof']],
                dof_files2nd = [['subject_id', 'second_runs', 'dof']])


    # Create a datasource node to get the task_mri and motion-noise files
    datasource = Node(DataGrabber(infields = ['subject_id','contrast', 
                                              'first_runs', 'second_runs'], 
                                  outfields = info.keys()), 
                                  name = 'datasource')
    datasource.inputs.template = '*'
    datasource.inputs.subject_id = subject_id
    datasource.inputs.base_directory = proj_dir
    datasource.inputs.field_template = dict(copes1st = 'frstlvl/wmaze_MRthesis/fixed_before_conditional/model3_1-3-2/%s/modelfit/contrasts/_estimate_model%s/cope*_%s.nii.gz',
                                            copes2nd = 'frstlvl/wmaze_MRthesis/fixed_before_conditional/model3_1-3-2/%s/modelfit/contrasts/_estimate_model%s/cope*_%s.nii.gz',
                                            varcopes1st = 'frstlvl/wmaze_MRthesis/fixed_before_conditional/model3_1-3-2/%s/modelfit/contrasts/_estimate_model%s/varcope*_%s.nii.gz',
                                            varcopes2nd = 'frstlvl/wmaze_MRthesis/fixed_before_conditional/model3_1-3-2/%s/modelfit/contrasts/_estimate_model%s/varcope*_%s.nii.gz',
                                            mask_file = 'preproc/%s/ref/_fs_threshold20/%s*_thresh.nii',
                                            dof_files1st = 'frstlvl/wmaze_MRthesis/fixed_before_conditional/model3_1-3-2/%s/modelfit/dofs/_estimate_model%s/%s',
                                            dof_files2nd = 'frstlvl/wmaze_MRthesis/fixed_before_conditional/model3_1-3-2/%s/modelfit/dofs/_estimate_model%s/%s')
    datasource.inputs.template_args = info
    datasource.inputs.sort_filelist = True
    datasource.inputs.ignore_exception = False
    datasource.inputs.raise_on_empty = True
    scndlvl_wf.connect(contrast_iterable, 'contrast', datasource, 'contrast')
    scndlvl_wf.connect(get_firstruns, 'first_runs', datasource, 'first_runs')
    scndlvl_wf.connect(get_secondruns, 'second_runs', datasource, 'second_runs')


    # Create an Inputspec node to deal with copes and varcopes doublelist issues
    fixedfx_inputspec = Node(IdentityInterface(fields = ['copes1st', 'copes2nd', 
                                                         'varcopes1st', 'varcopes2nd', 
                                                         'dof_files1st', 'dof_files2nd'],
                                               mandatory_inputs = True),
                             name = 'fixedfx_inputspec')
    scndlvl_wf.connect(datasource, ('copes1st',doublelist), fixedfx_inputspec, 'copes1st')
    scndlvl_wf.connect(datasource, ('copes2nd',doublelist), fixedfx_inputspec, 'copes2nd')
    scndlvl_wf.connect(datasource, ('varcopes1st',doublelist), fixedfx_inputspec, 'varcopes1st')
    scndlvl_wf.connect(datasource, ('varcopes2nd',doublelist), fixedfx_inputspec, 'varcopes2nd')
    scndlvl_wf.connect(datasource, 'dof_files1st', fixedfx_inputspec, 'dof_files1st')
    scndlvl_wf.connect(datasource, 'dof_files2nd', fixedfx_inputspec, 'dof_files2nd')
 

    # Create a Merge node to collect all of the RUN 1 COPES
    copemerge1st = MapNode(Merge(), iterfield = ['in_files'], 
                           name = 'copemerge1st')
    copemerge1st.inputs.dimension = 't'
    copemerge1st.inputs.environ = {'FSLOUTPUTTYPE': 'NIFTI_GZ'}
    copemerge1st.inputs.ignore_exception = False
    copemerge1st.inputs.output_type = 'NIFTI_GZ'
    copemerge1st.inputs.terminal_output = 'stream'
    scndlvl_wf.connect(fixedfx_inputspec, 'copes1st', copemerge1st, 'in_files')


    # Create a Merge node to collect all of the RUN 2 COPES
    copemerge2nd = MapNode(Merge(), iterfield = ['in_files'], 
                           name = 'copemerge2nd')
    copemerge2nd.inputs.dimension = 't'
    copemerge2nd.inputs.environ = {'FSLOUTPUTTYPE': 'NIFTI_GZ'}
    copemerge2nd.inputs.ignore_exception = False
    copemerge2nd.inputs.output_type = 'NIFTI_GZ'
    copemerge2nd.inputs.terminal_output = 'stream'
    scndlvl_wf.connect(fixedfx_inputspec, 'copes2nd', copemerge2nd, 'in_files')   


    # Create a Function node to generate a DOF volume for RUN 1
    gendofvolume1st = Node(Function(input_names = ['dof_files', 'cope_files'],
                                    output_names = ['dof_volume'],
                                    function = get_dofvolumes),
                           name = 'gendofvolume1st')
    gendofvolume1st.inputs.ignore_exception = False
    scndlvl_wf.connect(fixedfx_inputspec, 'dof_files1st', gendofvolume1st, 'dof_files')
    scndlvl_wf.connect(copemerge1st, 'merged_file', gendofvolume1st, 'cope_files')


    # Create a Function node to generate a DOF volume for RUN 2
    gendofvolume2nd = Node(Function(input_names = ['dof_files', 'cope_files'],
                                    output_names = ['dof_volume'],
                                    function = get_dofvolumes),
                           name = 'gendofvolume2nd')
    gendofvolume2nd.inputs.ignore_exception = False
    scndlvl_wf.connect(fixedfx_inputspec, 'dof_files2nd', gendofvolume2nd, 'dof_files')
    scndlvl_wf.connect(copemerge2nd, 'merged_file', gendofvolume2nd, 'cope_files')


    # Create a Merge node to collect all of the RUN 1 VARCOPES
    varcopemerge1st = MapNode(Merge(), iterfield=['in_files'], 
                              name = 'varcopemerge1st')
    varcopemerge1st.inputs.dimension = 't'
    varcopemerge1st.inputs.environ = {'FSLOUTPUTTYPE': 'NIFTI_GZ'}
    varcopemerge1st.inputs.ignore_exception = False
    varcopemerge1st.inputs.output_type = 'NIFTI_GZ'
    varcopemerge1st.inputs.terminal_output = 'stream'
    scndlvl_wf.connect(fixedfx_inputspec, 'varcopes1st', varcopemerge1st, 'in_files')


    # Create a Merge node to collect all of the RUN 2 VARCOPES
    varcopemerge2nd = MapNode(Merge(), iterfield = ['in_files'], 
                              name = 'varcopemerge2nd')
    varcopemerge2nd.inputs.dimension = 't'
    varcopemerge2nd.inputs.environ = {'FSLOUTPUTTYPE': 'NIFTI_GZ'}
    varcopemerge2nd.inputs.ignore_exception = False
    varcopemerge2nd.inputs.output_type = 'NIFTI_GZ'
    varcopemerge2nd.inputs.terminal_output = 'stream'
    scndlvl_wf.connect(fixedfx_inputspec, 'varcopes2nd', varcopemerge2nd, 'in_files')


    # Create a node to define the contrasts from the names of the copes
    getcontrasts = Node(Function(input_names = ['data_inputs'],
                                 output_names = ['contrasts'],
                                 function = get_contrasts),
                        name = 'getcontrasts')
    getcontrasts.inputs.ignore_exception = False
    scndlvl_wf.connect(datasource, ('copes1st',doublelist), getcontrasts, 'data_inputs')


    # Create a Function node to rename output files with something more meaningful
    getsubs = Node(Function(input_names = ['subject_id', 'cons'],
                            output_names = ['subs'],
                            function = get_subs),
                   name = 'getsubs')
    getsubs.inputs.ignore_exception = False
    getsubs.inputs.subject_id = subject_id
    scndlvl_wf.connect(getcontrasts, 'contrasts', getsubs, 'cons')


    # Create a l2model node for the Fixed Effects analysis (aka within subj across runs) for early and late learning
    l2model1st = Node(L2Model(), name = 'l2model1st')
    l2model1st.inputs.ignore_exception = False
    scndlvl_wf.connect(datasource, ('copes1st', num_copes), l2model1st, 'num_copes')


    # Create a l2model node for the Fixed Effects analysis (aka within subj across runs) for early and late learning
    l2model2nd = Node(L2Model(), name = 'l2model2nd')
    l2model2nd.inputs.ignore_exception = False
    scndlvl_wf.connect(datasource, ('copes2nd', num_copes), l2model2nd, 'num_copes')


    # Create a FLAMEO Node to run the fixed effects analysis for RUN 1
    flameo_fe1st = MapNode(FLAMEO(), iterfield = ['cope_file', 'var_cope_file'], 
                           name = 'flameo_fe1st')
    flameo_fe1st.inputs.environ = {'FSLOUTPUTTYPE': 'NIFTI_GZ'}
    flameo_fe1st.inputs.ignore_exception = False
    flameo_fe1st.inputs.log_dir = 'stats'
    flameo_fe1st.inputs.output_type = 'NIFTI_GZ'
    flameo_fe1st.inputs.run_mode = 'fe'
    flameo_fe1st.inputs.terminal_output = 'stream'
    scndlvl_wf.connect(varcopemerge1st, 'merged_file', flameo_fe1st, 'var_cope_file')
    scndlvl_wf.connect(l2model1st, 'design_mat', flameo_fe1st, 'design_file')
    scndlvl_wf.connect(l2model1st, 'design_con', flameo_fe1st, 't_con_file')
    scndlvl_wf.connect(l2model1st, 'design_grp', flameo_fe1st, 'cov_split_file')
    scndlvl_wf.connect(gendofvolume1st, 'dof_volume', flameo_fe1st, 'dof_var_cope_file')
    scndlvl_wf.connect(datasource, 'mask_file', flameo_fe1st, 'mask_file')
    scndlvl_wf.connect(copemerge1st, 'merged_file', flameo_fe1st, 'cope_file')


    # Create a FLAMEO Node to run the fixed effects analysis for RUN 2
    flameo_fe2nd = MapNode(FLAMEO(), iterfield = ['cope_file', 'var_cope_file'], 
                           name = 'flameo_fe2nd')
    flameo_fe2nd.inputs.environ = {'FSLOUTPUTTYPE': 'NIFTI_GZ'}
    flameo_fe2nd.inputs.ignore_exception = False
    flameo_fe2nd.inputs.log_dir = 'stats'
    flameo_fe2nd.inputs.output_type = 'NIFTI_GZ'
    flameo_fe2nd.inputs.run_mode = 'fe'
    flameo_fe2nd.inputs.terminal_output = 'stream'
    scndlvl_wf.connect(varcopemerge2nd, 'merged_file', flameo_fe2nd, 'var_cope_file')
    scndlvl_wf.connect(l2model2nd, 'design_mat', flameo_fe2nd, 'design_file')
    scndlvl_wf.connect(l2model2nd, 'design_con', flameo_fe2nd, 't_con_file')
    scndlvl_wf.connect(l2model2nd, 'design_grp', flameo_fe2nd, 'cov_split_file')
    scndlvl_wf.connect(gendofvolume2nd, 'dof_volume', flameo_fe2nd, 'dof_var_cope_file')
    scndlvl_wf.connect(datasource, 'mask_file', flameo_fe2nd, 'mask_file')
    scndlvl_wf.connect(copemerge2nd, 'merged_file', flameo_fe2nd, 'cope_file')


    # Create an outputspec node
    scndlvl_outputspec = Node(IdentityInterface(fields = ['res4d1st', 'copes1st', 'varcopes1st', 'zstats1st', 'tstats1st',
                                                          'res4d2nd', 'copes2nd', 'varcopes2nd', 'zstats2nd', 'tstats2nd'],
                                                mandatory_inputs = True),
                              name = 'scndlvl_outputspec')
    scndlvl_wf.connect(flameo_fe1st, 'res4d', scndlvl_outputspec, 'res4d1st')
    scndlvl_wf.connect(flameo_fe1st, 'copes', scndlvl_outputspec, 'copes1st')
    scndlvl_wf.connect(flameo_fe1st, 'var_copes', scndlvl_outputspec, 'varcopes1st')
    scndlvl_wf.connect(flameo_fe1st, 'zstats', scndlvl_outputspec, 'zstats1st')
    scndlvl_wf.connect(flameo_fe1st, 'tstats', scndlvl_outputspec, 'tstats1st')
    scndlvl_wf.connect(flameo_fe2nd, 'res4d', scndlvl_outputspec, 'res4d2nd')
    scndlvl_wf.connect(flameo_fe2nd, 'copes', scndlvl_outputspec, 'copes2nd')
    scndlvl_wf.connect(flameo_fe2nd, 'var_copes', scndlvl_outputspec, 'varcopes2nd')
    scndlvl_wf.connect(flameo_fe2nd, 'zstats', scndlvl_outputspec, 'zstats2nd')
    scndlvl_wf.connect(flameo_fe2nd, 'tstats', scndlvl_outputspec, 'tstats2nd')


    # Create a datasink node
    sinkd = Node(DataSink(), 
                 name = 'sinkd')
    sinkd.inputs.base_directory = sink_directory 
    sinkd.inputs.container = subject_id
    scndlvl_wf.connect(scndlvl_outputspec, 'copes1st', sinkd, 'fixedfx1sthalf.@copes')
    scndlvl_wf.connect(scndlvl_outputspec, 'varcopes1st', sinkd, 'fixedfx1sthalf.@varcopes')
    scndlvl_wf.connect(scndlvl_outputspec, 'tstats1st', sinkd, 'fixedfx1sthalf.@tstats')
    scndlvl_wf.connect(scndlvl_outputspec, 'zstats1st', sinkd, 'fixedfx1sthalf.@zstats')
    scndlvl_wf.connect(scndlvl_outputspec, 'res4d1st', sinkd, 'fixedfx1sthalf.@pvals')
    scndlvl_wf.connect(scndlvl_outputspec, 'copes2nd', sinkd, 'fixedfx2ndhalf.@copes')
    scndlvl_wf.connect(scndlvl_outputspec, 'varcopes2nd', sinkd, 'fixedfx2ndhalf.@varcopes')
    scndlvl_wf.connect(scndlvl_outputspec, 'tstats2nd', sinkd, 'fixedfx2ndhalf.@tstats')
    scndlvl_wf.connect(scndlvl_outputspec, 'zstats2nd', sinkd, 'fixedfx2ndhalf.@zstats')
    scndlvl_wf.connect(scndlvl_outputspec, 'res4d2nd', sinkd, 'fixedfx2ndhalf.@pvals')
    scndlvl_wf.connect(getsubs, 'subs', sinkd, 'substitutions')

    return scndlvl_wf

"""
Creates the full workflow
"""

def create_scndlvl_workflow(args, 
                            name = 'model3_1-3-2_split_scndlvl_wf'):

    kwargs = dict(subject_id = args.subject_id,
                  sink_directory = os.path.abspath(args.out_dir),
                  name = name)
    scndlvl_workflow = secondlevel_wf(**kwargs)
    return scndlvl_workflow


if __name__ == '__main__':
    from argparse import ArgumentParser
    parser = ArgumentParser(description = __doc__)
    parser.add_argument('-s', '--subject_id', dest = 'subject_id',
                        help = 'Current subject id', required = True)
    parser.add_argument('-o', '--output_dir', dest = 'out_dir',
                        help = 'Output directory base')
    parser.add_argument('-w', '--work_dir', dest = 'work_dir',
                        help = 'Working directory base')
    args = parser.parse_args()

    wf = create_scndlvl_workflow(args)

    if args.work_dir:
        work_dir = os.path.abspath(args.work_dir)
    else:
        work_dir = os.getcwd()

    wf.config['execution']['crashdump_dir'] = '/scratch/madlab/crash/wmaze_MRthesis/model3_1-3-2_split/lvl2'
    wf.base_dir = work_dir + '/' + args.subject_id
    wf.run(plugin = 'LSF', plugin_args = {'bsub_args': '-q PQ_madlab'})


