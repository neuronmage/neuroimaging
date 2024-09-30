#!/usr/bin/env python

"""
=====================================================================================
wmaze_fMRI: Mandy Thesis -- Fixed before Conditional -- Model 3 Version 1.3.2-B_split
=====================================================================================
Second level (Fixed Effects) workflow for UM GE 750 wmaze task data

- WMAZE Model 3 Version 1.3.2-B_split 
  - Use FSL ROI to recreate EPI data, removing last 3 volumes
  - Removed last 3 trials before EV creation
  - EV directory (Model 3) --- /home/data/madlab/data/mri/wmaze/scanner_behav/WMAZE_001/MRthesis/model3_1-3-2-B_split
  - Splits set into individual runs to compare early and late learning on B trials only


- python wmaze_lvl2.py -s WMAZE_001
                       -o /home/data/madlab/data/mri/wmaze/scndlvl/wmaze_MRthesis/fixed_before_conditional/model3_1-3-2-B_split
                       -w /home/data/madlab/scripts/wmaze/anal_MR_thesis/fixed_before_conditional/model3/model3_1-3-2-B_split/status/lvl2 
"""
import os
from nipype.pipeline.engine import Workflow, Node, MapNode
from nipype.interfaces.utility import IdentityInterface, Function
from nipype.utils.misc import getsource
from nipype.interfaces.fsl.model import L2Model, FLAMEO
from nipype.interfaces.io import DataGrabber,  DataSink
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

def secondlevel_wf(subject_id,
                   sink_directory,
                   name = 'model3_1-3-2-B_split_scndlvl_wf'):
    
    scndlvl_wf = Workflow(name = 'scndlvl_wf')

    # The contrast names
    contrasts = ['corr_B', 'incorr_B', 'all_remaining',  
                 'all_corr_minus_incorr', 'all_incorr_minus_corr']


    # Node: contrast_iterable
    contrast_iterable = Node(IdentityInterface(fields = ['contrast'], 
                                               mandatory_inputs = True), 
                             name = 'contrast_iterable')
    contrast_iterable.iterables = ('contrast', contrasts)

    info = dict(copes1st = [['subject_id', [0,2,4], 'contrast']],
                copes2nd = [['subject_id', [1,3,5], 'contrast']],
                varcopes1st = [['subject_id', [0,2,4], 'contrast']],
                varcopes2nd = [['subject_id', [1,3,5], 'contrast']],
                mask_file = [['subject_id', 'aparc+aseg_thresh']],
                dof_files1st = [['subject_id', [0,2,4], 'dof']],
                dof_files2nd = [['subject_id', [1,3,5], 'dof']])


    # Create a datasource node to get the task_mri and motion-noise files
    datasource = Node(DataGrabber(infields = ['subject_id','contrast'], 
                                  outfields = info.keys()), 
                                  name = 'datasource')
    datasource.inputs.template = '*'
    datasource.inputs.subject_id = subject_id
    datasource.inputs.base_directory = os.path.abspath('/home/data/madlab/data/mri/wmaze/')
    datasource.inputs.field_template = dict(copes1st = 'frstlvl/wmaze_MRthesis/fixed_before_conditional/model3_1-3-2-B/%s/modelfit/contrasts/_estimate_model%d/cope*_%s.nii.gz',
                                            copes2nd = 'frstlvl/wmaze_MRthesis/fixed_before_conditional/model3_1-3-2-B/%s/modelfit/contrasts/_estimate_model%d/cope*_%s.nii.gz',
                                            varcopes1st = 'frstlvl/wmaze_MRthesis/fixed_before_conditional/model3_1-3-2-B/%s/modelfit/contrasts/_estimate_model%d/varcope*_%s.nii.gz',
                                            varcopes2nd = 'frstlvl/wmaze_MRthesis/fixed_before_conditional/model3_1-3-2-B/%s/modelfit/contrasts/_estimate_model%d/varcope*_%s.nii.gz',
                                            mask_file = 'preproc/%s/ref/_fs_threshold20/%s*_thresh.nii',
                                            dof_files1st = 'frstlvl/wmaze_MRthesis/fixed_before_conditional/model3_1-3-2-B/%s/modelfit/dofs/_estimate_model%d/%s',
                                            dof_files2nd = 'frstlvl/wmaze_MRthesis/fixed_before_conditional/model3_1-3-2-B/%s/modelfit/dofs/_estimate_model%d/%s')
    datasource.inputs.template_args = info
    datasource.inputs.sort_filelist = True
    datasource.inputs.ignore_exception = False
    datasource.inputs.raise_on_empty = True
    scndlvl_wf.connect(contrast_iterable, 'contrast', datasource, 'contrast')


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
                            name = 'model3_1-3-2-B_split_scndlvl_wf'):

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

    wf.config['execution']['crashdump_dir'] = '/scratch/madlab/crash/wmaze_MRthesis/model3_1-3-2-B_split/lvl2'
    wf.base_dir = work_dir + '/' + args.subject_id
    wf.run(plugin = 'LSF', plugin_args = {'bsub_args': '-q PQ_madlab'})


