#!/usr/bin/env python

import os
import nipype.interfaces.io as nio
import nipype.interfaces.fsl as fsl
from nipype.interfaces.utility import Merge
import nipype.pipeline.engine as pe
import nipype.interfaces.freesurfer as fs


def pickfirst(files):
    if isinstance(files, list):
        return files[0]
    return files

def get_aparc_aseg(files):
    for name in files:
        if 'aparc+aseg' in name:
            return name
    raise ValueError('aparc+aseg.mgz not found')


def dwi_preproc_wf(subject_id, gpu_card, sink_directory, name='dwi_preproc'):
    dwi_preproc_wf = pe.Workflow(name='dwi_preproc_wf')
    
    os.system('export SUBJECTS_DIR=/home/data/madlab/Mattfeld_vCAT/derivatives/freesurfer')
    os.environ['CUDA_VISIBLE_DEVICES'] = gpu_card
    subjects_dir = '/home/data/madlab/Mattfeld_vCAT/derivatives/freesurfer'
    data_dir = '/home/data/madlab/Mattfeld_vCAT'
    sink_dir = '/home/data/madlab/Mattfeld_vCAT/derivatives/preproc/dmri'
    
    info = dict(dwi= [[subject_id]],
                fmap= [[subject_id]],
                bvec= [[subject_id]],
                bval= [[subject_id]])

    
    #create a datasource node to get the dwi data and bvec/bval files
    datasource = pe.Node(nio.DataGrabber(infields=['subject_id'], 
                                         outfields=list(info.keys())), 
                         name='datasource')
    datasource.inputs.template = '*'
    datasource.inputs.subject_id = subject_id
    datasource.inputs.base_directory = os.path.abspath(data_dir)
    datasource.inputs.field_template = dict(dwi= 'dset/sub-%s/dwi/*_dwi.nii.gz',
                                            fmap= 'dset/sub-%s/fmap/*_dwi_dir-AP.nii.gz',
                                            bvec= 'dset/sub-%s/dwi/*_dwi.bvec',
                                            bval= 'dset/sub-%s/dwi/*_dwi.bval')
    datasource.inputs.template_args = info
    datasource.inputs.sort_filelist = True
    datasource.inputs.ignore_exception = False
    datasource.inputs.raise_on_empty = True

    
    #register source file to fs space and create brain mask in source space
    fssource = pe.Node(nio.FreeSurferSource(), 
                       name='fssource')
    fssource.inputs.subject_id = "sub-" + subject_id
    fssource.inputs.subjects_dir = subjects_dir

    
    # THE FOLLOWING CODE CALCULATES REGISTRATION BETWEEN FS STRUCTURAL AND DWI SPACE USING BBREGISTER ALGORITHM

    #extract first b0 file from AP dwi data to register between dwi and structural space
    extractref = pe.Node(fsl.ExtractROI(t_size=1), 
                         name="extractref")
    extractref.inputs.t_min = 0
    dwi_preproc_wf.connect(datasource, 'dwi', extractref, 'in_file')

    
    #calculate transformation matrix from EPI to FreeSurfer space using BBReg 
    fs_register = pe.MapNode(fs.BBRegister(init='fsl'), 
                             iterfield=['source_file'], 
                             name='fs_register')
    fs_register.inputs.contrast_type = 't2'
    fs_register.inputs.out_fsl_file = True
    fs_register.inputs.subject_id = "sub-" + subject_id
    fs_register.inputs.subjects_dir = subjects_dir
    dwi_preproc_wf.connect(extractref, 'roi_file', fs_register, 'source_file')

    
    # THE FOLLOWING CODE BINARIZES APARC+ASEG FILE AND THEN TRANSFORMS TO DWI SPACE
    # USING TRANSFORMATIONS CALCULATED BY BBREG TO MOVE BINARY MASK TO DWI SPACE

    #extract aparc+aseg brain mask and binarize
    fs_threshold = pe.Node(fs.Binarize(min=0.5, out_type='nii'), 
                           name='fs_threshold')
    dwi_preproc_wf.connect(fssource, ('aparc_aseg', get_aparc_aseg), fs_threshold, 'in_file')

    
    #transform binarized aparc+aseg file to 1st volume of 1st run space
    fs_voltransform = pe.Node(fs.ApplyVolTransform(inverse=True), 
                              name='fs_transform')
    fs_voltransform.inputs.subjects_dir = subjects_dir
    dwi_preproc_wf.connect(extractref, 'roi_file', fs_voltransform, 'source_file')
    dwi_preproc_wf.connect(fs_register, ('out_reg_file', pickfirst), fs_voltransform, 'reg_file')
    dwi_preproc_wf.connect(fs_threshold, 'binary_file', fs_voltransform, 'target_file')

    
    #dilate binarized mask by 1 voxel that is now in dwi space
    fs_threshold2 = pe.Node(fs.Binarize(min=0.5, out_type='nii'), 
                            name='fs_threshold2')
    fs_threshold2.inputs.dilate = 1
    dwi_preproc_wf.connect(fs_voltransform, 'transformed_file', fs_threshold2, 'in_file')

    
    #THE FOLLOWING CODE EXTRACTS B0 FILES FOR TWO DIFFERENT PHASE ENCODE DIRECTIONS AND THEN COMBINES FOR INPUT TO TOPUP
    
    #extract first volume of first run as reference
    extract_b0 = pe.MapNode(fsl.ExtractROI(t_size=1), 
                            iterfield=['t_min'], 
                            name="extract_b0")
    extract_b0.inputs.t_min = [0, 1, 34, 48, 68, 88, 102]
    dwi_preproc_wf.connect(datasource, 'dwi', extract_b0, 'in_file')
    
    
    merge_with_fmap = pe.Node(Merge(2), 
                              name='merge_with_fmap')
    dwi_preproc_wf.connect(datasource, 'fmap', merge_with_fmap, 'in1')
    dwi_preproc_wf.connect(extract_b0, 'roi_file', merge_with_fmap, 'in2')
    
    
    merge_with_dwi = pe.Node(Merge(2), 
                             name='merge_with_dwi')
    dwi_preproc_wf.connect(datasource, 'fmap', merge_with_dwi, 'in1')
    dwi_preproc_wf.connect(datasource, 'dwi', merge_with_dwi, 'in2')

    
    merge_b0 = pe.Node(fsl.Merge(), 
                       name='merge_b0')
    merge_b0.inputs.dimension = 't'
    dwi_preproc_wf.connect(merge_with_fmap, 'out', merge_b0, 'in_files')

    
    merge_dwi = pe.Node(fsl.Merge(), 
                        name='merge_dwi')
    merge_dwi.inputs.dimension = 't'
    dwi_preproc_wf.connect(merge_with_dwi, 'out', merge_dwi, 'in_files')
    
    
    #run topup on merged b0s
    topup = pe.Node(fsl.TOPUP(), 
                    name='topup')
    topup.inputs.encoding_file = data_dir + '/code/dwi/acqparams.txt' #use emu version
    topup.inputs.config = data_dir + '/code/dwi/b02b0.cnf'
    dwi_preproc_wf.connect(merge_b0, 'merged_file', topup, 'in_file')

    
    #run Eddy correct on combined data using input from topup
    eddy_corr = pe.Node(fsl.Eddy(), 
                        name='eddy_corr')
    eddy_corr.inputs.in_acqp = data_dir + '/code/dwi/acqparams.txt'
    eddy_corr.inputs.in_index = data_dir + '/code/dwi/dwi.index' 
    eddy_corr.inputs.use_cuda = True
    eddy_corr.inputs.is_shelled = True
    eddy_corr.inputs.args = '--slspec=/home/data/madlab/Mattfeld_vCAT/code/dwi/sltiming.txt'
    eddy_corr.plugin_args = {'sbatch_args': ('-p GPU_madlab --qos pq_madlab -n 1 --cpus-per-task=8')}
    dwi_preproc_wf.connect(datasource, 'dwi', eddy_corr, 'in_file')
    dwi_preproc_wf.connect(datasource, 'bval', eddy_corr, 'in_bval') 
    dwi_preproc_wf.connect(datasource, 'bvec', eddy_corr, 'in_bvec')
    dwi_preproc_wf.connect(fs_threshold2, 'binary_file', eddy_corr, 'in_mask')
    dwi_preproc_wf.connect(topup, 'out_fieldcoef', eddy_corr, 'in_topup_fieldcoef')
    dwi_preproc_wf.connect(topup, 'out_movpar', eddy_corr, 'in_topup_movpar')

    
    #datasink node to save output
    datasink = pe.Node(nio.DataSink(), 
                       name='datasink')
    datasink.inputs.base_directory = os.path.abspath(sink_dir)
    datasink.inputs.container = "sub-" + subject_id
    dwi_preproc_wf.connect(eddy_corr, 'out_corrected', datasink, 'corrected_dwi.@corrected')
    dwi_preproc_wf.connect(eddy_corr, 'out_rotated_bvecs', datasink, 'corrected_dwi.@bvecs')

    return dwi_preproc_wf

    
#############################
### Creates full workflow ###
#############################

def create_dwi_preproc_workflow(args, name='dwi_preproc'):
    kwargs = dict(subject_id=args.subject_id, gpu_card=args.gpu_card, sink_directory=os.path.abspath(args.out_dir), name=name)
    dwi_preproc_workflow = dwi_preproc_wf(**kwargs)
    return dwi_preproc_workflow

if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("-s", "--subject_id", dest="subject_id", help="Current subject id", required=True)
    parser.add_argument("-g", "--gpu_card", dest="gpu_card", default=0, help="Which gpu card 0 or 1")
    parser.add_argument("-o", "--output_dir", dest="out_dir", help="Output directory base")
    parser.add_argument("-w", "--work_dir", dest="work_dir", help="Working directory base")
    args = parser.parse_args()
    wf = create_dwi_preproc_workflow(args)

    if args.work_dir:
        work_dir = os.path.abspath(args.work_dir)
    else:
        work_dir = os.getcwd()

    wf.config['execution']['crashdump_dir'] = '/scratch/madlab/crash/mandy/Mattfeld_vCAT/dmri'
    wf.base_dir = work_dir + '/' + args.subject_id
    wf.run(plugin ='SLURM',plugin_args={'sbatch_args':('-p centos7 --qos pq_madlab --account iacc_madlab -N 1 -n 1'), 'overwrite': True})
