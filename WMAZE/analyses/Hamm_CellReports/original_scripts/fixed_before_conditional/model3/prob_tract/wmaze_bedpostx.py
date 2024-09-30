#!/usr/bin/env python

"""
=====================================================================================
wmaze_probtrak: Mandy Thesis -- Fixed before Conditional -- Model 3 Version 1.3.2
=====================================================================================
Workflow for WMAZE diffusion weighted data

- WMAZE Model 3 Version 1.3.2
1 - Traditional general linear model
3 - Removed last 3 trials before EV creation
2 - Use FSL ROI to recreate EPI data, removing last 3 volumes
    - EV directory (Model 3) --- /home/data/madlab/data/mri/wmaze/scanner_behav/WMAZE_001/MRthesis/model3_1-3-2


Workflow components:

-Freesurfer
-FSL
-bedpostx5 GPU version

"""

import nipype.interfaces.io as nio
import nipype.interfaces.fsl as fsl
import nipype.interfaces.utility as util
import nipype.pipeline.engine as pe
import nipype.interfaces.freesurfer as fs
import os

def pickfirst(in_file):
    if isinstance(in_file, list):
        return in_file[0]
    else:
        return in_file

def pickmiddle(func):
    """Return the middle volume index."""
    from nibabel import load
    return [(load(f).get_shape()[3]/2)-1 for f in func]

def pickvol(filenames, fileidx, which):
    from nibabel import load
    import numpy as np
    if which.lower() == 'first':
        idx = 0
    elif which.lower() == 'middle':
        idx = int(np.ceil(load(filenames[fileidx]).get_shape()[3]/2))
    else:
        raise Exception('unkown value for volume selection : %s'%which)
    return idx

def get_aparc_aseg(files):
    for name in files:
        if 'aparc+aseg' in name:
            return name
    raise ValueError('aparc+aseg.mgz not found')


def bpX_wf(subject_id,
           gpu_card,
           sink_directory,
           name='hcp_bpX'):
    
    subjects_dir = '/home/data/madlab/surfaces/wmaze'

    os.environ['CUDA_VISIBLE_DEVICES'] = gpu_card

    wmaze_dwi_wf = pe.Workflow(name=name)

    info = dict(dwi = [['subject_id', 'subject_id']],
                bvecs = [['subject_id']],
                bvals = [['subject_id']])
   
    # Create a datasource node to get the dwi, bvecs, and bvals
    datasource = pe.Node(interface = nio.DataGrabber(infields = ['subject_id'],
                                                     outfields = info.keys()),
                                                     name = 'datasource')
    datasource.inputs.subject_id = subject_id
    datasource.inputs.template = '%s/%s'
    datasource.inputs.base_directory = os.path.abspath('/home/data/madlab/data/mri/wmaze')
    datasource.inputs.field_template = dict(dwi = '%s/dmri/dtiprep_nrrd/%s_dwi_QCed.nii.gz',
                                            bvecs = '%s/dmri/bvec_QCed',
                                            bvals = '%s/dmri/bval_QCed')
    datasource.inputs.template_args = info
    datasource.inputs.sort_filelist = True

    # Extract the first volume of the first run as the reference 
    extractref = pe.Node(fsl.ExtractROI(t_size = 1),
                         iterfield = ['in_file'],
                         name = "extractref")
    extractref.inputs.t_min = 0
    wmaze_dwi_wf.connect(datasource, ('dwi', pickfirst), extractref, 'in_file')

    # Register a source file to fs space and create a brain mask in source space
    fssource = pe.Node(nio.FreeSurferSource(),
                       name = 'fssource')
    fssource.inputs.subject_id = subject_id
    fssource.inputs.subjects_dir = subjects_dir

    # Extract aparc+aseg brain mask and binarize
    fs_threshold = pe.Node(fs.Binarize(min = 0.5, out_type = 'nii'),
                           name = 'fs_threshold')
    wmaze_dwi_wf.connect(fssource, ('aparc_aseg', get_aparc_aseg), fs_threshold, 'in_file')

    # Calculate the transformation matrix from EPI space to FreeSurfer space
    # using the BBRegister command
    fs_register = pe.MapNode(fs.BBRegister(init = 'fsl'),
                             iterfield = ['source_file'],
                             name = 'fs_register')
    fs_register.inputs.contrast_type = 't2'
    fs_register.inputs.out_fsl_file = True
    fs_register.inputs.subject_id = subject_id
    fs_register.inputs.subjects_dir = subjects_dir
    wmaze_dwi_wf.connect(extractref, 'roi_file', fs_register, 'source_file')
    
    # Transform the binarized aparc+aseg file to the 1st volume of 1st run space
    fs_voltransform = pe.MapNode(fs.ApplyVolTransform(inverse = True),
                                 iterfield = ['source_file', 'reg_file'],
                                 name = 'fs_transform')
    fs_voltransform.inputs.subjects_dir = subjects_dir
    wmaze_dwi_wf.connect(extractref, 'roi_file', fs_voltransform, 'source_file')
    wmaze_dwi_wf.connect(fs_register, 'out_reg_file', fs_voltransform, 'reg_file')
    wmaze_dwi_wf.connect(fs_threshold, 'binary_file', fs_voltransform, 'target_file')

    # Dilate the binarized mask by 1 voxel that is now in the DWI space
    fs_threshold2 = pe.MapNode(fs.Binarize(min = 0.5, out_type = 'nii'),
                               iterfield = ['in_file'],
                               name = 'fs_threshold2')
    fs_threshold2.inputs.dilate = 1
    wmaze_dwi_wf.connect(fs_voltransform, 'transformed_file', fs_threshold2, 'in_file')

    # Create bedpostx node
    bedp = pe.Node(interface = fsl.BEDPOSTX5(), 
                   name = 'bedp')
    bedp.inputs.model = 2
    bedp.inputs.cnlinear = True
    bedp.inputs.gradnonlin = False
    bedp.inputs.rician = False
    bedp.inputs.use_gpu = True
    bedp.inputs.n_fibres = 2
    bedp.plugin_args = {'bsub_args': ('-q PQ_madlab_gpu -n 8')}
    wmaze_dwi_wf.connect(datasource, 'dwi', bedp, 'dwi')
    wmaze_dwi_wf.connect(datasource, 'bvals', bedp, 'bvals')
    wmaze_dwi_wf.connect(datasource, 'bvecs', bedp, 'bvecs')
    wmaze_dwi_wf.connect(fs_threshold2, ('binary_file', pickfirst), bedp, 'mask')

    #ADD DATASINK NODE TO CAPTURE OUPUT FROM BEDPOSTX
    datasink = pe.Node(interface = nio.DataSink(),
                       name = 'datasink')
    datasink.inputs.base_directory = os.path.abspath(sink_directory)
    datasink.inputs.container = subject_id
    wmaze_dwi_wf.connect(bedp, 'merged_thsamples', datasink, 'wmazebpX.thsamples')
    wmaze_dwi_wf.connect(bedp, 'merged_phsamples', datasink, 'wmazebpX.phsamples')
    wmaze_dwi_wf.connect(bedp, 'merged_fsamples', datasink, 'wmazebpX.fsamples')
    wmaze_dwi_wf.connect(extractref, 'roi_file', datasink, 'wmazebpX.b0_file')
    wmaze_dwi_wf.connect(fs_threshold2, 'binary_file', datasink, 'wmazebpX.mask')

    return wmaze_dwi_wf

"""
Creates the full workflow
"""

def create_bpX_workflow(args, name = 'wmaze_bpX'):
    
    kwargs = dict(subject_id = args.subject_id,
                  gpu_card = args.gpu_card,
                  sink_directory = os.path.abspath(args.out_dir),
                  name = name)
    bpX_workflow = bpX_wf(**kwargs)
    return bpX_workflow

if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("-s", "--subject_id", dest = "subject_id",
                        help = "Current subject id", required = True)
    parser.add_argument("-g", "--gpu_card", dest = "gpu_card",
                        default = 0, help = "Which gpu card 0 or 1")
    parser.add_argument("-o", "--output_dir", dest = "out_dir",
                        help = "Output directory base")
    parser.add_argument("-w", "--work_dir", dest = "work_dir",
                        help = "Working directory base")
    args = parser.parse_args()

    wf = create_bpX_workflow(args)

    if args.work_dir:
        work_dir = os.path.abspath(args.work_dir)
    else:
        work_dir = os.getcwd()
    
    wf.config['execution']['crashdump_dir'] = '/scratch/madlab/wmaze/wmaze_bpX'
    wf.base_dir = work_dir + '/' + args.subject_id
    wf.run(plugin='LSF', plugin_args={'bsub_args': '-q PQ_madlab'})

