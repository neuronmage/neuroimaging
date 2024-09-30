#!/usr/bin/env python

import os
from nipype.pipeline.engine import Workflow
from nipype.pipeline.engine import Node
from nipype.interfaces.freesurfer.preprocess import MRIConvert
from nipype.interfaces.freesurfer.model import Binarize
from nipype.interfaces.freesurfer.utils import ApplyMask
from nipype.interfaces.utility import IdentityInterface
from nipype.interfaces.io import FreeSurferSource

def create_freesurfer_skullstrip_workflow(name = 'freesurfer_skullstrip'):

    # Workflow
    skull_strip = Workflow(name = name)

    
    #Set up node to define all inputs and outputs required for the freesurfer skull stripping workflow
    inputspec = Node(IdentityInterface(fields = ['subject_id','subjects_dir']),
                     name = 'inputspec')
    outputspec = Node(IdentityInterface(fields = ['mask_file','skullstripped_file']),
                      name = 'outputspec')


    #simple function to pick second
    picksecond = lambda x: x[1]
  
  

    # Node to source the FS files
    skull_strip_fssource = Node(FreeSurferSource(), 
                                name = 'skull_strip_fssource')
    skull_strip_fssource.inputs.hemi = 'both'
    skull_strip_fssource.inputs.ignore_exception = False
    skull_strip.connect(inputspec, 'subject_id', skull_strip_fssource, 'subject_id')
    skull_strip.connect(inputspec, 'subjects_dir', skull_strip_fssource, 'subjects_dir')



    # Node to define the segmentation and original brain to be stripped
    skull_strip_inputspec = Node(IdentityInterface(fields = ['brain', 'segmentation'], mandatory_inputs = True), 
                                 name = 'skull_strip_inputspec')
    skull_strip.connect(skull_strip_fssource, ('aparc_aseg', picksecond), skull_strip_inputspec, 'segmentation')
    skull_strip.connect(skull_strip_fssource, 'orig', skull_strip_inputspec, 'brain')



    # Node to convert the segmentation file to NIFTI 
    skull_strip_aparcaseg_2nii = Node(MRIConvert(), 
                                      name = "skull_strip_aparcaseg_2nii")
    skull_strip_aparcaseg_2nii.inputs.ignore_exception = False
    skull_strip_aparcaseg_2nii.inputs.out_type = 'nii'
    skull_strip.connect(skull_strip_inputspec, 'segmentation', skull_strip_aparcaseg_2nii, 'in_file')
    skull_strip.connect(inputspec, 'subjects_dir', skull_strip_aparcaseg_2nii, 'subjects_dir')



    # Node to create a binarized mask of the segmentation file -- 1 where brain, 0 where non-brain
    skull_strip_create_mask = Node(Binarize(), 
                                   name = "skull_strip_create_mask")
    skull_strip_create_mask.inputs.dilate = 1
    skull_strip_create_mask.inputs.ignore_exception = False
    skull_strip_create_mask.inputs.min = 1.0
    skull_strip_create_mask.inputs.out_type = 'nii'
    skull_strip.connect(skull_strip_aparcaseg_2nii, 'out_file', skull_strip_create_mask, 'in_file')
    skull_strip.connect(inputspec, 'subjects_dir', skull_strip_create_mask, 'subjects_dir')
    # Take the output from the binarize node as the mask output file for output node
    skull_strip.connect(skull_strip_create_mask, 'binary_file', outputspec, 'mask_file')



    # Node to convert non-stripped brain to NIFTI format (from .mgz)
    skull_strip_brain_2nii = Node(MRIConvert(), 
                                  name = "skull_strip_brain_2nii")
    skull_strip_brain_2nii.inputs.ignore_exception = False
    skull_strip_brain_2nii.inputs.out_type = 'nii'
    skull_strip.connect(skull_strip_inputspec, 'brain', skull_strip_brain_2nii, 'in_file')
    skull_strip.connect(inputspec, 'subjects_dir', skull_strip_brain_2nii, 'subjects_dir')



    # Node to apply the NIFTI mask to the NIFT brain -- removes non-brain (0) values
    skull_strip_apply_mask = Node(ApplyMask(), 
                                  name = "skull_strip_apply_mask")
    skull_strip_apply_mask.inputs.ignore_exception = False
    skull_strip.connect(inputspec, 'subjects_dir', skull_strip_apply_mask, 'subjects_dir')
    skull_strip.connect(skull_strip_create_mask, 'binary_file', skull_strip_apply_mask, 'mask_file')
    skull_strip.connect(skull_strip_brain_2nii, 'out_file', skull_strip_apply_mask, 'in_file')
    # Take the output from the applymask node as the skullstripped output file for the output node
    skull_strip.connect(skull_strip_apply_mask, 'out_file', outputspec, 'skullstripped_file')

    return skull_strip

