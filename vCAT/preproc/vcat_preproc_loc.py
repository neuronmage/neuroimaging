#!/usr/bin/env python
"""
================================================================
vCAT_fMRI: AFNI, FS, FSL, NiPy
================================================================
A preprocessing workflow for Siemens vCAT task data

"""

from warnings import warn
import os
import numpy as np
import scipy as sp
import nibabel as nb

from nipype.algorithms.confounds import TSNR
from nipype.algorithms.rapidart import ArtifactDetect
from nipype.interfaces.afni.preprocess import Despike
from nipype.interfaces.c3 import C3dAffineTool
from nipype.interfaces.freesurfer.model import Binarize
from nipype.interfaces.freesurfer.preprocess import ApplyVolTransform, BBRegister
from nipype.interfaces.fsl.utils import ImageMaths, ImageStats, ExtractROI, PlotMotionParams
from nipype.interfaces.fsl.preprocess import MCFLIRT, SUSAN
from nipype.interfaces.io import DataGrabber, DataSink, FreeSurferSource
from nipype.interfaces.nipy import SpaceTimeRealigner
from nipype.interfaces.utility.base import IdentityInterface, Merge, Rename 
from nipype.interfaces.utility.wrappers import Function
from nipype.pipeline.engine import Workflow, Node, MapNode 
from nipype.utils.filemanip import filename_to_list


imports = ['import os',
           'import nibabel as nb',
           'import numpy as np',
           'import scipy as sp',
           'from nipype.utils.filemanip import filename_to_list, list_to_filename, split_filename',
           'from scipy.special import legendre']


def pickfirst(func): #get the first item in a list
    if isinstance(func, list):
        return func[0]
    else:
        return func


def pickmiddle(func): #return middle volume index
    from nibabel import load
    return [(load(f).get_shape()[3]/2)-1 for f in func] #load func file, grab the last volume after halving 4D


def pickvol(filenames, fileidx, which): #get selected reference volume
    from nibabel import load
    import numpy as np
    if which.lower() == 'first': #if first volume indicated
        idx = 0 #get first volume
    elif which.lower() == 'middle': #if middle volume is indicated
        idx = int(np.ceil(load(filenames[fileidx]).get_shape()[3]/2))
    else:
        raise Exception('unknown value for volume selection : %s' % which) #if error in volume indication
    return idx


def motion_regressors(motion_params, order = 0, derivatives = 1): #motion regressor function to include 1st and 2nd derivatives later in the script 
    out_files = []
    for idx, filename in enumerate(filename_to_list(motion_params)):
        params = np.genfromtxt(filename)
        out_params = params
        for d in range(1, derivatives + 1):
            cparams = np.vstack((np.repeat(params[0, :][None, :], d, axis = 0),
                                 params))
            out_params = np.hstack((out_params, np.diff(cparams, d, axis = 0)))
        out_params2 = out_params
        for i in range(2, order + 1):
            out_params2 = np.hstack((out_params2, np.power(out_params, i)))
        filename = os.path.join(os.getcwd(), "motion_regressor%02d.txt" % idx)
        np.savetxt(filename, out_params2, fmt = "%.10f")
        out_files.append(filename)
    return out_files


"""
Parameters 
----------
motion_params: a text file containing motion parameters and its derivatives
comp_norm: a text file containing the composite norm
outliers: a text file containing 0-based outlier indices
detrend_poly: number of polynomials to add to detrend

Returns
-------
components_file: a text file containing all the regressors
"""
#builds regressor set with motion parameters, composite norm, and outliers
#outliers added as a single time point column for each outlier
def build_filter1(motion_params, comp_norm, outliers, detrend_poly = None):
    out_files = []
    for idx, filename in enumerate(filename_to_list(motion_params)): #iterate through motion params list
        params = np.genfromtxt(filename) #generate current motion params from .txt
        norm_val = np.genfromtxt(filename_to_list(comp_norm)[idx]) #generate parallel compositve norm file
        out_params = np.hstack((params, norm_val[:, None])) #horizontally stack params and norm_val, while adding second axis to norm_val
        if detrend_poly:
            timepoints = out_params.shape[0] 
            X = np.ones((timepoints, 1))
            for i in range(detrend_poly):
                X = np.hstack((X, np.array(legendre(i + 1)(np.linspace(-1, 1, timepoints)))[:, None]))
            out_params = np.hstack((out_params, X))
        try:
            outlier_val = np.genfromtxt(filename_to_list(outliers)[idx])
        except IOError:
            outlier_val = np.empty((0))
        for index in np.atleast_1d(outlier_val):
            outlier_vector = np.zeros((out_params.shape[0], 1))
            outlier_vector[int(index)] = 1 #python 2 --> 3 required explicit integer
            out_params = np.hstack((out_params, outlier_vector))
        filename = os.path.join(os.getcwd(), "filter_regressor%02d.txt" % idx)
        np.savetxt(filename, out_params, fmt = "%.10f")
        out_files.append(filename)
    return out_files


"""
Parameters
----------
realigned_file: a 4D Nifti file containing realigned volumes
mask_file: a 3D Nifti file containing white matter + ventricular masks
num_components: number of components to use for noise decomposition
extra_regressors: additional regressors to add

Returns
----------
components_file: a text file containing the noise components
"""
def extract_noise_components(realigned_file, mask_file, num_components = 5, extra_regressors = None): #derive components most reflective of physiological noise
    imgseries = nb.load(realigned_file)
    components = None
    for filename in filename_to_list(mask_file):
        mask = nb.load(filename).get_data()
        if len(np.nonzero(mask > 0)[0]) == 0:
            continue
        voxel_timecourses = imgseries.get_data()[mask > 0]
        voxel_timecourses[np.isnan(np.sum(voxel_timecourses, axis = 1)), :] = 0
        #remove mean and normalize by variance
        X = voxel_timecourses.T
        stdX = np.std(X, axis = 0)
        stdX[stdX == 0] = 1.
        stdX[np.isnan(stdX)] = 1.
        stdX[np.isinf(stdX)] = 1.
        X = (X - np.mean(X, axis = 0))/stdX
        u, _, _ = sp.linalg.svd(X, full_matrices = False)
        if components is None:
            components = u[:, :num_components]
        else:
            components = np.hstack((components, u[:, :num_components]))
    if extra_regressors:
        regressors = np.genfromtxt(extra_regressors)
        components = np.hstack((components, regressors))
    components_file = os.path.join(os.getcwd(), 'noise_components.txt')
    np.savetxt(components_file, components, fmt = "%.10f")
    return components_file


"""
Parameters
----------
files: list of 4d nifti files
lowpass_freq: cutoff frequency for the low pass filter (in Hz)
highpass_freq: cutoff frequency for the high pass filter (in Hz)
fs: sampling rate (in Hz)
"""
def bandpass_filter(files, lowpass_freq, highpass_freq, fs): #bandpass filter input files
    out_files = []
    for filename in filename_to_list(files):
        path, name, ext = split_filename(filename)
        out_file = os.path.join(os.getcwd(), name + '_bp' + ext)
        img = nb.load(filename)
        timepoints = img.shape[-1]
        F = np.zeros((timepoints))
        lowidx = timepoints//2 + 1 #python 2 --> 3 required integer division (default is float)
        if lowpass_freq > 0:
            lowidx = int(np.round(lowpass_freq / fs * timepoints)) #python 2 --> 3 required explicit integer
        highidx = 0
        if highpass_freq > 0:
            highidx = int(np.round(highpass_freq / fs * timepoints)) #python 2 --> 3 required explicit integer
        F[highidx:lowidx] = 1
        F = ((F + F[: :-1]) > 0).astype(int)
        data = img.get_data()
        if np.all(F == 1):
            filtered_data = data
        else:
            filtered_data = np.real(np.fft.ifftn(np.fft.fftn(data) * F))
        img_out = nb.Nifti1Image(filtered_data, img.get_affine(), img.get_header())
        img_out.to_filename(out_file)
        out_files.append(out_file)
    return list_to_filename(out_files)


def getmeanscale(medianvals): #get scale value to set grand mean of timeseries ~10000
    return ['-mul %.10f' % (10000. / val) for val in medianvals]


def getbtthresh(medianvals): #brightness threshold for SUSAN
    return [0.75 * val for val in medianvals]


def getusans(inlist): #return usans at right threshold
    return [[tuple([val[0], 0.75 * val[1]])] for val in inlist]


def calc_fslbp_sigmas(tr, highpass_freq, lowpass_freq): #return highpass and lowpass sigmas for fslmaths -bptf filter
    if highpass_freq < 0:
        highpass_sig = -1
    else:
        highpass_sig = 1 / (2 * tr * highpass_freq)
    if lowpass_freq < 0:
        lowpass_sig = -1
    else:
        lowpass_sig = 1 / (2 * tr * lowpass_freq)
    return highpass_sig, lowpass_sig 


def chooseindex(fwhm):
    if fwhm < 1:
        return [0]
    else:
        return [1]


def get_aparc_aseg(files):
    for name in files:
        if 'aparc+aseg' in name:
            return name
    raise ValueError('aparc+aseg.mgz not found')


tolist = lambda x: [x]

highpass_operand = lambda x: '-bptf {} {}'.format(x[0], x[1])


################################
#### Preprocessing Workflow ####
################################

def create_workflow(func_runs,
                    subject_id,
                    subjects_dir,
                    fwhm,
                    slice_times,
                    highpass_frequency,
                    lowpass_frequency,
                    TR,
                    sink_directory,
                    use_fsl_bp,
                    num_components,
                    whichvol,
                    name = 'vcat'):
    
    wf = Workflow(name = name)

    #define outputs for preprocessing workflow
    output_fields = ['reference', #extracted reference volume
                     'motion_parameters',
                     'motion_parameters_plusDerivs', 
                     'motionandoutlier_noise_file',
                     'noise_components',
                     'realigned_files', #realigned func volumes
                     'motion_plots',
                     'mask_file',
                     'smoothed_files',
                     'bandpassed_files',
                     'reg_file', #freesurfer registration file
                     'reg_cost', #freesurfer registration minumum cost file
                     'reg_fsl_file', #FLIRT-style registration file
                     'artnorm_files',
                     'artoutlier_files',
                     'artdisplacement_files',
                     'tsnr_file'] #TSNR img files   
    outputnode = Node(IdentityInterface(fields = output_fields),
                      name = 'outputspec')


    datasource = Node(DataGrabber(infields = ['subject_id', 'subject_id', 'run'], outfields = ['func']),
                      name = 'datasource')
    datasource.inputs.subject_id = subject_id
    datasource.inputs.run = func_runs
    datasource.inputs.template = '/home/data/madlab/Mattfeld_vCAT/dset/%s/func/%s_task-ROI_loc-%02d.nii.gz' #INPUT: func files
    datasource.inputs.sort_filelist = True
   

    #rename files in case identically named
    name_unique = MapNode(Rename(format_string = 'vcat_%(run)02d'), #new filename format
                          iterfield = ['in_file', 'run'],
                          name = 'rename')
    name_unique.inputs.keep_ext = True
    name_unique.inputs.run = func_runs
    wf.connect(datasource, 'func', name_unique, 'in_file') ##INPUT: datasource(func files) ----> name_unique(files to be renamed)

    
    #convert functional images to float representation
    img2float = MapNode(ImageMaths(out_data_type = 'float', op_string = '', suffix = '_dtype'),
                        iterfield = ['in_file'],
                        name = 'img2float')
    wf.connect(name_unique, 'out_file', img2float, 'in_file') ##INPUT: name_unique(renamed nii.gz files) ----> img2float(nii.gz imgs to be converted)

    

    #run AFNI's despike. This is always run, however, whether fed to realign depends on input configuration
    despiker = MapNode(Despike(outputtype = 'NIFTI_GZ'), #output as nii.gz files
                       iterfield = ['in_file'],
                       name = 'despike')
    num_threads = 4 # setting number of threads (2 cores and 4 threads -- every 2 threads share same core)
    despiker.inputs.environ = {'OMP_NUM_THREADS': '%d' % num_threads} #slurm threads
    despiker.plugin_args = {'sbatch_args': '--ntasks %d' % num_threads} #slurm  
    despiker.plugin_args = {'sbatch_args': '--tasks-per-node=1'} #slurm tasks per node
    wf.connect(img2float, 'out_file', despiker, 'in_file') ##INPUT: img2float(converted float func nii.gz imgs) ----> despiker(float func nii.gz imgs to be despiked)

    

    #extract first volume of first run as reference 
    extractref = Node(ExtractROI(t_size = 1), #single volume
                      iterfield = ['in_file'],
                      name = "extractref")
    wf.connect(despiker, ('out_file', pickfirst), extractref, 'in_file') ##INPUT: despiker(despiked float func nii.gz) ----> extractref(nii.gz file containing vols)
    wf.connect(despiker, ('out_file', pickvol, 0, whichvol), extractref, 't_min') ##INPUT: despiker(selected vol position) ----> extractref(timepoint to extract)
    wf.connect(extractref, 'roi_file', outputnode, 'reference') ##OUTPUT: extractref(reference nii.gz) ----> outputnode

    if slice_times is not None:    
        motion_correct = Node(SpaceTimeRealigner(), #simultaneous motion/slice timing correction with Nipy algorithm
                              name = 'motion_correct')
        motion_correct.inputs.tr = TR
        motion_correct.inputs.slice_times = slice_times
        motion_correct.inputs.slice_info = 2
        motion_correct.plugin_args = {'sbatch_args': '-n %s' % os.environ['MKL_NUM_THREADS']} #slurm assignment of nodes
        motion_correct.plugin_args = {'sbatch_args': '--tasks-per-node=1'} #slurm assignment of tasks per node
        wf.connect(despiker, 'out_file', motion_correct, 'in_file') ##INPUT: despiker(despiked float func nii.gz) ----> motion_correct(func files to realign)
        wf.connect(motion_correct, 'par_file', outputnode, 'motion_parameters') ##OUTPUT: motion_correct(motion parameter files) ----> outputnode
        wf.connect(motion_correct, 'out_file', outputnode, 'realigned_files') ##OUTPUT: motion_correct(realigned func files) ----> outputnode
    else:       
        motion_correct =  MapNode(MCFLIRT(save_mats = True, #correct functional runs to reference img (1st vol/1st run)
                                              save_plots = True, 
                                              interpolation = 'sinc'),
                                     name = 'motion_correct',
                                     iterfield = ['in_file'])
        wf.connect(despiker, 'out_file', motion_correct, 'in_file') ##INPUT: despiker(despiked float func nii.gz) ----> motion_correct(func files to realign)
        wf.connect(extractref, 'roi_file', motion_correct, 'ref_file') ##INPUT: extractref(reference volume) ----> motion_correct(target img for correction)
        wf.connect(motion_correct, 'par_file', outputnode, 'motion_parameters') ##OUTPUT: motion_correct(motion parameter files) ----> outputnode
        wf.connect(motion_correct, 'out_file', outputnode, 'realigned_files') ##OUTPUT: motion_correct(realigned func files) ----> outputnode


        
    #compute time-course SNR on realigned data time series
    tsnr = MapNode(TSNR(regress_poly = 2), #polynomial regression limit (2nd order)
                   iterfield = ['in_file'], 
                   name = 'tsnr')
    wf.connect(motion_correct, 'out_file', tsnr, 'in_file') ##INPUT: motion_correct(realigned func files) ----> tsnr(time-series for SNR calculation)
    wf.connect(tsnr, 'tsnr_file', outputnode, 'tsnr_file') ##OUTPUT: tsnr(TSNR image files) ----> outputnode

    

    #use fsl_tsplot to plot estimated motion parameters for realigned func data
    plot_motion = MapNode(PlotMotionParams(in_source = 'fsl'), #program generating img (FSL/SPM)
                          name = 'plot_motion',
                          iterfield = ['in_file'])
    plot_motion.iterables = ('plot_type', ['rotations', 'translations']) #FSL order (*SPM reverses order)
    wf.connect(motion_correct, 'par_file', plot_motion, 'in_file') ##INPUT: motion_correct(motion parameter files) ----> plot_motion(file containing motion params)
    wf.connect(plot_motion, 'out_file', outputnode, 'motion_plots') ##OUTPUT: plot_motion(motion plot img) ----> outputnode

    

    #register source file to freesurfer space & create brain mask in source space 
    #generates fs info from respective directories (metric fuckton of output file types)
    fssource = Node(FreeSurferSource(), 
                    name = 'fssource')
    fssource.inputs.subjects_dir = subjects_dir
    fssource.inputs.subject_id = subject_id
    

    
    #extract and binarize aparc+aseg brain mask 
    fs_threshold = Node(Binarize(min = 0.5, #minimum threshold
                                 out_type = 'nii'), 
                        name = 'fs_threshold')
    wf.connect(fssource, ('aparc_aseg', get_aparc_aseg), fs_threshold, 'in_file') ##INPUT: fssource(aparc/aseg file) ----> fs_threshold(volume to binarize)


    
    #calculate transformation matrix from EPI space to FreeSurfer anatomical space
    fs_register = MapNode(BBRegister(init = 'fsl'),
                          iterfield = ['source_file'],
                          name = 'fs_register')
    fs_register.inputs.contrast_type = 't2' #which contrast type to 
    fs_register.inputs.out_fsl_file = True #write in FSL format
    fs_register.inputs.subject_id = subject_id
    fs_register.inputs.subjects_dir = subjects_dir
    wf.connect(extractref, 'roi_file', fs_register, 'source_file') ##INPUT: extractref(reference volume) ----> fs_register(file to be registered to FS space)
    wf.connect(fs_register, 'out_reg_file', outputnode, 'reg_file') ##OUTPUT: fs_register(registration file) ----> outputnode
    wf.connect(fs_register, 'min_cost_file', outputnode, 'reg_cost') ##OUTPUT: fs_register(registration minimum cost file) ----> outputnode
    wf.connect(fs_register, 'out_fsl_file', outputnode, 'reg_fsl_file') ##OUTPUT: fs_register(FLIRT-style registration file) ----> outputnode

    

    #extract white matter (wm) + cerebral spinal fluid (csf) brain masks by eroding freesurfer labels
    wmcsf = MapNode(Binarize(), 
                    iterfield = ['match', 'binary_file', 'erode'], 
                    name = 'wmcsfmask')
    wmcsf.inputs.match = [[2, 41], [4, 5, 14, 15, 24, 31, 43, 44, 63]] #match labels instead of using a threshold
    wmcsf.inputs.binary_file = ['wm.nii.gz', 'csf.nii.gz'] #binary output volumes
    wmcsf.inputs.erode = [2, 2] #erode binarization in 3D
    wf.connect(fssource, ('aparc_aseg', get_aparc_aseg), wmcsf, 'in_file') #INPUT: fssource(apar/aseg file) ----> wmcsf(volume to binarize)



    #transform wm and csf masks from FS to EPI space using reference volume
    wmcsftransform = MapNode(ApplyVolTransform(inverse = True, #sample from target source
                                               interp = 'nearest'), #interpolation method
                             iterfield = ['target_file'],
                             name = 'wmcsftransform')
    wmcsftransform.inputs.subjects_dir = subjects_dir
    wf.connect(extractref, 'roi_file', wmcsftransform, 'source_file') #INPUT: extractref(reference file) ----> wmcsftransform(volume space/type to match)
    wf.connect(fs_register, ('out_reg_file', pickfirst), wmcsftransform, 'reg_file') #INPUT: fs_register(fs registration file) ----> wmcsftransform(tkregister matrix)
    wf.connect(wmcsf, 'binary_file', wmcsftransform, 'target_file') #INPUT: wmcsf(binarized wmcsf mask) ----> wmcsftransform(volume to transform)



    #transform binarized aparc+aseg file from FS to EPI space using reference volume
    fs_voltransform = MapNode(ApplyVolTransform(inverse = True), #sample from target source
                              iterfield = ['source_file', 'reg_file'],
                              name = 'fs_transform')
    fs_voltransform.inputs.subjects_dir = subjects_dir
    wf.connect(extractref, 'roi_file', fs_voltransform, 'source_file') #INPUT: extractref(reference file) ----> fs_voltransform(volume space/type to match)
    wf.connect(fs_register, 'out_reg_file', fs_voltransform, 'reg_file') #INPUT: fs_register(fs registration file) ----> fs_voltransform(tkregister matrix)
    wf.connect(fs_threshold, 'binary_file', fs_voltransform, 'target_file') #INPUT: wmcsf(binarized wmcsf mask) ----> fs_voltransform(volume to transform)



    #dilate binarized mask (that is now in EPI space) by 1 voxel 
    fs_threshold2 = MapNode(Binarize(min = 0.5, #minimum threshold
                                     out_type = 'nii'),
                            iterfield = ['in_file'],
                            name = 'fs_threshold2')
    fs_threshold2.inputs.dilate = 1 #dilate binarization in 3D
    wf.connect(fs_voltransform, 'transformed_file', fs_threshold2, 'in_file') #INPUT: fs_voltransform(EPI space- binarized mask) ----> fs_threshold2(mask to dilate)
    wf.connect(fs_threshold2, 'binary_file', outputnode, 'mask_file') #OUTPUT: fs_threshold2(dilated binarized aparc/aseg mask) ----> outputnode


    
    #use RapidART to detect motion/intensity outliers
    art = MapNode(ArtifactDetect(use_differences = [True, False], #use differences between successive motion (1st) and intensity params (2nd) estimates to find outliers
                                 use_norm = True, #use composite of motion params to determine outliers
                                 zintensity_threshold = 3, #intensity z-threshold to detect images that deviate from mean
                                 norm_threshold = 1, #threshold to detect motion-related outliers when composite motion is used
                                 bound_by_brainmask = True, #use the brainmask to determine bounding box for composite norm
                                 mask_type = "file"), #type used to mask functional data
                  iterfield = ["realignment_parameters","realigned_files"], 
                  name = "art")
    if slice_times is not None: #if slice-timing correction is performed
        art.inputs.parameter_source = "NiPy"
    else:
        art.inputs.parameter_source = "FSL"
    wf.connect(motion_correct, 'par_file', art, 'realignment_parameters') ##INPUT: motion_correct(parameters file) ----> art(realignment parameter file)
    wf.connect(motion_correct, 'out_file', art, 'realigned_files') ##INPUT: motion_correct(realigned func files) ----> art(files for detection of outliers)
    wf.connect(fs_threshold2, ('binary_file', pickfirst), art, 'mask_file') ##INPUT: fs_threshold2(dilated binarized mask file) ----> art(mask file to be used)
    wf.connect(art, 'norm_files', outputnode, 'artnorm_files') #OUTPUT: art(file containing composite norm) ----> outputnode
    wf.connect(art, 'outlier_files', outputnode, 'artoutlier_files') #OUTPUT: art(list of 0-based indices corresponding to outlier volumes) ----> outputnode
    wf.connect(art, 'displacement_files', outputnode, 'artdisplacement_files') #OUTPUT: art(file containing the voxel displacement timeseries) ----> outputnode



    #compute motion regressors (save file w/ 1st & 2nd derivatives)
    motreg = Node(Function(input_names = ['motion_params', 'order', 'derivatives'],
                           output_names = ['out_files'],
                           function = motion_regressors, #see function to details
                           imports = imports),
                  name = 'getmotionregress')
    wf.connect(motion_correct, 'par_file', motreg, 'motion_params') ##INPUT: motion_correct(motion parameters) ----> motreg(params to be combined w/ derivatives) 
    wf.connect(motreg, 'out_files', outputnode, 'motion_parameters_plusDerivs') #OUTPUT: motreg(combined motion + derivatives file) ----> outputnode



    #create filter text file to remove motion (+ derivatives), art confounds, and 1st, 2nd, & 3rd order legendre polynomials
    createfilter1 = Node(Function(input_names = ['motion_params', 'comp_norm', 'outliers', 'detrend_poly'],
                                  output_names = ['out_files'],
                                  function = build_filter1, #builds regressor set with motion parameters, composite norm, and outliers
                                  imports = imports), #unique imports required within function
                         name = 'makemotionbasedfilter')
    createfilter1.inputs.detrend_poly = 3 #determines what order of polynomials to include
    wf.connect(motreg, 'out_files', createfilter1, 'motion_params') ##INPUT: motreg() ----> createfilter1()
    wf.connect(art, 'norm_files', createfilter1, 'comp_norm') ##INPUT: art(file containing composite norm) ----> createfilter1(composite norm list)
    wf.connect(art, 'outlier_files', createfilter1, 'outliers') ##INPUT: art(list of 0-based indices corresponding to outlier volumes) ----> createfilter1(outlier list)
    wf.connect(createfilter1, 'out_files', outputnode, 'motionandoutlier_noise_file') #OUTPUT: createfilter1(complete regressor set) -----> outputnode



    #create filter to remove noise components based on white matter and CSF
    createfilter2 = MapNode(Function(input_names = ['realigned_file', 'mask_file', 'num_components', 'extra_regressors'],
                                     output_names = ['out_files'],
                                     function = extract_noise_components, #derive components most reflective of physiological noise
                                     imports = imports), #unique imports required within function
                            iterfield = ['realigned_file', 'extra_regressors'],
                            name = 'makecompcorrfilter')
    createfilter2.inputs.num_components = num_components
    wf.connect(createfilter1, 'out_files', createfilter2, 'extra_regressors') ##INPUT: createfilter1(complete regressor set) ----> createfilter2()
    wf.connect(motion_correct, 'out_file', createfilter2, 'realigned_file') ##INPUT: motion_correct(realigned func files) ----> createfilter2(realigned files)
    wf.connect(wmcsftransform, 'transformed_file', createfilter2, 'mask_file') ##INPUT: wmcsftransform(transformed wmcsf volume) ----> createfilter2(mask)
    wf.connect(createfilter2, 'out_files', outputnode, 'noise_components') #OUTPUT: createfilter2(filter to remove wmcsf noise) ----> outputnode



    #mask functional runs with extracted mask
    maskfunc = MapNode(ImageMaths(suffix = '_bet', op_string = '-mas'), #-mas = use to mask current image
                       iterfield = ['in_file'],
                       name = 'maskfunc')
    wf.connect(motion_correct, 'out_file', maskfunc, 'in_file') ##INPUT: motion_correct(corrected func files) ----> maskfunc(files to mask)
    wf.connect(fs_threshold2, ('binary_file', pickfirst), maskfunc, 'in_file2') ##INPUT: fs_threshold2(first volume extracted mask) ----> maskfunc(mask to use)
 

   
    #smooth each run using SUSAN with brightness threshold set to 50% of median value for each run and mask constituting mean functional
    smooth_median = MapNode(ImageStats(op_string = '-k %s -p 50'), #use image for masking, output 50th percentile
                            iterfield = ['in_file'],
                            name = 'smooth_median')
    wf.connect(maskfunc, 'out_file', smooth_median, 'in_file') ##INPUT: maskfunc(masked func runs) ----> smooth_median(infile to perform stats upon)
    wf.connect(fs_threshold2, ('binary_file', pickfirst), smooth_median, 'mask_file') ##INPUT: fs_threshold2(1st binarized volume) ----> smooth_median(mask file -k)
    


    #obtain temporal mean of maskfunc
    smooth_meanfunc = MapNode(ImageMaths(op_string = '-Tmean',  suffix = '_mean'), #calculate mean across time
                              iterfield = ['in_file'],
                              name = 'smooth_meanfunc')
    wf.connect(maskfunc, 'out_file', smooth_meanfunc, 'in_file') ##INPUT: maskfunc(masked func runs) ----> smooth_meanfunc(image to be temporally averaged)


    
    #merge smooth median and smooth temporal mean
    smooth_merge = Node(Merge(2, axis = 'hstack'),
                        name = 'smooth_merge')
    wf.connect(smooth_meanfunc, 'out_file', smooth_merge, 'in1') ##INPUT:
    wf.connect(smooth_median, 'out_stat', smooth_merge, 'in2') ##INPUT:


    #FSL noise reduction algorithm using nonlinear filtering 
    smooth = MapNode(SUSAN(),
                     iterfield = ['in_file', 'brightness_threshold', 'usans'],
                     name = 'smooth')
    smooth.inputs.fwhm = fwhm
    wf.connect(maskfunc, 'out_file', smooth, 'in_file') ##INPUT:
    wf.connect(smooth_median, ('out_stat', getbtthresh), smooth, 'brightness_threshold') ##INPUT:
    wf.connect(smooth_merge, ('out', getusans), smooth, 'usans') ##INPUT:
 

   
    #mask smoothed data with dilated mask
    maskfunc2 = MapNode(ImageMaths(suffix = '_mask', op_string = '-mas'),
                        iterfield = ['in_file'],
                        name = 'maskfunc2')
    wf.connect(smooth, 'smoothed_file', maskfunc2, 'in_file') ##INPUT:
    wf.connect(fs_threshold2, ('binary_file', pickfirst), maskfunc2, 'in_file2') ##INPUT:
    wf.connect(maskfunc2, 'out_file', outputnode, 'smoothed_files') #OUTPUT:



    #band-pass filter timeseries
    if use_fsl_bp == 'True':
        determine_bp_sigmas = Node(Function(input_names = ['tr', 'highpass_freq', 'lowpass_freq'],
                                            output_names = ['out_sigmas'],
                                            function = calc_fslbp_sigmas),
                                   name = 'determine_bp_sigmas')
        determine_bp_sigmas.inputs.tr = float(TR)
        determine_bp_sigmas.inputs.highpass_freq = float(highpass_frequency)
        determine_bp_sigmas.inputs.lowpass_freq = float(lowpass_frequency)

        bandpass = MapNode(ImageMaths(suffix = '_tempfilt'),
                           iterfield = ["in_file"],
                           name = "bandpass")
        wf.connect(determine_bp_sigmas, ('out_sigmas', highpass_operand), bandpass, 'op_string') ##INPUT:
        wf.connect(maskfunc2, 'out_file', bandpass, 'in_file') ##INPUT:
        wf.connect(bandpass, 'out_file', outputnode, 'bandpassed_files') #OUTPUT:
    else:
        bandpass = Node(Function(input_names = ['files', 'lowpass_freq', 'highpass_freq', 'fs'],
                                 output_names = ['out_files'],
                                 function = bandpass_filter,
                                 imports = imports),
                           name = 'bandpass')
        bandpass.inputs.fs = 1./TR
        if highpass_frequency < 0:
            bandpass.inputs.highpass_freq = -1
        else:
            bandpass.inputs.highpass_freq = highpass_frequency
        if lowpass_frequency < 0:
            bandpass.inputs.lowpass_freq = -1
        else:
            bandpass.inputs.lowpass_freq = lowpass_frequency
        wf.connect(maskfunc2, 'out_file', bandpass, 'files') ##INPUT:
        wf.connect(bandpass, 'out_files', outputnode, 'bandpassed_files') #OUTPUT:

        
        
    #save relevant data into output directory
    datasink = Node(DataSink(), 
                    name = "datasink")
    datasink.inputs.base_directory = sink_directory
    datasink.inputs.container = subject_id
    wf.connect(outputnode, 'reference', datasink, 'ref')
    wf.connect(outputnode, 'motion_parameters', datasink, 'motion')
    wf.connect(outputnode, 'realigned_files', datasink, 'func.realigned')
    wf.connect(outputnode, 'motion_plots', datasink, 'motion.@plots')
    wf.connect(outputnode, 'mask_file', datasink, 'ref.@mask')
    wf.connect(outputnode, 'smoothed_files', datasink, 'func.smoothed_fullspectrum')
    wf.connect(outputnode, 'bandpassed_files', datasink, 'func.smoothed_bandpassed')
    wf.connect(outputnode, 'reg_file', datasink, 'bbreg.@reg')
    wf.connect(outputnode, 'reg_cost', datasink, 'bbreg.@cost')
    wf.connect(outputnode, 'reg_fsl_file', datasink, 'bbreg.@regfsl')
    wf.connect(outputnode, 'artnorm_files', datasink, 'art.@norm_files')
    wf.connect(outputnode, 'artoutlier_files', datasink, 'art.@outlier_files')
    wf.connect(outputnode, 'artdisplacement_files', datasink, 'art.@displacement_files')
    wf.connect(outputnode, 'motion_parameters_plusDerivs', datasink, 'noise.@motionplusDerivs')
    wf.connect(outputnode, 'motionandoutlier_noise_file', datasink, 'noise.@motionplusoutliers')
    wf.connect(outputnode, 'noise_components', datasink, 'compcor')
    wf.connect(outputnode, 'tsnr_file', datasink, 'tsnr')    

    return wf

"""
Creates the full workflow; gets information from dicom files
"""

def create_preproc_workflow(args, name = 'vcat'):
    import numpy as np
    TR = args.tr
    if args.do_slice_times == 'True':
       n_slices = 66
       slice_order = list(range(0, n_slices, 2)) + list(range(1, n_slices, 2))
       slice_order = np.argsort(slice_order)
       sliceTimes = (slice_order * TR/n_slices).tolist()
    else:
       sliceTimes = None

    kwargs = dict(func_runs = list(map(int, args.func_runs)),
                  subject_id = args.subject_id,
                  subjects_dir = args.subjects_dir,
                  fwhm = args.fwhm,
                  slice_times = sliceTimes,
                  highpass_frequency = args.highpass_freq,
                  lowpass_frequency = args.lowpass_freq,
                  TR = TR,
                  sink_directory = os.path.abspath(args.sink),
                  use_fsl_bp = args.use_fsl_bp,
                  num_components = int(args.num_components),
                  whichvol = args.whichvol,
                  name = name)
    wf = create_workflow(**kwargs)
    return wf

if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser(description = __doc__)
    parser.add_argument("-r", "--runs", dest = "func_runs", nargs = "+", help = "List of functional runs to preprocess", required = True)
    parser.add_argument("-s", "--subject_id", dest = "subject_id", help = "FreeSurfer subject id", required = True)
    parser.add_argument("-d", "--subjects_dir", dest = "subjects_dir", help = "FreeSurfer subject dir", required = True)
    parser.add_argument("-u", "--highpass_freq", dest = "highpass_freq", default = -1, type = float, help = "High pass frequency in Hz")
    parser.add_argument("-l", "--lowpass_freq", dest = "lowpass_freq", default = -1, type = float, help = "Low pass frequency in Hz")
    parser.add_argument("--do_slice_times", dest = "do_slice_times", default = False, help = "Slice times in seconds")
    parser.add_argument("--use_fsl_bp", dest = "use_fsl_bp", default = True, help = "Bandpass filter using FSLmaths")
    parser.add_argument("-t", "--repetition_time", dest = "tr", default = 1.76, type = float, help = "TR of functional data in seconds")
    parser.add_argument("-k", "--fwhm_kernel", dest = "fwhm", default = 5.0, type = float, help = "FWHM smoothing kernel (mm)")
    parser.add_argument('-n', "--num_comp", dest = "num_components", default = 3, help = "Number of components to extract for aCompCor")
    parser.add_argument('-v', "--which_vol", dest = "whichvol", default = 'first', help = "Which volume for ref registration")
    parser.add_argument("-o", "--output_dir", dest = "sink", help = "Output directory base")
    parser.add_argument("-w", "--work_dir", dest = "work_dir", help = "Output directory base")
    args = parser.parse_args()

    wf = create_preproc_workflow(args)

    if args.work_dir:
        work_dir = os.path.abspath(args.work_dir)
    else:
        work_dir = os.getcwd()

    wf.config['execution']['crashdump_dir'] = '/scratch/madlab/vcat/crash'
    wf.base_dir = work_dir + '/' + args.subject_id
    wf.run(plugin='SLURM', plugin_args={'sbatch_args': ('-p IB_16C_96G --qos pq_madlab --account iacc_madlab'), 'overwrite': True})    
