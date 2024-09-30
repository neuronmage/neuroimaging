#!/bin/bash


echo ""



if [ $# -lt 2 ]
then
    echo "ABCD processing script that combines several features from different pipelines. "
    echo "Usage: ABCD_DWI_processing options"
    echo ""
    echo "DATA INPUTS: (no space after the = sign)"
    echo "    -u=full_path_to_the_main_NIFTI_file   OR --up_data=full_path_to_the_main_NIFTI_file (required)"
    echo "    -ub=main_data_bvals_file  OR --up_bvals=main_data_bvals_file (OPTIONAL). If not provided, the script will use bvals file from the main data folder"
    echo "    -uv=main_data_bvecs_file  OR --up_bvecs=main_data_bvecs_file (OPTIONAL). If not provided, the script will use bvecs file from the main data folder"    
    echo "    -d=full_path_to_the_reverse_PE_NIFTI_file   OR --down_data=full_path_to_the_reverse_PE_NIFTI_file (required)"
    echo "    -db=down_data_bvals_file  OR --down_bvals=down_data_bvals_file (OPTIONAL). If not provided, the script will create a dummy bvals file with zeros"
    echo "    -dv=down_data_bvecs_file  OR --down_bvecs=down_data_bvecs_file (OPTIONAL). If not provided, the script will create a dummy bvecs file with zeros"    


    echo " "
    echo "    -T1=full_path_to_T1W_structural_image  (REQUIRED)"
    echo "    -T2=full_path_to_T2W_structural_image  (OPTIONAL)"
    echo ""
    echo "    --gw_coeffs=full_path_to_gradient_nonlinearity_coefficients_file (OPTIONAL). If provided, gradwarp will be applied and the gradient deviation matrix image L will be output in HCP style"
    echo "    --gw_field_full=full_path_to_3D_gradient_nonlinearity_field in ITK displacement field format. (OPTIONAL)"
    echo "    --gw_field_L=full_path_to_gradient_nonlinearity_field_L_fiel_in_MGH_format (OPTIONAL). All three must be provided if one of them is."
    echo "    --gw_field_P=full_path_to_gradient_nonlinearity_field_P_fiel_in_MGH_format (OPTIONAL). All three must be provided if one of them is."
    echo "    --gw_field_H=full_path_to_gradient_nonlinearity_field_H_fiel_in_MGH_format (OPTIONAL). All three must be provided if one of them is."
    echo "    --gw_mode=2d/3d/1d . (OPTIONAL). Flag to determine if the gradwarp will be applied in 2D or 3D or 1d (which will only be applied along the Z axis). Default:3d"
    echo "    --isocenter=0 or 1 .  (OPTIONAL). Default:0 . When 1, the image center is assumed to be the isocenter of the magnet. When 0, the image header is used for the z axis as well."

    
    echo ""
    echo "TEMPORARY DATA options:"
    echo "     -p=temporary_processing_folder OR --proc_dir=temporary_processing_folder. (OPTIONAL). If not provided, the script will create a proc folder in the folder containing the up NIFTI file"


    echo ""
    echo "OUTPUT options:"
    echo "     -o=output_folder OR --output_dir=output_folder (OPTIONAL). If not provided, the script will create a final folder in the folder containing the up NIFTI file."
    echo "    --output_res=final_DWI_and_structural_resolution    (OPTIONAL. Default:native_resolution. Example usage: --output_res=1.5 (for isotropic) OR --output_res=1.5x1.5x3 (for anisotropic))"
    echo "    --output_size=final_DWI_and_structural_size    (OPTIONAL. Default:Size computed from structural FOV and desired resolution. Example usage: --output_size=150x160x140 "
    echo "    --output_orientation={LPS,RAI,LAI, etc..}  (Default:LPS , i.e. from Right-to-Left (from left of the image to right), Anterior-to-Posterior (from top of the image to bottom), Inferior-to-superior"

                  
    echo ""
    echo " PROCESSING OPTIONS: (no space after the = sign)"
    echo "    --denoising={0,1} (OPTIONAL. Default:0)"
    echo "    --gibbs={0,1}     (OPTIONAL. Default:0)"
    echo "    --use_T2={0,1}    (OPTIONAL. Default:1. T2W image, if present, will be used along with T1W for susceptibility distortion correction)"
    
    
    echo " By M. Okan Irfanoglu 10/20/2021"
    echo " "
    
    exit
fi

cwd=`pwd`
script_folder=$( dirname $( realpath $0 ))
echo SCRIPT FOLDER:  ${script_folder}

export FSLOUTPUTTYPE=NIFTI_GZ
EDDY_EXECUTABLE=eddy_openmp
chmod -R a+x ${script_folder}

module load singularity-3.5.3 #NOTE ATM 2.18.22: TRY singularity-3.5.3  and then try singularity-3
export CUDA_VISIBLE_DEVICES=""


export PATH=${script_folder}/ANTS:${PATH}
export PATH=${script_folder}/FREESURFER:${PATH}
export PATH=${script_folder}/FSL:${PATH}
export PATH=${script_folder}/SYNB0-DISCO:${PATH}
export PATH=${script_folder}/TORTOISE/DIFFCALC:${PATH}
export PATH=${script_folder}/TORTOISE/DIFFPREP/bin/bin:${PATH}
export PATH=${script_folder}/TORTOISE/DRBUDDI/bin:${PATH}

Ncores=`nproc --all`
export OMP_NUM_THREADS=${Ncores}


if [ -z "${LD_LIBRARY_PATH}"  ]  
then
    export LD_LIBRARY_PATH=${script_folder}/FSL/lib
else
    export LD_LIBRARY_PATH=${script_folder}/FSL/lib:${LD_LIBRARY_PATH}
fi



#ALL THE MAIN VARIABLES

up_nifti=''
up_bvals=''
up_bvecs=''

down_nifti=''
down_bvals=''
down_bvecs=''

T1W_structural=''
T2W_structural=''

gw_coeffs_file=''
gw_field_full=''
gw_field_L_file=''
gw_field_P_file=''
gw_field_H_file=''
gw_mode='3d'
isocenter=0

denoising=0
gibbs=0
use_T2=1
output_res=""
output_size=""
output_orientation=""

proc_dir=""
output_dir=""



#INPUT PARSING

cmd=""

for i in "$@"
do
   cmd="$cmd $i"

case $i in
    -u=*|--up_data=*)
        up_nifti="${i#*=}"
        shift 
    ;;
    -ub=*|--up_bvals=*)
        up_bvals="${i#*=}"
        shift 
    ;;   
    -uv=*|--up_bvecs=*)
        up_bvecs="${i#*=}"
        shift 
    ;;   
    -d=*|--down_data=*)
        down_nifti="${i#*=}"
        shift 
    ;;
    -db=*|--down_bvals=*)
        down_bvals="${i#*=}"
        shift 
    ;;   
    -dv=*|--down_bvecs=*)
        down_bvecs="${i#*=}"
        shift 
    ;;    
    -T1=*)
        T1W_structural="${i#*=}"
        shift 
    ;;          
    -T2=*)
        T2W_structural="${i#*=}"
        shift 
    ;;     
    --gw_coeffs=*)
        gw_coeffs_file="${i#*=}"
        shift 
    ;;  
    --gw_field_full=*)
        gw_field_full="${i#*=}"
        shift 
    ;;   
    --gw_field_L=*)
        gw_field_L_file="${i#*=}"
        shift 
    ;;   
    --gw_field_P=*)
        gw_field_P_file="${i#*=}"
        shift 
    ;;   
    --gw_field_H=*)
        gw_field_H_file="${i#*=}"
        shift 
    ;;   
    --gw_mode=*)
        gw_mode="${i#*=}"
        shift 
    ;;   
    --isocenter=*)
        isocenter="${i#*=}"
        shift 
    ;;   
    --denoising=*)
        denoising="${i#*=}"
        shift 
    ;;    
    --gibbs=*)
        gibbs="${i#*=}"
        shift 
    ;;    
    --use_T2=*)
        use_T2="${i#*=}"
        shift 
    ;;    
    --output_res=*)
        output_res="${i#*=}"
        shift 
    ;;   
    --output_size=*)
        output_size="${i#*=}"
        shift 
    ;;  
    --output_orientation=*)
        output_orientation="${i#*=}"
        shift 
    ;;  
    -p=*|--proc_dir=*)
        proc_dir="${i#*=}"
        shift 
    ;;   
    -o=*|--output_dir=*)
        output_dir="${i#*=}"
        shift 
    ;;   
    *)
        echo Unrecognized command line option ${i}.  Exiting
        exit
    ;;
    *)
            # unknown option
    ;;
esac
done


#INPUT ERROR CHECKS

if [ ! -e "${up_nifti}" ]
then
    echo Up NIFTI: ${up_nifti} does not exist. Exiting....
    exit
fi
if [ ! -e "${down_nifti}" ]
then
    echo Down NIFTI: ${down_nifti} does not exist. Exiting....
    exit
fi
if [ ! -e "${T1W_structural}" ]
then
    echo T1W_structural: ${T1W_structural} does not exist. Exiting....
    exit
fi


if [ ! -e "${T1W_structural}" ]
then
    echo T1W_structural: ${T1W_structural} does not exist. Exiting....
    exit
fi



#HANDLING OF BVECS/BVALS WHEN THEY ARE NOT ENTERED

if [ -z "${up_bvecs}"  ]  
then
    up_bvecs=$(find  $(dirname ${up_nifti})  -maxdepth 1  -name "*.bvec*")  
fi
if [ -z "${up_bvals}"  ]  
then
    up_bvals=$(find  $(dirname ${up_nifti})  -maxdepth 1  -name "*.bval*")  
fi




#CREATION OF PROCESSING FOLDERS AND FILES

up_folder=`dirname ${up_nifti}`

if [ -z "${proc_dir}"  ]  
then
   up_proc_folder=${up_folder}/proc
else
   up_proc_folder=${proc_dir}
fi

if [ ! -e "${up_proc_folder}" ]
then
    mkdir -p ${up_proc_folder}
fi


# OUTPUT THE ENTERED COMMAND TO LOG
cmd_file=${up_proc_folder}/cmd.log
if [  -e "${cmd_file}" ]
then
   rm ${cmd_file}
fi
echo ${cmd} > ${cmd_file} 





#up_ext="${up_nifti#*.}"
#down_ext="${down_nifti#*.}"


if [[ "$up_nifti" == *".nii.gz"* ]]; then
  up_ext=nii.gz
else
  up_ext=nii  
fi

if [[ "$down_nifti" == *".nii.gz"* ]]; then
  down_ext=nii.gz
else
  down_ext=nii  
fi



up_nifti_proc=${up_proc_folder}/`basename ${up_nifti} .${up_ext}`_proc.${up_ext}
echo "ABCD MESSAGE:  COPYING THE UP FILES TO THE PROC FOLDER"
cp ${up_nifti} ${up_nifti_proc}
cp ${up_bvecs} ${up_proc_folder}
cp ${up_bvals} ${up_proc_folder}

if [ ${up_ext} = "nii.gz" ]
then 
    gunzip -f ${up_nifti_proc}
fi

up_nifti_proc=${up_proc_folder}/`basename ${up_nifti} .${up_ext}`_proc.nii
up_bvecs_proc=${up_proc_folder}/`basename ${up_bvecs}`
up_bvals_proc=${up_proc_folder}/`basename ${up_bvals}`




down_nifti_proc=${up_proc_folder}/`basename ${down_nifti} .${down_ext}`_proc.${down_ext}
echo "ABCD MESSAGE:  COPYING THE DOWN FILES TO THE PROC FOLDER"
cp ${down_nifti} ${down_nifti_proc}
if [ ${down_ext} = "nii.gz" ]
then 
    gunzip -f ${down_nifti_proc}
fi
down_nifti_proc=${up_proc_folder}/`basename ${down_nifti} .${down_ext}`_proc.nii


down_nvol=`fslinfo ${down_nifti_proc}| grep ^dim4 | awk '{print $2}'`
if [ -z "${down_bvals}"  ]
then
    down_folder=`dirname ${down_nifti_proc}`
    down_temp_bval=${down_folder}/dummy.bval
    
    if [ -e "${down_temp_bval}" ]
    then
        rm -rf ${down_temp_bval}  
    fi
    
    zeros=""
    for ((i=1; i<=${down_nvol}; i+=1)); do zeros="0 ${zeros}"; done
    echo ${zeros} > ${down_temp_bval}
    
    down_bvals=${down_temp_bval}
else    
    cp ${down_bvals} ${up_proc_folder}
    down_bvals=${up_proc_folder}/`basename ${down_bvals}`   
fi



if [ -z "${down_bvecs}"  ]  
then
    down_folder=`dirname ${down_nifti_proc}`
    down_temp_bvec=${down_folder}/dummy.bvec
    
    if [ -e "${down_temp_bvec}" ]
    then
        rm -rf ${down_temp_bvec}  
    fi
    zeros=""
    for ((i=1; i<=${down_nvol}; i+=1)); do zeros="0 ${zeros}"; done
    echo ${zeros} > ${down_temp_bvec}
    echo ${zeros} >> ${down_temp_bvec}
    echo ${zeros} >> ${down_temp_bvec}
    
    down_bvecs=${down_temp_bvec}
    
else
    cp ${down_bvecs} ${up_proc_folder}
    down_bvecs=${up_proc_folder}/`basename ${down_bvecs}`
fi








#OPTIONAL  DENOSING AND GIBBS RINGING

if [ ${denoising} -eq 1 ]   
then
    echo "ABCD MESSAGE:  DENOISING UP DATA"
    DWIDenoise ${up_nifti_proc} ${up_nifti_proc} ${up_proc_folder}/`basename ${up_nifti_proc} .nii`_noise.nii      
fi

if [ ${gibbs} -eq 1 ]   
then
    echo "ABCD MESSAGE:  GIBBS RINGING CORRECTION UP DATA"
    UnRing ${up_nifti_proc} ${up_nifti_proc}
fi









#FSL EDDY

fslroi ${up_nifti_proc} ${up_proc_folder}/up_b0 0 1
bet2 ${up_proc_folder}/up_b0  ${up_proc_folder}/up_b0_mask -m -f 0.3
printf "0 -1 0 0.1\n0 1 0 0.1" >  ${up_proc_folder}/acqparams.txt


nvol=`fslinfo ${up_nifti_proc}| grep ^dim4 | awk '{print $2}'`
indx=""
for ((i=1; i<=${nvol}; i+=1)); do indx="$indx 1"; done
echo $indx > ${up_proc_folder}/index.txt

slspec=${script_folder}/my_slspec.txt

echo "ABCD MESSAGE:  FSL EDDY"
${EDDY_EXECUTABLE} --imain=${up_nifti_proc}  --mask=${up_proc_folder}/up_b0_mask_mask  --acqp=${up_proc_folder}/acqparams.txt --index=${up_proc_folder}/index.txt --bvecs=${up_bvecs_proc} --bvals=${up_bvals_proc} --out=${up_proc_folder}/`basename ${up_nifti_proc} .nii`_eddy_corrected_data --niter=8 --fwhm=10,8,4,2,0,0,0,0 --repol  --mporder=6 --slspec=${slspec} --s2v_niter=5 --s2v_lambda=1 --s2v_interp=trilinear  --data_is_shelled  --verbose --nvoxhp=50000 --initrand=1


eddy_corrected_data=${up_proc_folder}/`basename ${up_nifti_proc} .nii`_eddy_corrected_data.nii.gz
#eddy_outlier_free_data=${up_proc_folder}/`basename ${up_nifti_proc} .nii`_eddy_corrected_data.eddy_outlier_free_data.nii.gz
eddy_outlier_free_data=${eddy_corrected_data}



#EDDY SOMETIMES CRASHES WITH ABCD DATA. IF THE ABOVE RUN CRASHED, RUN IT WITH MORE RELAXED PARAMETERS
if [ ! -f "$eddy_corrected_data" ]
then
    ${EDDY_EXECUTABLE} --imain=${up_nifti_proc}  --mask=${up_proc_folder}/up_b0_mask_mask  --acqp=${up_proc_folder}/acqparams.txt --index=${up_proc_folder}/index.txt --bvecs=${up_bvecs_proc} --bvals=${up_bvals_proc} --out=${up_proc_folder}/`basename ${up_nifti_proc} .nii`_eddy_corrected_data --niter=8 --fwhm=10,8,4,2,0,0,0,0 --repol  --mporder=6 --slspec=${slspec} --s2v_niter=5 --s2v_lambda=1 --s2v_interp=trilinear  --data_is_shelled   --ol_ec=2  --verbose --nvoxhp=50000 --initrand=1
fi

#EDDY SOMETIMES CRASHES WITH ABCD DATA. IF THE ABOVE RUN CRASHED, RUN IT WITH MORE RELAXED PARAMETERS.
if [ ! -f "$eddy_corrected_data" ]
then
    ${EDDY_EXECUTABLE} --imain=${up_nifti_proc}  --mask=${up_proc_folder}/up_b0_mask_mask  --acqp=${up_proc_folder}/acqparams.txt --index=${up_proc_folder}/index.txt --bvecs=${up_bvecs_proc} --bvals=${up_bvals_proc} --out=${up_proc_folder}/`basename ${up_nifti_proc} .nii`_eddy_corrected_data --niter=8 --fwhm=10,8,4,2,0,0,0,0 --repol  --mporder=6 --slspec=${slspec} --s2v_niter=5 --s2v_lambda=10 --s2v_interp=trilinear  --data_is_shelled --ol_ec=2   --ol_nstd=5 --verbose --nvoxhp=50000 --initrand=1
fi
up_bvecs_proc=${up_proc_folder}/`basename ${eddy_outlier_free_data} .nii.gz`.eddy_rotated_bvecs



# PROCESS THE STRUCTURAL IMAGES

echo "ABCD MESSAGE:  COPYING THE T1W IMAGE"
cp ${T1W_structural} ${up_proc_folder}
T1W_structural_proc=${up_proc_folder}/`basename ${T1W_structural}`
if [[ "$T1W_structural_proc" == *".nii.gz"*  ]]
then
    gunzip -f ${T1W_structural_proc}
    T1W_structural_proc=`dirname ${T1W_structural_proc}`/`basename ${T1W_structural_proc} .nii.gz`.nii
fi    
echo "ABCD MESSAGE:  REORIENTING THE T1W IMAGE TO AXIAL"
ReorientImage3D -i ${T1W_structural_proc}
T1W_structural_proc_to_use=${up_proc_folder}/`basename ${T1W_structural_proc} .nii`_reoriented.nii
    


if [ ! -z "${T2W_structural}"  ]  
then
    echo "ABCD MESSAGE:  COPYING THE T2W IMAGE"
    cp ${T2W_structural} ${up_proc_folder}
    T2W_structural_proc=${up_proc_folder}/`basename ${T2W_structural}`
    if [[ "$T2W_structural_proc" == *".nii.gz"*  ]]
    then
        gunzip -f ${T2W_structural_proc}
        T2W_structural_proc=`dirname ${T2W_structural_proc}`/`basename ${T2W_structural_proc} .nii.gz`.nii
    fi
    
    echo "ABCD MESSAGE:  REORIENTING THE T2W IMAGE TO AXIAL"
    ReorientImage3D -i ${T2W_structural_proc}
    T2W_structural_proc_to_use=${up_proc_folder}/`basename ${T2W_structural_proc} .nii`_reoriented_weighted.nii

    echo "ABCD MESSAGE:  SMI MASKING THE T2W IMAGE"
    SmoothTransitionHCPStructural ${up_proc_folder}/`basename ${T2W_structural_proc} .nii`_reoriented.nii 100 ${T2W_structural_proc_to_use}
fi


mkdir ${up_proc_folder}/docker
mkdir ${up_proc_folder}/docker/INPUTS
mkdir ${up_proc_folder}/docker/OUTPUTS


gunzip -f ${eddy_outlier_free_data}
eddy_outlier_free_data=${up_proc_folder}/`basename ${eddy_outlier_free_data} .gz`
ExtractImage -i ${eddy_outlier_free_data} -v 0

echo "ABCD MESSAGE:  COPYING THE B0 IMAGE TO SYNB0-DISCO/INPUTS"
cp ${up_proc_folder}/`basename ${eddy_outlier_free_data} .nii`_V000.nii ${up_proc_folder}/docker/INPUTS/b0_lowres.nii
cp ${T1W_structural_proc_to_use} ${up_proc_folder}/docker/INPUTS/T1.nii 
gzip -f ${up_proc_folder}/docker/INPUTS/T1.nii 



echo "ABCD MESSAGE:  RIGIDLY REGISTERING THE B0 IMAGE TO T1W IMAGE"
antsRegistration -d 3 -o \[${up_proc_folder}/docker/INPUTS/trans,${up_proc_folder}/docker/INPUTS/b0.nii.gz\] -n Bspline\[3\] -r \[${up_proc_folder}/docker/INPUTS/T1.nii.gz,${up_proc_folder}/docker/INPUTS/b0_lowres.nii,0\] -m MI\[${up_proc_folder}/docker/INPUTS/T1.nii.gz,${up_proc_folder}/docker/INPUTS/b0_lowres.nii,1,50\] -t Rigid\[0.25\]  -f 6x4x2x1 -c 100x100x100x100 -s 3x2x1x0vox




echo "ABCD MESSAGE: SYNB0-DISCO"
#T1 to T2 conversion with SynB0
cd ${up_proc_folder}/docker
singularity run -e \
-B INPUTS/:/INPUTS \
-B OUTPUTS/:/OUTPUTS \
-B ${script_folder}/license.txt:/extra/freesurfer/license.txt \
${script_folder}/SYNB0-DISCO/synb0_latest.sif
cd ${cwd}


fslmaths ${up_proc_folder}/docker/OUTPUTS/b0_all ${up_proc_folder}/docker/OUTPUTS/b0_all_float -odt float
gunzip -f ${up_proc_folder}/docker/OUTPUTS/b0_all_float.nii.gz


echo "ABCD MESSAGE:  GENERATING THE SYNTHETIC T2W from STNB0-DISCO OUTPUT"
ExtractImage -i ${up_proc_folder}/docker/OUTPUTS/b0_all_float.nii -v 1
cp ${up_proc_folder}/docker/OUTPUTS/b0_all_float_V001.nii ${up_proc_folder}/T2W_synth_str.nii



# IMPORT DWIs TO TORTOISE

echo "ABCD MESSAGE:  IMPORTING THE UP AND DOWN DATA TO TORTOISE"
ImportNIFTI -i ${eddy_outlier_free_data} -p vertical -b ${up_bvals_proc} -v ${up_bvecs_proc}
ImportNIFTI -i ${down_nifti_proc} -p vertical -b ${down_bvals} -v ${down_bvecs}


up_list=${up_proc_folder}/`basename ${eddy_outlier_free_data} .nii`_proc/`basename ${eddy_outlier_free_data} .nii`.list
up_proc_list=${up_list}

down_list=${up_proc_folder}/`basename ${down_nifti_proc} .nii`_proc/`basename ${down_nifti_proc} .nii`.list
down_proc_list=${down_list}

if [ ${down_nvol} -ge 2 ]
then
    echo "ABCD MESSAGE:  RUNNING MOTION CORRECTION ON DOWN DATA"
    DIFFPREP -i ${down_list} --will_be_drbuddied 1 --gibbs_ringing_correction off 
    down_proc_list=${up_proc_folder}/`basename ${down_nifti_proc} .nii`_proc/`basename ${down_nifti_proc} .nii`_proc.list
fi


# NONLINEARITY FIELD CONVERSION,  SUSCEPTIBILITY DISTORTION CORRECTION, ALIGNMENT TO THE T1W IMAGE


if [ -z "${output_dir}"  ]  
then
   up_final_folder=${up_folder}/final_data
else
   up_final_folder=${output_dir}
fi

if [ ! -e "${up_final_folder}" ]
then
    mkdir -p ${up_final_folder}
fi

up_final_list=${up_final_folder}/`basename ${up_nifti_proc} .nii`_final.list
up_final_nii=${up_final_folder}/`basename ${up_nifti_proc} .nii`_final.nii

down_final_list=${up_final_folder}/`basename ${down_proc_list} .list`_DRBUDDI_final.list


up_data_parameter="--up_data ${up_proc_list}"
down_data_parameter="--down_data ${down_proc_list}"
structural_parameter="--structural ${up_proc_folder}/T2W_synth_str.nii"
if [ ${use_T2} -eq 1 ] && [ ! -z "${T2W_structural}"  ] 
then
    structural_parameter="${structural_parameter} --structural ${T2W_structural_proc_to_use}"
fi
output_parameter="--output ${up_final_list}"
bval_parameter="--DWI_bval_tensor_fitting 1200"

if [  -z "${output_res}"  ]   
then
    res_parameter=""

else
    if [[ "${output_res}" == *"x"*  ]]
    then
        res=(${output_res//x/ })
        res_parameter="--res ${res[0]} ${res[1]} ${res[2]}"       
    else
        res_parameter="--res ${output_res} ${output_res} ${output_res}"    
    fi      
fi


if [  -z "${output_size}"  ]   
then
    size_parameter=""

else
    if [[ "${output_size}" == *"x"*  ]]
    then
        size=(${output_size//x/ })
        size_parameter="--nvoxels ${size[0]} ${size[1]} ${size[2]}"          
    fi      
fi





gradnonlin_parameter="--gradnonlin_output_type grad_dev"
if [ ${isocenter} -eq 1 ]
then
    gradnonlin_parameter="${gradnonlin_parameter} --gradwarp_mode ge_${gw_mode}  "
else  
    gradnonlin_parameter="${gradnonlin_parameter} --gradwarp_mode ${gw_mode}  "  
fi


if [ ! -z "${gw_coeffs_file}"  ]
then
    gradnonlin_parameter="${gradnonlin_parameter} --gradnonlin_coeff_file ${gw_coeffs_file}"
else
    if [ ! -z "${gw_field_full}"  ]
    then
        gradnonlin_parameter="${gradnonlin_parameter} --grad_nonlin_distortion_field ${gw_field_full}"  
    else
	    if [ ! -z "${gw_field_L_file}"  ]
	    then
		export FREESURFER_HOME=${script_folder}/FREESURFER
		echo "ABCD MESSAGE:  CONVERTING THE FREESURFER FIELDS TO TORTOISE FIELDS"
		mri_convert ${gw_field_L_file}  ${up_proc_folder}/field_L.nii
		mri_convert ${gw_field_P_file}  ${up_proc_folder}/field_P.nii
		mri_convert ${gw_field_H_file}  ${up_proc_folder}/field_H.nii     
		
		ImageMath 3  ${up_proc_folder}/field_L.nii m ${up_proc_folder}/field_L.nii -1  
		ImageMath 3  ${up_proc_folder}/field_P.nii m ${up_proc_folder}/field_P.nii -1
		        
		CreateDisplacementField 3 0 ${up_proc_folder}/field_L.nii ${up_proc_folder}/field_P.nii ${up_proc_folder}/field_H.nii ${up_proc_folder}/field.nii
		
		InvertTransformation ${up_proc_folder}/field.nii
		ApplyTransformationToScalar ${up_proc_folder}/T2W_synth_str.nii ${up_proc_folder}/field_inv.nii ${up_proc_folder}/T2W_synth_str.nii ${up_proc_folder}/T2W_synth_str.nii
		        
		gradnonlin_parameter="${gradnonlin_parameter} --grad_nonlin_distortion_field ${up_proc_folder}/field.nii"               
	    fi
    fi
fi

echo "ABCD MESSAGE:  DRBUDDI"
DR_BUDDIC_withoutGUI_011921 ${up_data_parameter} ${down_data_parameter} ${structural_parameter} ${output_parameter} ${bval_parameter} ${gradnonlin_parameter} ${res_parameter} ${size_parameter}
cp ${up_final_folder}/structural.nii ${up_final_folder}/T2W_synth_structural.nii

# THIS next part is because we want untouched T1W. So we reoutput the DWIs in the T1W space.
structural_parameter2="--structural ${T1W_structural_proc_to_use}"
DR_BUDDIC_withoutGUI_011921 ${up_data_parameter} ${down_data_parameter} ${structural_parameter2} ${output_parameter} ${bval_parameter} ${gradnonlin_parameter} ${res_parameter} ${size_parameter} --step 3
cp ${up_final_folder}/structural.nii ${up_final_folder}/T1W_structural.nii


if [ ${output_orientation} !=  "LPS" ]
then 

    echo "ABCD MESSAGE:  REORIENTING THE IMAGES TO DESIRED ORIENTATION"
    ReorientImage3D -i ${up_final_folder}/T1W_structural.nii -d ${output_orientation} -o ${up_final_folder}/T1W_structural.nii
    ReorientImage -i ${up_final_list} -d ${output_orientation} -o ${up_final_list}
    ReorientImage -i ${down_final_list} -d ${output_orientation} -o ${down_final_list}
fi



echo "ABCD MESSAGE:  RIGIDLY REGISTERING THE T1W IMAGE TO THE B0 IMAGE"
antsRegistration -d 3 -o \[${up_final_folder}/rigidtrans_T2,${up_final_folder}/T2W_synth_structural.nii\] -n Bspline\[3\] -r \[${up_final_folder}/T1W_structural.nii,${up_final_folder}/T2W_synth_structural.nii,0\] -m MI\[${up_final_folder}/T1W_structural.nii,${up_final_folder}/T2W_synth_structural.nii,1,50\] -t Rigid\[0.25\]  -f 6x4x2x1 -c 100x100x100x100 -s 3x2x1x0vox 

if [ ! -z "${T2W_structural}"  ] 
then
    antsRegistration -d 3 -o \[${up_final_folder}/rigidtrans_T22,${up_final_folder}/T2W_structural.nii\] -n Bspline\[3\] -r \[${up_final_folder}/T1W_structural.nii,${T2W_structural},0\] -m MI\[${up_final_folder}/T1W_structural.nii,${T2W_structural},1,50\] -t Rigid\[0.25\]  -f 6x4x2x1 -c 100x100x100x100 -s 3x2x1x0vox 
fi



# FINAL TOUCHES

rm ${up_final_folder}/structural_used*.nii
rm ${up_final_folder}/blip_up*.nii*
rm ${up_final_folder}/blip_down*
rm ${up_final_folder}/b0_corrected_final.nii
#rm ${up_final_folder}/b0_to_str_rigid.txt
#rm ${up_final_folder}/deformation*




ExtractImage -i ${up_final_nii} -v 0
bet2 ${up_final_folder}/`basename ${up_final_nii} .nii`_V000.nii ${up_final_folder}/mask -m -f 0.3
EstimateTensorNLLS -i ${up_final_list} -m ${up_final_folder}/mask_mask.nii.gz -b 1200 -L ${up_final_folder}/`basename ${up_final_nii} .nii`_graddev.nii
ComputeDECMap -i ${up_final_folder}/`basename ${up_final_nii} .nii`_N1_DT.nii
ComputeFAMap  ${up_final_folder}/`basename ${up_final_nii} .nii`_N1_DT.nii
ComputeTRMap  ${up_final_folder}/`basename ${up_final_nii} .nii`_N1_DT.nii

TORTOISEBmatrixToFSLBVecs  ${up_final_folder}/`basename ${up_final_nii} .nii`.bmtxt


rm ${up_final_folder}/mask*



