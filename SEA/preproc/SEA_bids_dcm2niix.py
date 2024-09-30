#!/usr/bin/env python
import os
import os.path as op
import json
import shutil
import tarfile
import subprocess as sp
from glob import glob
import pydicom as dcm
import pandas as pd
import numpy as np
import shutil 

#### FUNCTIONS ####

def convert_nifti(target_folder, destination): #conversion function, removes superfluous directory tree
    os.makedirs(destination, exist_ok=True) #make new directory in permanent storage
    src_dcms = '{0}/resources/DICOM/files/'.format(target_folder) 
    cmd = "dcm2niix -o {0} {1}".format(destination, src_dcms) #convert and save to respective scan type folder
    sp.Popen(cmd, shell=True).wait() #subprocess executes conversion command, waits if necessary
    for i in glob(destination + "/*"): #grabs all files
        idx = i.index(".") #returns index of "." in current path/filename to identify beginning of extension
        shutil.move(i, destination + i[idx:]) #move converted NIFTI from scratch to data location with extension
    shutil.rmtree(destination) #remove old tree
    
###################
       
df = pd.read_csv('/home/data/madlab/Pruden_SEA/code/preproc/heuristics.csv') #open heuristics file as pandas dataframe  

for tar_file in glob('/home/data/madlab/Pruden_SEA/sourcedata/Pruden_SEA_*-S1.tar.gz'): #iterate through Ss tars
#for tar_file in glob('/home/data/madlab/Pruden_SEA/sourcedata/Pruden_SEA_1008EH-S1.tar.gz'): #used for debugging
    sub = os.path.basename(tar_file).split('_')[-1].split('-')[0][:4] #gets only Ss ID 
    pain_name = os.path.basename(tar_file).split('_')[-1].split('-')[0][:] #gets only Ss ID
       
    if os.path.exists('/home/data/madlab/Pruden_SEA/dset/sub-{0}/'.format(sub)): #check to see if Ss is new
        continue        
    os.makedirs('/scratch/madlab/bidsify_sea/SEA-{0}/'.format(sub), exist_ok=True) #creates subdirectories
    
    #unzipping file into new directories without deletion 
    sp.Popen('tar -xf {0} -C /scratch/madlab/bidsify_sea/SEA-{1}'.format(tar_file, sub), shell=True).wait()
    #defines current working directory for subject data
    curr_dir = '/scratch/madlab/bidsify_sea/SEA-{0}/Pruden_SEA_{1}-S1/scans/'.format(sub, pain_name)
       
    #grabs relevant directories (excluding setter, 32ch, and MPRG) 
    all_dirs = sorted([x for x in glob(curr_dir + "*") if not "setter" in x and not "32ch" in x and not "MPRG" in x])
    for x in ['T1', 't2', 'pd_tse']: #deals with multiple scans of the same type (uses last scan acquired)
        type_nums = []
        type_idxs = []
        for idx, scan_type in enumerate(all_dirs): #iterate through the filtered scans     
            if x in scan_type: #iterate through targeted scan types
                type_nums.append(int(all_dirs[idx].split('/')[-1].split('-')[0])) #grab the scan CIS assigned scan #
                type_idxs.append(idx) #grab the all_dir index of targeted scans
                if len(type_nums) == 2: #when a second scan of the target type is found                  
                    all_dirs.pop(type_idxs[type_nums.index(min(type_nums))]) #remove the first of that type from all_dirs   
       
    for i in ['anat', 'dwi', 'fmap']: #iterate through folder types      
        os.makedirs('/home/data/madlab/Pruden_SEA/dset/sub-{0}/{1}/'.format(sub, i), exist_ok=True) #create dir for each
    completed = []    
    for curr_scan in all_dirs: #iterate through filtered scans
        for row, scan_type in enumerate(df['init_name']): #iterate through scan types defined in heuristics file
            if scan_type in curr_scan: #if init name matches current scan type 
                if not scan_type in completed:
                    break #use row where break stops
        #print(curr_scan, '/home/data/madlab/Pruden_SEA/dset/sub-{0}/'.format(sub)+df['dest_name'][row].format(sub))
        #use custom function to run dcm2niix on current dicoms
        convert_nifti(curr_scan,'/home/data/madlab/Pruden_SEA/dset/sub-{0}/'.format(sub)+df['dest_name'][row].format(sub)) 
        completed.append(scan_type)
        print(sub, curr_scan) #prints current Ss identity failure
