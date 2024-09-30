#!/usr/bin/env python
import os
import os.path as op
import json
import shutil
import tarfile
import subprocess as sp
from glob import glob
import dcmstack as ds
import dicom as dcm
import pandas as pd
import numpy as np
import shutil 

def convert_nifti(target_folder, destination):
    os.makedirs(destination, exist_ok=True)
    src_dcms = '{0}/resources/DICOM/files/'.format(target_folder)
    cmd = "dcm2niix -o {0} {1}".format(destination, src_dcms)
    sp.Popen(cmd, shell=True).wait()
    for i in glob(destination + "/*"):
        idx = i.index(".")
        shutil.move(i, destination + i[idx:])
    shutil.rmtree(destination)
    
def save_json(filename, data): 
    with open(filename, 'w') as fp:
        json.dump(data, fp, sort_keys=True, indent=4)
        
df = pd.read_csv('heuristics.csv')        
        
with open('Reeb_EBC.json', 'r') as _:
    config = json.load(_)
    
#used to extract the desired header info
info_keys = ["ConversionSoftware", "ConversionSoftwareVersion", "EchoTime", "FlipAngle", 
             "MRAcquisitionType", "MagneticFieldStrength", "Manufacturer", "ManufacturersModelName",
             "Modality", "PixelBandwidth", "ProtocolName", "RepetitionTime", "SAR", "ScanningSequence",
             "SequenceName", "SequenceVariant", "SeriesDescription", "SeriesNumber",  "SliceThickness"]   

for tar_file in glob('/home/data/madlab/Reeb_EBC/sourcedata/Reeb_EBC-000-*D-S1.tar.gz'):
    check = os.path.basename(tar_file).split('_')[-1].split('-')[2]
    sub = check[:-1] #gets only subject ID
    if os.path.exists('/home/data/madlab/Reeb_EBC/dset/sub-{0}/'.format(sub)): #check to see if this Ss is new
        continue       
    os.makedirs('/scratch/madlab/bidsify_ebc/EBC-{0}/'.format(sub), exist_ok=True) #creates directories
    #unzipping the file without deletion into new directories
    sp.Popen('tar -xf {0} -C /scratch/madlab/bidsify_ebc/EBC-{1}'.format(tar_file, sub), shell=True).wait()
    long_dir='/scratch/madlab/bidsify_ebc/EBC-{0}/home/data/madlab/Reeb_EBC/sourcedata/Reeb_EBC-000-{0}D-S1/scans/'.format(sub)
    short_dir = '/scratch/madlab/bidsify_ebc/EBC-{0}/Reeb_EBC-000-{0}D-S1/scans/'.format(sub)
    if os.path.exists(long_dir):
        curr_dir = long_dir #for old downloadeds
    else:
        curr_dir = short_dir #for new downloads
    all_dirs = [x for x in glob(curr_dir + "*") if not "setter" in x and not "32ch" in x and not "MPRG" in x]
    corr_scan = {} #contains only corrected scans 
    corr_idx = {} #contains indices of correct scans within all_dirs
    for idx, curr_scan in enumerate(all_dirs): #dealing with the inconsistent numbering of scans
        i = int(curr_scan.split("/")[-1].split("-")[0]) #get scan number
        for row, scan in enumerate(df['init_name']): #scan names in heuristics file
            if scan in curr_scan: #if one of the scans we care about
                if not scan in corr_scan or corr_scan[scan] < i: 
                    corr_idx[scan] = idx
                    corr_scan[scan] = i 
    temp = []
    for key in corr_scan:
        temp.append(all_dirs[corr_idx[key]])
    all_dirs = temp
        
    for i in ['anat', 'dwi', 'fmap']:       
        os.makedirs('/home/data/madlab/Reeb_EBC/dset/sub-{0}/{1}/'.format(sub, i), exist_ok=True)
    for curr_scan in all_dirs:
        for row, scan in enumerate(df['init_name']): #heuristics file columns
            if scan in curr_scan:
                break 
        convert_nifti(curr_scan,'/home/data/madlab/Reeb_EBC/dset/sub-{0}/'.format(sub)+df['dest_name'][row].format(sub)) 
        print(sub)
        dest_dir = '/home/data/madlab/Reeb_EBC/dset/sub-{0}/'.format(sub)
        saving_dict = {}
        for d in glob('{0}/*/DICOM/files/*.dcm'.format(curr_scan))[:1]:
            mw = ds.wrapper_from_data(dcm.read_file(d, force=True))
            for key in info_keys:
                if key in mw.dcm_data:
                    saving_dict[key] = getattr(mw.dcm_data,key)
                else:
                    if "Version" in key:
                        saving_dict[key] = "1.0.20170314"
                    else:
                        saving_dict[key] = "dcm2niix"
        save_json(dest_dir + df['dest_name'][row].format(sub)+'.json', saving_dict)       