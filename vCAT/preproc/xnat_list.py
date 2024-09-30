#!/usr/bin/env python
import os
import json
import argparse
import re
import shutil
import subprocess as sp
from datetime import datetime
import xnat

def json_load(filename): #quick function to load json from file in one line   
    with open(filename, 'r') as read_json:
        data = json.load(read_json)
    return data

def add_to_tar(tar_path, input_item): #adds path to tar file
    sp.Popen(['tar', '-czf', tar_path, input_item]).wait() 

def xget_file(credentials=None, project=None, #gets XNAT files
              _filter=None, dicom_dir=None): 
    if not _filter:
        _filter = r'[\w\w]'
    work_dir = os.getcwd()
    #produces dictionary of subject and session labels already grabbed from xnat,
    #and if json from previous run exists loads them for exclusion
    previous_subjs = {}
    os.makedirs(dicom_dir, exist_ok=True)
    subjs_json = dicom_dir + '/downloaded_subjects.json'
    credentials = json_load(credentials) #xnatpassfile
    previous_subjs = {}
    if os.path.isfile(subjs_json): #if not first retrieval
        previous_subjs = json_load(subjs_json) 
    #creates XNAT session and looks for new scans, ends session once completed
    session = xnat.connect(credentials['server'], #uses info from xnatpassfile to login
                           user=credentials['user'],
                           password=credentials['password'])
    project_data = session.projects[project] #selects my specific project
    for subject in project_data.subjects: #for each subject within my project
        subject_data = project_data.subjects[subject]
        if not re.search(_filter, subject_data.label): #does this do anything?!
            continue
        ##new subjects##
        if subject_data.label not in previous_subjs.keys(): #is the subject new?
            previous_subjs[subject_data.label] = []
        for exp in subject_data.experiments: #presumably iterates through sessions
            exp_data = subject_data.experiments[exp]
            if exp_data.label not in previous_subjs[subject_data.label]:
                previous_subjs[subject_data.label].append(exp_data.label)
                exp_data.download_dir(dicom_dir) #downloads the XNAT data
                os.chdir(dicom_dir) #change dir to reduce number of tar dirs
                add_to_tar(exp_data.label + '.tar.gz', exp_data.label + '/scans/') #adds directory to tar
                shutil.rmtree(exp_data.label) #removes original uncompressed files
    session.disconnect() 
    os.chdir(work_dir) #switch back to working dir
    with open(dicom_dir + 'pull_list.txt', 'a+') as rf: #appends date and time
        rf.write('\nRun completed on ' + str(datetime.now()))    
    with open(subjs_json, 'w') as dump_file: #dumps json with any new scans grabbed and tarred from XNAT
        json.dump(previous_subjs, dump_file, indent=4)

if __name__ == '__main__':
    #new set of command line arguments --config is the credentials file
    parser = argparse.ArgumentParser('Arguments required to pull files')
    parser.add_argument('-c', '--credentials', dest='credentials', required=True, help='Path to credentials file for logging in to XNAT')
    parser.add_argument('-p', '--project', dest='project', required=True, help='Name of project on XNAT')
    parser.add_argument('-d', dest='dicom_dir', required=True, help='location of main directory for storing dicoms')
    parser.add_argument('-f', '--filter', dest='filter')
    args = parser.parse_args()
    xget_file(credentials=args.credentials, project=args.project, dicom_dir=args.dicom_dir, _filter=args.filter)
