#!/bin/bash
#sbatch --qos qos_download
#sbatch --partition download
source ~/.bashrc;
vcat_env;

python /home/data/madlab/Mattfeld_vCAT/code/preproc/xnat_list.py --credentials /home/arenf001/xnatpassfile -p Mattfeld_REVL -f Mattfeld_REVL-000-vCAT -d /home/data/madlab/Mattfeld_vCAT/sourcedata/
python /home/data/madlab/Mattfeld_vCAT/code/preproc/bids_dcm2niix.py
