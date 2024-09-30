
#!/bin/bash -l
source ~/.bashrc;
sea_env;
#sbatch --qos qos_download
#sbatch --partition download



python /home/data/madlab/Pruden_SEA/code/preproc/xnat_list.py --credentials /home/arenf001/xnatpassfile -p Pruden_SEA -f Pruden_SEA -d /home/data/madlab/Pruden_SEA/sourcedata/
python /home/data/madlab/Pruden_SEA/code/preproc/SEA_bids_dcm2niix.py
