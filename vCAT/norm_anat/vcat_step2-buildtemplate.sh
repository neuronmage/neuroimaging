#!/usr/bin/env python
#SBATCH --nodelist="IB_40C_512G"
#SBATCH -n 20

buildtemplateparallel.sh -d 3 -n 1 -g 0.200000 -i 4 -m 100x70x50x20 -o vcat_antsTMPL_ -c 2 -j 20 -r 1 -s CC -t GR *.nii.gz

