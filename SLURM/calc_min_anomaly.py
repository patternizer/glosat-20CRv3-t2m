#!/usr/bin/env python
    
#------------------------------------------------------------------------------
# PROGRAM: calc_min_anomaly.py
#------------------------------------------------------------------------------
# Version 0.2
# 12 March, 2021
# Michael Taylor
# https://patternizer.github.io
# patternizer AT gmail DOT com
# michael DOT a DOT taylor AT uea DOT ac DOT uk
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# IMPORT PYTHON LIBRARIES
#------------------------------------------------------------------------------
import numpy as np
import os
import stat
import subprocess

#------------------------------------------------------------------------------
# SLURM SHELL SCRIPT GENERATING FUNCTION
#------------------------------------------------------------------------------

def make_shell_command(year):
         
    job_id = 'SBATCH --job-name=min_anomaly.{0}\n'.format(year)      
#   job_normals = 'ncea ensemble-min/{1961..1990}_ensemble_min.nc ensemble_min_normals.nc\n'
    job_nco = 'ncdiff ensemble-min/{0}_ensemble_min.nc ensemble_mean_normals.nc {0}_ensemble_min_anomaly.nc\n'.format(year)
    job_file = 'run.{0}.min_anomaly.sh'.format(year)
    with open(job_file,'w') as fp:
        fp.write('#!/bin/bash\n')
        fp.write('SBATCH --partition=short-serial\n')
        fp.write(job_id)
        fp.write('SBATCH -o %j.out\n')
        fp.write('SBATCH -e %j.err\n')
        fp.write('SBATCH --time=05:00\n')
        fp.write('module load jaspy\n')
#       fp.write(job_normals)
        fp.write(job_nco)

    # Make script executable

    os.chmod(job_file,stat.S_IRUSR | stat.S_IWUSR | stat.S_IXUSR)

    # Submit script to SLURM

    job = ['sbatch',job_file]
    subprocess.call(job)

if __name__ == "__main__":
    
    year1 = 1806
    year2 = 2021
    for year in range(year1,year2):
        file_out = '{0}_ensemble_min_anomaly.nc'.format(year)
        if not os.path.exists(file_out):
            make_shell_command(year)

