#!/usr/bin/env python
    
#------------------------------------------------------------------------------
# PROGRAM: calc_anomaly_20CRv3_t2m_all_years.py
#------------------------------------------------------------------------------
# Version 0.1
# 8 March, 2021
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

def make_shell_command(year,member):
         
    job_id = 'SBATCH --job-name=20CRv3-t2m.{0}.{1:03d}\n'.format(year,member)        
    job_tar = 'tar xvf ensemble-tar/TMP2m_{0}_mnmean.tar\n'.format(year) 
    job_nco = 'ncdiff {0}/TMP2m.{0}.mnmean_mem{1:03d}.nc ensemble-mean/{0}_ensemble_mean.nc {0}_{1:03d}_ensemble_anomaly.nc\n'.format(year,member)
    job_file = 'run.{0}.{1:03d}.sh'.format(year,member)
    with open(job_file,'w') as fp:
        fp.write('#!/bin/bash\n')
        fp.write('SBATCH --partition=short-serial\n')
        fp.write(job_id)
        fp.write('SBATCH -o %j.out\n')
        fp.write('SBATCH -e %j.err\n')
        fp.write('SBATCH --time=05:00\n')
        fp.write('module load jaspy\n')
        fp.write(job_tar)
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
        for member in range(1,81):
            file_out = '{0}_{1:03d}_ensemble_anomaly.nc'.format(year,member)
            if not os.path.exists(file_out):
                make_shell_command(year,member)

