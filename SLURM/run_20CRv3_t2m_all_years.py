#!/usr/bin/env python
    
#------------------------------------------------------------------------------
# PROGRAM: run_20CRv3_t2m_all_years.py
#------------------------------------------------------------------------------
# Version 0.3
# 5 March, 2021
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
         
    job_id = 'SBATCH --job-name=20CRv3-t2m.{0}\n'.format(year)        
    job_str = 'wget https://portal.nersc.gov/archive/home/projects/incite11/www/20C_Reanalysis_version_3/everymember_anal_netcdf/mnmean/TMP2m/TMP2m_{0}_mnmean.tar\n'.format(year) 
    job_file = 'run.{0}.sh'.format(year)
    with open(job_file,'w') as fp:
        fp.write('#!/bin/bash\n')
        fp.write('SBATCH --partition=short-serial\n')
        fp.write(job_id)
        fp.write('SBATCH -o %j.out\n')
        fp.write('SBATCH -e %j.err\n') 
        fp.write('SBATCH --time=05:00\n')
        #SBATCH --tasks-per-node=2
        #SBATCH --nodes=4
        fp.write(job_str)    

    # Make script executable

    os.chmod(job_file,stat.S_IRUSR | stat.S_IWUSR | stat.S_IXUSR)

    # Submit script to SLURM

    job = ['sbatch',job_file]
    subprocess.call(job)

if __name__ == "__main__":
    
    year1 = 1806
    year2 = 2021
    for year in range(year1,year2):
        file_out = 'TMP2m_{0}_mnmean.tar'.format(year)                       
        if not os.path.exists(file_out):
            make_shell_command(year)       

