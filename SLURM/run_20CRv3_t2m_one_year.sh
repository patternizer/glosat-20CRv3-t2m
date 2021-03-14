#!/bin/bash
SBATCH --partition=short-serial 
SBATCH --job-name=20CRv3-t2m-1806
SBATCH -o %j.out 
SBATCH -e %j.err 
SBATCH --time=05:00

wget https://portal.nersc.gov/archive/home/projects/incite11/www/20C_Reanalysis_version_3/everymember_anal_netcdf/mnmean/TMP2m/TMP2m_1806_mnmean.tar

