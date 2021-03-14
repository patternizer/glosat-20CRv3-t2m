#!/bin/bash
SBATCH --partition=short-serial 
SBATCH --job-name=pctl_25-20CRv3-t2m-all-years
SBATCH -o %j.out 
SBATCH -e %j.err 
SBATCH --time=05:00

module load jaspy
python calc_pctl_25_20CRv3_t2m_all_years.py

