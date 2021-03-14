#!/bin/bash
SBATCH --partition=short-serial 
SBATCH --job-name=20CRv3-t2m-all-years
SBATCH -o %j.out 
SBATCH -e %j.err 
SBATCH --time=05:00

#SBATCH --tasks-per-node=2
#SBATCH --nodes=4

module load jaspy
python run_20CRv3_t2m_all_years.py
