#!/bin/bash
SBATCH --partition=short-serial 
SBATCH --job-name=pctl_10_anomaly
SBATCH -o %j.out 
SBATCH -e %j.err 
SBATCH --time=05:00

module load jaspy
python calc_pctl_10_anomaly.py
