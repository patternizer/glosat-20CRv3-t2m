#!/bin/bash
SBATCH --partition=short-serial 
SBATCH --job-name=pctl_05_anomaly
SBATCH -o %j.out 
SBATCH -e %j.err 
SBATCH --time=05:00

module load jaspy
python calc_pctl_05_anomaly.py
