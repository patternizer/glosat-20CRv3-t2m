#!/bin/bash
SBATCH --partition=short-serial 
SBATCH --job-name=mean_anomaly
SBATCH -o %j.out 
SBATCH -e %j.err 
SBATCH --time=05:00

module load jaspy
python calc_mean_anomaly.py
