#!/bin/bash
SBATCH --partition=short-serial 
SBATCH --job-name=max_anomaly
SBATCH -o %j.out 
SBATCH -e %j.err 
SBATCH --time=05:00

module load jaspy
python calc_max_anomaly.py
