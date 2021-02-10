#!/bin/bash
#SBATCH --job-name=sens_GP
#SBATCH --account=smith
#SBATCH -o /scicore/home/penny/burlyd00/smc_lai/JOB_OUT/sens_GP.out
#SBATCH --mem=200G
#SBATCH --qos=6hours
#SBATCH --cpus-per-task=4

#######################################
# script for sensitivity analysis
#
# created 17.06.2019
# monica.golumbeanu@unibas.ch
######################################
ml purge
ml R/3.6.0-foss-2018b

Rscript sens_GP_pppy_lowHL.R 
