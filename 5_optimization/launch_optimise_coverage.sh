#!/bin/bash                                                                                                   

#SBATCH --job-name=optimisation_LAI                   
#SBATCH --cpus-per-task=1                  
#SBATCH --mem-per-cpu=1G                   
#SBATCH --time=05:00:00        
#SBATCH --qos=6hours        
#SBATCH --output=/scicore/home/penny/laagmi01/smc_lai/JOB_OUT/optimise_coverage.out           
#SBATCH --error=/scicore/home/penny/laagmi01/smc_lai/JOB_OUT/optimise_coverage.er



ID=$(expr ${SLURM_ARRAY_TASK_ID})
echo "$ID"
ml R/3.6.0-foss-2018b
Rscript optimise_coverage.R $ID
