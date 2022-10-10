#!/bin/bash                                                                                                   

#SBATCH --job-name=optimisation_LAI                   
#SBATCH --cpus-per-task=1                  
#SBATCH --mem-per-cpu=1G                   
#SBATCH --time=05:00:00        
#SBATCH --qos=6hours        
#SBATCH --output=/scicore/home/penny/laagmi01/smc_lai/JOB_OUT/optimise_coverage.out           
#SBATCH --error=/scicore/home/penny/laagmi01/smc_lai/JOB_OUT/optimise_coverage.er


gp_file=$1
ranges_file=$2
results_folder=$3


ml R/3.6.0-foss-2018b
Rscript optimise_coverage.R $gp_file $ranges_file $results_folder
