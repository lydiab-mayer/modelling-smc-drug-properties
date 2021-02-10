#!/bin/bash                                                                                                   

#SBATCH --job-name=train_gp                  
#SBATCH --cpus-per-task=1                  
#SBATCH --mem-per-cpu=1G
#SBATCH --time=00:15:00                   
#SBATCH --qos=30min    
#SBATCH --output=/scicore/home/penny/burlyd00/smc_lai/JOB_OUT/trainGP.o%j            
#SBATCH --error=/scicore/home/penny/burlyd00/smc_lai/JOB_OUT/trainGP.e%j  


SEASONALITY=$1
DECAY=$2
ACESS=$3
AGE=$4
EIR=$5


ml R/3.6.0-foss-2018b
Rscript train_gp_implementation_setting_pppy.R $SEASONALITY $DECAY $ACESS $AGE $EIR
