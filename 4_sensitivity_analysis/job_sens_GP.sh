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

GP_DIR=$1
PARAM_RANGES_FILE=$2
SENS_DEST_DIR=$3

# IMPORTANT: the number of files must equal to the task array length (index starts at 0)
gp_files=(${GP_DIR}*.RData)
#i=0
#echo ${split_files[$i]}


# Select scenario file in array
ID=$(expr ${SLURM_ARRAY_TASK_ID} - 1)
gp_file=${gp_files[$ID]}
echo "Postprocessing for $gp_file"

Rscript sens_GP_pppy_lowHL.R $gp_file $PARAM_RANGES_FILE $SENS_DEST_DIR
