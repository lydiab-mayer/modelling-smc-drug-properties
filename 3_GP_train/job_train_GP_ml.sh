#!/bin/bash
#SBATCH --job-name=train_GP_prev_red
#SBATCH --account=smith
#SBATCH -o /scicore/home/smith/laagmi01/smc_lai/train_GP.out
#SBATCH --mem=2G
#SBATCH --time=00:15:00
#SBATCH --qos=30min
#SBATCH --cpus-per-task=1

#######################################
# Script for training a Gaussian Process emulator using
# a database of OpenMalaria simulation results.
# INPUT:
#       INPUT_DIR = folder containing postprocessing results for each setting
#       DEST_DIR = folder where the trained models will be saved
#       PREDICTED = name of the output measure to be predicted, must match a
#                   column name in the postprocessing data frame
#       RANGES_FILE = file with parameter names and ranges
#
# OUTPUT:
#       The script creates in the specified DEST_DIR a .RData file corresponding
#       to each trained GP emulator.
#
# SYNTHAX: 
#       sbatch --array=1-$NUM job_train_GP.sh $INPUT_DIR $DEST_DIR $PREDICTED
#
# created 17.06.2019
# monica.golumbeanu@unibas.ch
######################################
ml purge
ml R/3.6.0-foss-2018b

INPUT_DIR=$1
DEST_DIR=$2
PREDICTED=$3
RANGES_FILE=$4

# IMPORTANT: the number of files must equal to the task array length (index starts at 0)
setting_postprocessing_results=(${INPUT_DIR}seeds_*.txt)

# Select scenario file in array
ID=$(expr ${SLURM_ARRAY_TASK_ID} - 1)
setting_postprocessing_result=${setting_postprocessing_results[$ID]}
echo "Postprocessing for $setting_postprocessing_result"

Rscript train_GP_ml.R $setting_postprocessing_result $DEST_DIR $PREDICTED $RANGES_FILE
