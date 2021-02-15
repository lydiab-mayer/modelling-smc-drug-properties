#!/bin/bash
#SBATCH --job-name=as_GP
#SBATCH --account=smith
#SBATCH -o /scicore/home/smith/burlyd00/smc_lai/JOB_OUT/as_gp.out
#SBATCH --mem=5G
#SBATCH --qos=6hours
#SBATCH --cpus-per-task=1

#######################################
# Script for training a Gaussian Process emulator using
# a database of OpenMalaria simulation results.
# INPUT:
#       INPUT_DIR = folder containing postprocessing results for each setting
#       DEST_DIR = folder where the trained models will be saved
#       RANGES_FILE = file with parameter names and ranges
#       FOLLOW_UP = integer representing the survey index to consider for 
#                   evaluating intervention impact
#       PREDICTED = name of the output measure to be predicted, must match a
#                   column name in the postprocessing data frame
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

SCAFFOLD_XML=$1
GP_TRAINING_DIR=$2 
AS_DIR=$3
RANGES_FILE=$4
FOLLOW_UP=$5
PREDICTED=$6


# IMPORTANT: the number of files must equal to the task array length (index starts at 0)
GP_models=(${GP_TRAINING_DIR}*.RData)

# Select GP model file in array
ID=$(expr ${SLURM_ARRAY_TASK_ID} - 1)
GP_model_file=${GP_models[$ID]}

# Create the folder for adaptive sampling resources (scenarios and simulations)
AS_RUN_DIR=$AS_DIR"as_model_"$ID"/"
mkdir -p $AS_RUN_DIR
scp $SCAFFOLD_XML $AS_RUN_DIR

# Rscript adaptive_design_horizon_om.R $AS_RUN_DIR $GP_model_file $RANGES_FILE $FOLLOW_UP $PREDICTED
