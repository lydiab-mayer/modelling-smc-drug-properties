#!/bin/bash
#
##############################
# Main script for running a set of optimization problems on a GP model. 
# INPUT:
#       SIM_FOLDER: simulation folder with a predefined structure
#       FOLLOW_UP: year when the output is evaluated (can be 4 or 6)
# OUTPUT:
#       The script saves the optimisation results for each GP model in a predefined folder 
# SYNTHAX: 
#       bash optim_workflow.sh SIM_FOLDER FOLLOW_UP
# 
#
# created 14.09.2019
# monica.golumbeanu@unibas.ch
#############################

SIM_FOLDER=$1

PARAM_RANGES_FILE=$SIM_FOLDER"param_ranges.RData"
GP_FOLDER1=$SIM_FOLDER"gp_5/trained/UL/"
GP_FOLDER2=$SIM_FOLDER"gp_5/trained/CI_high_HR/"
GP_AS_FOLDER=$GP_FOLDER1"as/"
OPT_SETUP_FILE=$SIM_FOLDER"gp_5/optimisation/non_inf/opt_setup_file_lhs.txt"
OPT_DEST_DIR=$SIM_FOLDER"gp_5/optimisation/non_inf/"

# create destination directory
mkdir -p $OPT_DEST_DIR

# Submit GP optimization analysis array job
gp_files1=(${GP_FOLDER1}*.RData)
NUM=${#gp_files1[@]}
echo $NUM

OPT_SETUP=$(wc -l < $OPT_SETUP_FILE)

for row_n in `seq 1 1 $OPT_SETUP`
do
    sbatch --array=1 job_optim_GP_NI.sh $GP_FOLDER1 $GP_FOLDER2 $PARAM_RANGES_FILE $OPT_DEST_DIR $OPT_SETUP_FILE $row_n
   # sbatch --array=1-$NUM job_optim_GP_NI.sh $GP_FOLDER1 $GP_FOLDER2 $PARAM_RANGES_FILE $OPT_DEST_DIR $OPT_SETUP_FILE $row_n
done
