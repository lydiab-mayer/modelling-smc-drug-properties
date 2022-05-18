#!/bin/bash
#SBATCH --job-name=adaptive_sampling
#SBATCH --account=penny
#SBATCH -o /scicore/home/penny/brauna0000/M3TPP/Experiments/iTPP3_bloodstage_mixed/JOB_OUT/as_submission.out
#SBATCH --mem=2G
#SBATCH --qos=6hours
#SBATCH --cpus-per-task=1
#
##############################
# Main script for refining a Gaussian process emulator with adaptive sampling
# INPUT:
#       SIM_FOLDER: simulation folder containing all the input files as well as
#                   files and folders created by the various steps of the workflow.
#       PREDICTED = name of the output measure to be predicted, must match a
#                   column name in the postprocessing data frame
#       SCALE = TRUE/FALSE indicator for whether inputs have/should be scaled to c(0, 1)
#
# OUTPUT:
#
# SYNTAX: 
#       bash GP_as_workflow.sh SIM_FOLDER
#       bash GP_as_workflow.sh ~/MMC/TPP/simulations/MAB_once_3years_avg_prev/scaffold.xml ~/MMC/TPP/simulations/MAB_once_3years_avg_prev/gp_4/trained/ ~/MMC/TPP/simulations/MAB_once_3years_avg_prev/gp_4/as/ ~/MMC/TPP/simulations/MAB_once_3years_avg_prev/param_ranges.RData 4 prev_red
# 
#
# created 02.10.2019, adapted 18.05.2022
# monica.golumbeanu@unibas.ch
# lydia.braunack-mayer@swisstph.ch
#############################

SIM_FOLDER=$1
PREDICTED=$2
SCALE=$3

SCAFFOLD_FILE=$SIM_FOLDER"scaffold.xml"
GP_FOLDER=$SIM_FOLDER"gp/"$PREDICTED"/"
GP_TRAINING_FOLDER=$GP_FOLDER"trained/"
GP_AS_FOLDER=$GP_FOLDER"as/"

# Debug
echo $SIM_FOLDER
echo $PREDICTED
echo $SCALE
echo $SCAFFOLD_FILE
echo $GP_FOLDER
echo $GP_TRAINING_FOLDER
echo $GP_AS_FOLDER

mkdir -p $GP_AS_FOLDER

# Submit as array job
GP_models=(${GP_TRAINING_FOLDER}*.RData)
NUM=${#GP_models[@]}

for (( ID=0; ID<$NUM; ID++ ))
do  
    $GP_MODEL_FILE=${GP_models[$ID]}
    echo "Adaptive sampling for model $GP_MODEL_FILE"
    AS_RUN_DIR=$GP_AS_FOLDER"as_model_"$ID"/"
    mkdir -p $AS_RUN_DIR
    scp $SCAFFOLD_FILE $AS_RUN_DIR
    Rscript adaptive_design_om.R $AS_RUN_DIR $GP_MODEL_FILE $SCALE
done

# Older version where this was submitting an array job:
# sbatch -W --array=1-$NUM job_as_GP.sh $SCAFFOLD_XML $GP_TRAINING_DIR $AS_DIR $RANGES_FILE $FOLLOW_UP $PREDICTED


