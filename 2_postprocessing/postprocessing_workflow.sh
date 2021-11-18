#!/bin/bash
#SBATCH --mem=2G
#SBATCH --qos=6hours
#SBATCH --cpus-per-task=1
#
##############################
# Main script for postprocessing OpenMalaria simulations on the cluster. 
# INPUT:
#       SIM_FOLDER = folder with simulation results
#       PARAM_TAB = file with the simulation parameters
#       OUTPUT_FOLDER = folder with postprocessing results
#       FOLLOW_UP = integer representing the survey index to consider for 
#                   evaluating intervention impact
# OUTPUT:
#       The script creates the following folders in the specified OUTPUT_FOLDER:
#               split/ = folder containing scenario parameters for each setting
#               processed/ = folder containing processing results
#
# SYNTHAX: 
#       bash postprocessing_workflow.sh SIM_FOLDER FOLLOW_UP
# 
#
# created 14.09.2019
# monica.golumbeanu@unibas.ch
#############################
ml purge
ml R/3.6.0-foss-2018b

SIM_FOLDER=$1
DATE=$2
FMONTH=$3
MONTHS=$4
YEAR_COUNTERFACTUAL=$5
YEAR_INTERVENTION=$6
MIN_INT=$7

PARAM_TAB=$SIM_FOLDER"param_tab.txt"
PARAM_CAT=$SIM_FOLDER"param_ranges_cat.RData"
OM_FOLDER=$SIM_FOLDER"om/"
OUTPUT_FOLDER=$SIM_FOLDER"postprocessing/"
ERROR_FOLDER=$OUTPUT_FOLDER"err/"
SPLIT_FOLDER=$OUTPUT_FOLDER"split/"

# create the necessary folders
mkdir -p $OUTPUT_FOLDER
mkdir -p $SPLIT_FOLDER
mkdir -p $ERROR_FOLDER

# split the parameter table by setting
Rscript ../../../analysisworkflow/2_postprocessing/split_param_table.R $PARAM_TAB $SPLIT_FOLDER $PARAM_CAT

echo $PARAM_TAB
echo $SPLIT_FOLDER
echo $PARAM_CAT

# Submit postprocessing array job
split_files=(${SPLIT_FOLDER}*.txt)
NUM=${#split_files[@]}
sbatch -W --array=1-$NUM job_postprocessing.sh $SPLIT_FOLDER $OM_FOLDER $DATE $FMONTH $MONTHS $YEAR_COUNTERFACTUAL $YEAR_INTERVENTION $MIN_INT
