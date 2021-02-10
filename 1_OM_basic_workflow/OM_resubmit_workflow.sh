#!/bin/bash
#
##############################
# Main script for resubmitting OpenMalaria simulations on the cluster. 
# INPUT:
#	SIM_FOLDER = simulation  folder with all necessary input files:
#			param_tab.txt = file where each line is a parameter configuration to simulate
#			scaffold.xml = xml file containing @parameter@ wildcards for varied values across simulations
# OUTPUT:
#	The script creates the following folders in the specified SIM_FOLDER:
#		base/ = folder containing replacement patterns and base xml files
#		scenarios/ = folder with scenario files
#		om/ = folder with OpenMalaria simulations
#
# SYNTHAX: 
#	bash OM_resubmit_workflow.sh SIM_FOLDER
# 
#
# created 14.09.2019
# monica.golumbeanu@unibas.ch
#############################

SIM_FOLDER=$1
PARAM_TABLE_FILE=$SIM_FOLDER"param_tab.txt"
SCENARIOS_FOLDER=$SIM_FOLDER"scenarios/"
OM_FOLDER=$SIM_FOLDER"om/"
LOG_FILE=$SIM_FOLDER"resubmit.txt"

# Extract the number of lines in the parameter table
NUM=$(wc -l < $PARAM_TABLE_FILE)
echo "" > $LOG_FILE
# Run OpenMalaria for each scenario
echo "Running OpenMalaria simulations ..."
sbatch -W --array=1-$NUM resubmit_OM.sh $OM_FOLDER $SCENARIOS_FOLDER
