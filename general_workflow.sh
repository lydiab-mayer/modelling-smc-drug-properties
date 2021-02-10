#!/bin/bash
#
################################################
# Main TPP workflow script 
#
# INPUT: 
#       SIM_FOLDER = simulation  folder with all necessary input files:
#			    param_tab.txt = file where each line is a parameter configuration to simulate
#			    scaffold.xml = xml file containing @parameter@ wildcards for varied values across simulations
#       FOLLOW_UP: index of the year when the follow up is done
#
# OUTPUT:
#	The script creates the following folders in the specified SIM_FOLDER:
#		base/ = folder containing replacement patterns and base xml files 
#               (will be deleted after scenario creation)
#		scenarios/ = folder with scenario files
#		om/ = folder with OpenMalaria simulations
#
# SYNTHAX: 
#	bash OM_base_workflow.sh SIM_FOLDER
#
# created 02.10.2019
# monica.golumbeanu@unibas.ch
###############################################
 
# Define variables
SIM_FOLDER=$1
OM_FOLDER=$SIM_FOLDER"om/"
LOG_FILE=$SIM_FOLDER"LOG.txt"
JOB_OUTPUTS=~/MMC/TPP/JOB_OUT/

# Create scenarios and run OpenMalaria simulations
echo "Step1: Creating scenarios and running OpenMalaria" > $LOG_FILE
cd 1_OM_basic_workflow/
bash OM_base_workflow.sh $SIM_FOLDER
echo "Number of simulation files:" >> $LOG_FILE
ls -l $OM_FOLDER | wc -l >> $LOG_FILE
cd -

# Clear job output files
rm -r $JOB_OUTPUTS
mkdir $JOB_OUTPUTS

# Resubmit unsuccessful jobs
cd 1_OM_basic_workflow/
echo "Resubmit unfinished jobs" >> $LOG_FILE
bash OM_resubmit_workflow.sh $SIM_FOLDER
cd -

# Clear job output files
rm -r $JOB_OUTPUTS
mkdir $JOB_OUTPUTS

# # Run analysis with followup year 4
# bash main_analysis $SIM_FOLDER 4
# 
# # Run analysis with followup year 6
# bash main_analysis $SIM_FOLDER 6

