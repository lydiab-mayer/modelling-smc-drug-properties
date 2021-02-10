#!/bin/bash
#
################################################
# Main TPP simulation workflow script 
#
# INPUT: 
#       SIM_FOLDER = simulation  folder with all necessary input files:
#			    param_tab.txt = file where each line is a parameter configuration to simulate
#			    scaffold.xml = xml file containing @parameter@ wildcards for varied values across simulations
#
# OUTPUT:
#	Scenario and simulation files along with postprocessing results.
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
SCENARIOS_FOLDER=$SIM_FOLDER"scenarios/"
LOG_FILE=$SIM_FOLDER"LOG.txt"

# Destination of the job output files
JOB_OUTPUTS=~/MMC/TPP/JOB_OUT/
    
# Create scenarios and run OpenMalaria simulations
echo "Step1: Creating scenarios and running OpenMalaria" > $LOG_FILE
cd 1_OM_basic_workflow/
bash OM_base_workflow.sh $SIM_FOLDER
echo "Number of simulation files:" >> $LOG_FILE
ls -l $OM_FOLDER | wc -l >> $LOG_FILE
cd -

# Resubmit unsuccessful jobs
cd 1_OM_basic_workflow/
echo "Resubmit unfinished jobs"
echo "Resubmit unfinished jobs" >> $LOG_FILE
bash OM_resubmit_workflow.sh $SIM_FOLDER
cd -
# TO DO: if number of outputs different than parameter table length, stop here
# if(ls -l $OM_FOLDER | wc -l /2 == wc -l $PARAM_TAB - 1 ) else write the two numbers in the LOG file
    
# Postprocess OpenMalaria simulation results
POSTPROCESSING_FOLDER=$SIM_FOLDER"postprocessing_4/"
echo "Step2: Postprocessing"
echo "Step2: Postprocessing" >> $LOG_FILE
cd 2_postprocessing/
bash postprocessing_workflow.sh $SIM_FOLDER 4
echo "Number of postprocessed files:" >> $LOG_FILE
ls -l $POSTPROCESSING_FOLDER | wc -l >> $LOG_FILE
cd -

# Postprocess OpenMalaria simulation results
POSTPROCESSING_FOLDER=$SIM_FOLDER"postprocessing_6/"
echo "Step2: Postprocessing"
echo "Step2: Postprocessing" >> $LOG_FILE
cd 2_postprocessing/
bash postprocessing_workflow.sh $SIM_FOLDER 6
echo "Number of postprocessed files:" >> $LOG_FILE
ls -l $POSTPROCESSING_FOLDER | wc -l >> $LOG_FILE
cd -

# TO DO: add condition to remove files
rm -r $OM_FOLDER
rm -r $SCENARIOS_FOLDER


