##############################
# Main script for running OpenMalaria simulations on the cluster. 
# 
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
#	bash OM_base_workflow.sh SIM_FOLDER
# 
#
# created 03.05.2020
#lydia.burgert@unibas.ch adapted from theresa.reiker@unibas.ch
#############################
setwd("~/M3TPP/analysis_workflow/1_OM_basic_workflow")
exp = "xxx"



GROUP = "/scicore/home/penny/GROUP/M3TPP/"

SIM_FOLDER=paste0(GROUP,exp,"/")

# Extract the number of lines in the parameter table
sink(paste0("submission_new_",exp,".sh"))
cat("#!/bin/bash\n")

  cat("PARAM_TABLE_FILE=",SIM_FOLDER,"param_tab_new.txt","\n\n", sep ="")
  cat("SCAFFOLD_FILE=",SIM_FOLDER,"scaffold.xml","\n", sep ="")
  cat("BASE_FOLDER=",SIM_FOLDER,"base_new/","\n", sep ="")
  cat("SCENARIOS_FOLDER=",SIM_FOLDER,"scenarios_new/","\n", sep ="")
  cat("OM_FOLDER=",SIM_FOLDER,"om/","\n", sep ="")
  
  cat("NUM=$(wc -l < $PARAM_TABLE_FILE)","\n", sep ="")
 
  cat(" echo \"Creating $NUM-1 scenarios ...\" ","\n", sep ="")
  
  cat("sbatch -W ../job_create_scenarios.sh $PARAM_TABLE_FILE $SCAFFOLD_FILE $BASE_FOLDER $SCENARIOS_FOLDER","\n\n", sep ="")
  
  cat("rm -r $BASE_FOLDER","\n", sep ="")
  cat(" echo \"Running OM simulations... \" ","\n", sep ="")
  
  cat("sbatch -W --array=1-$NUM run_OM.sh $SCENARIOS_FOLDER $OM_FOLDER","\n\n", sep ="")
  

  # Close the sink!
  sink()
  # either submit run to cluster here or from terminal 
 #submission <<- system(paste0("bash submission",sprintf("%02i",j),exp,".sh"), intern=TRUE)  
  
