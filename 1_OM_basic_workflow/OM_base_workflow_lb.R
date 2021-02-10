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
setwd("~/smc_lai/analysis_workflow/1_OM_basic_workflow")
exp = "E3_SMCSMCdisc"
user= "burlyd00"



GROUP = "/scicore/home/penny/GROUP/smc_lai/"

SIM_FOLDER=paste0(GROUP,exp,"/")

# Extract the number of lines in the parameter table
no.commands=as.numeric(system(paste0("wc -l < /scicore/home/penny/GROUP/smc_lai/",exp, "/param_tab.txt"), intern = TRUE))

  no.bats = no.commands %/% 37500
  if(no.commands %% 37500 >0){no.bats = no.bats+1}
  for(j in 0:(no.bats-1)){ #"split" counting automatically starts at 0, so myst start counting from 0 here
    sink(paste0("/scicore/home/penny/",user,"/smc_lai/analysis_workflow/1_OM_basic_workflow/",exp,"/submission",sprintf("%02i",j),exp,".sh"))
    cat("#!/bin/bash\n")

  cat("PARAM_TABLE_FILE=",SIM_FOLDER,"param_tab_",j,".txt","\n\n", sep ="")
  cat("SCAFFOLD_FILE=",SIM_FOLDER,"scaffold.xml","\n", sep ="")
  cat("BASE_FOLDER=",SIM_FOLDER,"base_",j,"/","\n", sep ="")
  cat("SCENARIOS_FOLDER=",SIM_FOLDER,"scenarios_",j,"/","\n", sep ="")
  cat("OM_FOLDER=",SIM_FOLDER,"om/","\n", sep ="")
  
  cat("NUM=$(wc -l < $PARAM_TABLE_FILE)","\n", sep ="")
 
  cat(" echo \"Creating $NUM-1 scenarios ...\" ","\n", sep ="")
  
  cat("sbatch -W ../job_create_scenarios.sh $PARAM_TABLE_FILE $SCAFFOLD_FILE $BASE_FOLDER $SCENARIOS_FOLDER","\n\n", sep ="")
  
  cat(" echo \"Running OM simulations... \" ","\n", sep ="")
  
  cat("sbatch -W --array=1-$NUM ../run_OM.sh $SCENARIOS_FOLDER $OM_FOLDER","\n\n", sep ="")
  

  # Close the sink!
  sink()
  # either submit run to cluster here or from terminal 
 #submission <<- system(paste0("bash submission",, intern=TRUE)  
  } 

  exp = "E4_SMCLAIdisc"
  user= "burlyd00"
  

  
  GROUP = "/scicore/home/penny/GROUP/smc_lai/"
  
  SIM_FOLDER=paste0(GROUP,exp,"/")
  
  # Extract the number of lines in the parameter table
  no.commands=as.numeric(system(paste0("wc -l < /scicore/home/penny/GROUP/smc_lai/",exp, "/param_tab.txt"), intern = TRUE))
  
  no.bats = no.commands %/% 37500
  if(no.commands %% 37500 >0){no.bats = no.bats+1}
  for(j in 0:(no.bats-1)){ #"split" counting automatically starts at 0, so myst start counting from 0 here
    sink(paste0("/scicore/home/penny/",user,"/smc_lai/analysis_workflow/1_OM_basic_workflow/",exp,"/submission",sprintf("%02i",j),exp,".sh"))
    cat("#!/bin/bash\n")
    
    cat("PARAM_TABLE_FILE=",SIM_FOLDER,"param_tab_",j,".txt","\n\n", sep ="")
    cat("SCAFFOLD_FILE=",SIM_FOLDER,"scaffold.xml","\n", sep ="")
    cat("BASE_FOLDER=",SIM_FOLDER,"base_",j,"/","\n", sep ="")
    cat("SCENARIOS_FOLDER=",SIM_FOLDER,"scenarios_",j,"/","\n", sep ="")
    cat("OM_FOLDER=",SIM_FOLDER,"om/","\n", sep ="")
    
    cat("NUM=$(wc -l < $PARAM_TABLE_FILE)","\n", sep ="")
    
    cat(" echo \"Creating $NUM-1 scenarios ...\" ","\n", sep ="")
    
    cat("sbatch -W ../job_create_scenarios.sh $PARAM_TABLE_FILE $SCAFFOLD_FILE $BASE_FOLDER $SCENARIOS_FOLDER","\n\n", sep ="")
    
    cat(" echo \"Running OM simulations... \" ","\n", sep ="")
    
    cat("sbatch -W --array=1-$NUM ../run_OM.sh $SCENARIOS_FOLDER $OM_FOLDER","\n\n", sep ="")
    
    
    # Close the sink!
    sink()
    # either submit run to cluster here or from terminal 
    #submission <<- system(paste0("bash submission",, intern=TRUE)  
  } 
  
  

  
  exp = "E5_2_CT_LAI"
  user= "burlyd00"
  
  
  GROUP = "/scicore/home/penny/GROUP/smc_lai/"
  
  SIM_FOLDER=paste0(GROUP,exp,"/")
  
  # Extract the number of lines in the parameter table
  no.commands=as.numeric(system(paste0("wc -l < /scicore/home/penny/GROUP/smc_lai/",exp, "/param_tab.txt"), intern = TRUE))
  
  no.bats = no.commands %/% 37500
  if(no.commands %% 37500 >0){no.bats = no.bats+1}
  for(j in 0:(no.bats-1)){ #"split" counting automatically starts at 0, so myst start counting from 0 here
    sink(paste0("/scicore/home/penny/",user,"/smc_lai/analysis_workflow/1_OM_basic_workflow/",exp,"/submission",sprintf("%02i",j),exp,".sh"))
    cat("#!/bin/bash\n")
    
    cat("PARAM_TABLE_FILE=",SIM_FOLDER,"param_tab_",j,".txt","\n\n", sep ="")
    cat("SCAFFOLD_FILE=",SIM_FOLDER,"scaffold.xml","\n", sep ="")
    cat("BASE_FOLDER=",SIM_FOLDER,"base_",j,"/","\n", sep ="")
    cat("SCENARIOS_FOLDER=",SIM_FOLDER,"scenarios_",j,"/","\n", sep ="")
    cat("OM_FOLDER=",SIM_FOLDER,"om/","\n", sep ="")
    
    cat("NUM=$(wc -l < $PARAM_TABLE_FILE)","\n", sep ="")
    
    cat(" echo \"Creating $NUM-1 scenarios ...\" ","\n", sep ="")
    
    cat("sbatch -W ../job_create_scenarios.sh $PARAM_TABLE_FILE $SCAFFOLD_FILE $BASE_FOLDER $SCENARIOS_FOLDER","\n\n", sep ="")
    
    cat(" echo \"Running OM simulations... \" ","\n", sep ="")
    
    cat("sbatch -W --array=1-$NUM ../run_OM.sh $SCENARIOS_FOLDER $OM_FOLDER","\n\n", sep ="")
    
    
    # Close the sink!
    sink()
    # either submit run to cluster here or from terminal 
    #submission <<- system(paste0("bash submission",, intern=TRUE)  
  } 
  
  
  exp = "E5_1_CT_SMC"
  user= "burlyd00"
  
  
  GROUP = "/scicore/home/penny/GROUP/smc_lai/"
  
  SIM_FOLDER=paste0(GROUP,exp,"/")
  
  # Extract the number of lines in the parameter table
  no.commands=as.numeric(system(paste0("wc -l < /scicore/home/penny/GROUP/smc_lai/",exp, "/param_tab.txt"), intern = TRUE))
  
  no.bats = no.commands %/% 37500
  if(no.commands %% 37500 >0){no.bats = no.bats+1}
  for(j in 0:(no.bats-1)){ #"split" counting automatically starts at 0, so myst start counting from 0 here
    sink(paste0("/scicore/home/penny/",user,"/smc_lai/analysis_workflow/1_OM_basic_workflow/",exp,"/submission",sprintf("%02i",j),exp,".sh"))
    cat("#!/bin/bash\n")
    
    cat("PARAM_TABLE_FILE=",SIM_FOLDER,"param_tab_",j,".txt","\n\n", sep ="")
    cat("SCAFFOLD_FILE=",SIM_FOLDER,"scaffold.xml","\n", sep ="")
    cat("BASE_FOLDER=",SIM_FOLDER,"base_",j,"/","\n", sep ="")
    cat("SCENARIOS_FOLDER=",SIM_FOLDER,"scenarios_",j,"/","\n", sep ="")
    cat("OM_FOLDER=",SIM_FOLDER,"om/","\n", sep ="")
    
    cat("NUM=$(wc -l < $PARAM_TABLE_FILE)","\n", sep ="")
    
    cat(" echo \"Creating $NUM-1 scenarios ...\" ","\n", sep ="")
    
    cat("sbatch -W ../job_create_scenarios.sh $PARAM_TABLE_FILE $SCAFFOLD_FILE $BASE_FOLDER $SCENARIOS_FOLDER","\n\n", sep ="")
    
    cat(" echo \"Running OM simulations... \" ","\n", sep ="")
    
    cat("sbatch -W --array=1-$NUM ../run_OM.sh $SCENARIOS_FOLDER $OM_FOLDER","\n\n", sep ="")
    
    
    # Close the sink!
    sink()
    # either submit run to cluster here or from terminal 
    #submission <<- system(paste0("bash submission",, intern=TRUE)  
  } 
  
  
  
  
  