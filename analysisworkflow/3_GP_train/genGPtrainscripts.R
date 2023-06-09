########################################
# script genGPtrainscripts.R
#
# Contains function that generates scripts for training the emulator and submits the jobs to cluster 
# INPUTS:
#   exp: experiment name
#   chunk_size: batch size for simulation submission

# OUTPUTS:
#	  OM scenario xml files and simulations in GROUP folder
#
# UPDATES:
# lydia.braunack-mayer@swisstph.ch
# 08.10.2021
#
########################################

genGPtrainscripts <- function(exp, predicted, lower, upper, scale, hetGP){
  
  # Set up
  user <- strsplit(getwd(), "/", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]][5]
  GROUP <- "/scicore/home/penny/GROUP/M3TPP/"
  ERROR_FOLDER = paste0(GROUP,exp,"/gp/trained/err/")
  SIM_FOLDER <- paste0(GROUP,exp,"/")
  
  # Generate .sh file to train each emulator
  sink(paste0("/scicore/home/penny/",user,"/M3TPP/Experiments/",exp,"/OM_JOBS/job_train_GP.sh"))
  
  cat("#!/bin/bash","\n", sep ="")
  cat("#SBATCH --job-name=GPtrain","\n", sep ="")
  cat("#SBATCH --account=penny","\n", sep ="")
  cat("#SBATCH -o ", ERROR_FOLDER, "%A_%a.out","\n", sep ="")
  cat("#SBATCH --mem=2G","\n", sep ="")
  cat("#SBATCH --qos=30min","\n", sep ="")
  cat("#SBATCH --cpus-per-task=1","\n", sep ="")
  
  cat("###########################################","\n", sep ="")
  
  cat("ml purge","\n", sep ="")
  cat("ml R/3.6.0-foss-2018b","\n", sep ="")
  
  cat("INPUT_DIR=$1","\n", sep ="")
  cat("DEST_DIR=$2","\n", sep ="")
  cat("PREDICTED=$3","\n", sep ="")
  cat("RANGES_FILE=$4","\n", sep ="")
  cat("LOWER=$5","\n", sep ="")
  cat("UPPER=$6","\n", sep ="")
  cat("SCALE=$7","\n", sep ="")
  cat("HETGP=$8", "\n", sep = "")
  
  cat("# IMPORTANT: the number of files must equal to the task array length (index starts at 0)","\n", sep ="")
  cat("if [ $HETGP = \"TRUE\" ]", "\n", sep = "")
  cat("then", "\n", sep = "")
  cat("setting_postprocessing_results=(${INPUT_DIR}seeds_*.txt)","\n", sep ="")
  cat("else", "\n", sep = "")
  cat("setting_postprocessing_results=(${INPUT_DIR}agg_*.txt)","\n", sep ="")
  cat("fi", "\n", sep = "")
  
  cat("# Select scenario file in array","\n", sep ="")
  cat("ID=$(expr ${SLURM_ARRAY_TASK_ID} - 1)","\n", sep ="")
  cat("setting_postprocessing_result=${setting_postprocessing_results[$ID]}","\n", sep ="")
  cat("echo \"Postprocessing $PREDICTED for $setting_postprocessing_result\"","\n", sep ="")
  
  cat("Rscript ../../../analysisworkflow/3_GP_train/train_GP.R $setting_postprocessing_result $DEST_DIR $PREDICTED $RANGES_FILE $LOWER $UPPER $SCALE $HETGP","\n", sep ="")
  
  sink()
  
  # Prepare .sh file to submit job_train_GP.sh to cluster
  file.copy(paste0("/scicore/home/penny/",user,"/M3TPP/analysisworkflow/3_GP_train/GP_train_workflow.sh"), 
            paste0("/scicore/home/penny/",user,"/M3TPP/Experiments/",exp,"/OM_JOBS/GP_train_workflow.sh"),overwrite=TRUE)
  
  setwd(paste0("/scicore/home/penny/",user,"/M3TPP/Experiments/",exp,"/OM_JOBS/"))
  
  sys_command <- paste("sbatch GP_train_workflow.sh", SIM_FOLDER, predicted, lower, upper, scale, hetGP)
  
  # Submit job to cluster
  system(sys_command)
  
}
