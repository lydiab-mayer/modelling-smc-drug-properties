########################################
# script genGPtrainscripts.R
#
# creates scripts for running OM and submits the jobs to cluster 
# INPUTS:
#   exp: experiment name
#   chunk_size: batch size for simulation submission
#
# OUTPUTS:
#	- OM scenario xml files and simulations in GROUP folder
#
# Updates
# narimane.nekkab@swisstph.ch
# lydia.braunack-mayer@swisstph.ch
# 15.10.2021
#
########################################

genOMpostprocscripts <- function(exp, date, timesteps){
  
  user <- strsplit(getwd(), "/", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]][5]
  
  if(!dir.exists(paste0("/scicore/home/penny/", user, "/M3TPP/Experiments/", exp, "/JOB_OUT"))){
    dir.create(paste0("/scicore/home/penny/", user, "/M3TPP/Experiments/", exp, "/JOB_OUT"))
  }
  
  GROUP <- "/scicore/home/penny/GROUP/M3TPP/"
  ERROR_FOLDER <- paste0(GROUP, exp, "/postprocessing/err/")
  SIM_FOLDER <- paste0(GROUP, exp, "/")

  file.copy(paste0("/scicore/home/penny/", user, "/M3TPP/analysisworkflow/2_postprocessing/postprocessing_workflow_iTPP3_bloodstage_replication.sh"), 
                   paste0("/scicore/home/penny/", user, "/M3TPP/Experiments/", exp, "/OM_JOBS/postprocessing_workflow_iTPP3_bloodstage_replication.sh"), overwrite = TRUE)
    
  sink(paste0("/scicore/home/penny/", user, "/M3TPP/Experiments/", exp, "/OM_JOBS/job_postprocessing.sh"))
    
  cat("#!/bin/bash","\n", sep ="")
  cat("#SBATCH --job-name=OM_pp","\n", sep ="")
  cat("#SBATCH --account=penny","\n", sep ="")
  cat("#SBATCH -o ",ERROR_FOLDER,"%A_%a.out","\n", sep ="")
  cat("#SBATCH --mem=2G","\n", sep ="")
  cat("#SBATCH --qos=1day","\n", sep ="")
  cat("#SBATCH --cpus-per-task=1","\n", sep ="")
  cat("###########################################","\n", sep ="")
  cat("ml purge","\n", sep ="")
  cat("ml R/3.6.0-foss-2018b","\n", sep ="")
    
  cat("INPUT_DIR=$1", "\n", sep = "")
  cat("OM_RESULTS_DIR=$2", "\n", sep = "")
  cat("DATE=$3", "\n", sep = "")
  cat("TIMESTEPS=$4", "\n", sep = "")
    
  cat("echo \"Input: $INPUT_DIR\"", "\n", sep = "")
  cat("echo \"OM dir: $OM_RESULTS_DIR\"", "\n", sep = "")
  cat("echo \"Date: $DATE\"", "\n", sep = "")
  cat("echo \"Time steps: $TIMESTEPS\"", "\n", sep = "")
    
  # IMPORTANT: the number of files must equal to the task array length (index starts at 0)","\n", sep ="")
  cat("split_files=(${INPUT_DIR}*.txt)","\n", sep ="")
  
  cat("# Select scenario file in array","\n", sep ="")
  cat("ID=$(expr ${SLURM_ARRAY_TASK_ID} - 1)","\n", sep ="")
  cat("split_file=${split_files[$ID]}","\n", sep ="")
  cat("echo \"Postprocessing for $split_file\"","\n", sep ="")
    
  cat("Rscript ../../../analysisworkflow/2_postprocessing/calc_sim_outputs_iTPP3_bloodstage_replication.R $OM_RESULTS_DIR $split_file $DATE $TIMESTEPS","\n", sep ="")
    
  sink()
    
    setwd(paste0("/scicore/home/penny/",user,"/M3TPP/Experiments/",exp,"/OM_JOBS/"))
    
    sys_command = paste("sbatch postprocessing_workflow_iTPP3_bloodstage_replication.sh", SIM_FOLDER, date, timesteps)
    
    # Run  command
    system(sys_command)
  
}
 
