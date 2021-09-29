########################################
# script genoptimizationscripts.R
#
# creates scripts for running sensitvity analysis and submits the jobs to cluster 
# INPUTS:
#   exp: experiment name
#   chunk_size: batch size for simulation submission

# OUTPUTS:
#	- OM scenario xml files and simulations in GROUP folder

########################################

genoptimizationscripts <- function(exp, predicted, reductions, optimized, n_gridpoints, scale){
  

  user <- strsplit(getwd(), "/", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]][5]
  
  
  GROUP = "/scicore/home/penny/GROUP/M3TPP/"
  
  SIM_FOLDER=paste0(GROUP,exp,"/")
  ERROR_FOLDER=paste0(GROUP,exp,"/gp/optimisation/err/")

  opt_folder <- paste0(SIM_FOLDER,"gp/optimisation/")
  dir.create(opt_folder)
  dir.create(paste0(opt_folder,predicted,"/"))
  dir.create(ERROR_FOLDER)
  
  param_ranges_file = paste0(SIM_FOLDER, "param_ranges.RData")
  table_file = paste0(SIM_FOLDER, "opt_setup.txt")
  load(param_ranges_file)
  
  Param_opt = c(optimized)
  Param_vec = t(param_ranges_cont[,1])
  opt_setup = Reduce(merge, list(as.data.frame(reductions), as.data.frame( Param_vec ), as.data.frame(Param_opt)))
  
  table_file <- paste0(opt_folder,predicted, "/",Param_opt,"_",reductions,"_opt_setup.txt")
 
   write.table(opt_setup, table_file, sep = "\t", quote = FALSE, col.names = TRUE,
              row.names = FALSE)
  
  
  
 file.copy(paste0("/scicore/home/penny/",user,"/M3TPP/analysisworkflow/5_optimization/optimise_workflow.sh"), 
            paste0("/scicore/home/penny/",user,"/M3TPP/Experiments/",exp,"/OM_JOBS/optimise_workflow.sh"),overwrite=TRUE)
  
  
  
sink(paste0("/scicore/home/penny/",user,"/M3TPP/Experiments/",exp,"/OM_JOBS/job_optimise_parameter.sh"))
  
cat("#!/bin/bash","\n", sep ="")
cat("#SBATCH --job-name=optimisation","\n", sep ="")
cat("#SBATCH --account=penny","\n", sep ="")
cat("#SBATCH -o ",ERROR_FOLDER,"%A_%a.out","\n", sep ="")
# cat("#SBATCH -o /scicore/home/penny/",user,"/M3TPP/Experiments/",exp,"/JOB_OUT/5_optimization.out","\n", sep ="")
cat("#SBATCH --mem=200G","\n", sep ="")
cat("#SBATCH --qos=1day","\n", sep ="")
cat("#SBATCH --cpus-per-task=4","\n", sep ="")
cat("###########################################","\n", sep ="")
cat("ml purge","\n", sep ="")
cat("ml R/3.6.0-foss-2018b","\n", sep ="")

cat("GP_DIR=$1","\n", sep ="")
cat("PARAM_RANGES_FILE=$2","\n", sep ="")
cat("SENS_DEST_DIR=$3","\n", sep ="")
cat("opt_setup_file=$4","\n", sep ="")
cat("n_gridpoints=$5","\n", sep ="")
cat("SCALE=$6","\n", sep ="")

cat("gp_files=(${GP_DIR}*.RData)","\n", sep ="")

cat("# Select scenario file in array","\n", sep ="")
cat("ID=$(expr ${SLURM_ARRAY_TASK_ID} - 1)","\n", sep ="")
cat("gp_file=${gp_files[$ID]}","\n", sep ="")
cat("echo $gp_file","\n", sep ="")


cat("ml R/3.6.0-foss-2018b","\n", sep ="")
cat("Rscript ../../../analysisworkflow/5_optimization/optimise_parameter.R $gp_file $PARAM_RANGES_FILE $SENS_DEST_DIR $opt_setup_file $n_gridpoints $SCALE","\n", sep ="")

sink()


setwd(paste0("/scicore/home/penny/",user,"/M3TPP/Experiments/",exp,"/OM_JOBS/"))

  

sys_command = paste("bash optimise_workflow.sh", SIM_FOLDER ,predicted, table_file, n_gridpoints, scale)

# Run  command
system(sys_command)

}
