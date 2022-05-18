########################################
# script genOMsimscripts.R
#
# creates scripts for running OM and submits the jobs to cluster 
# INPUTS:
#   exp: experiment name
#   QOS: queue length for simulation submission

# OUTPUTS:
#	- OM scenario xml files and simulations in GROUP folder

########################################

adaptive_design_genOMsimscripts <- function(exp, QOS){
  
  user <- strsplit(getwd(), "/", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]][5]
  user_dir = paste0("/scicore/home/penny/", user, "/M3TPP")
  
  if (!dir.exists(paste0("/scicore/home/penny/", user, "/M3TPP/Experiments/", exp, "/JOB_OUT"))) dir.create(paste0("/scicore/home/penny/", user, "/M3TPP/Experiments/", exp, "/JOB_OUT"))

  GROUP = "/scicore/home/penny/GROUP/M3TPP/"
  
  SIM_FOLDER = paste0(GROUP, exp, "/gp/as/")
  ERROR_FOLDER = paste0(SIM_FOLDER, "err/")
  
  if (!dir.exists(SIM_FOLDER)) dir.create(SIM_FOLDER)
  if (!dir.exists(ERROR_FOLDER)) dir.create(ERROR_FOLDER)
  
  file.copy(paste0(user_dir, "/Experiments/", exp, "/OM_JOBS/scaffold.xml"), paste0(SIM_FOLDER, "/scaffold.xml"), overwrite = TRUE)
  
  #Generate script to create scenarios and run OM simulations
  sink(paste0("/scicore/home/penny/", user, "/M3TPP/Experiments/", exp, "/OM_JOBS/run_OM_as.sh"))
  
  cat("#!/bin/bash","\n", sep ="")
  cat("#SBATCH --job-name=OMSimulation","\n", sep ="")
  cat("#SBATCH --account=penny","\n", sep ="")
  #cat("#SBATCH --error=/dev/null","\n", sep ="") # UNCOMMENT TO HAVE ERROR FILES GENERATED
  #cat("#SBATCH --output=/dev/null","\n", sep ="") # UNCOMMENT TO HAVE ERROR FILES GENERATED
  cat("#SBATCH -o ",ERROR_FOLDER,"%A_%a.out","\n", sep ="") # COMMENT TO HAVE NO ERROR FILES GENERATED
  cat("#SBATCH --mem=2G","\n", sep ="")
  cat("#SBATCH --qos=30min","\n", sep ="")
  cat("#SBATCH --cpus-per-task=1","\n", sep ="")
  cat("###########################################","\n", sep ="")
  cat("# %A_%a.out","\n", sep ="")
  cat("# Script for submitting OpenMalaria simulations as jobs to the cluster","\n", sep ="")
  cat("# Arguments:","\n", sep ="")
  cat("#INPUT_DIR: directory containing the scenario files (.xml)","\n", sep ="")
  cat("#DEST_DIR: directory where two output files will be created per OpenMalaria simulation","\n", sep ="")
  cat("# Calling the script:","\n", sep ="")
  cat("#sbatch --array=1-NUM run_OM_as.sh PARAM_TABLE_FILE SCAFFOLD_FILE BASE_FOLDER INPUT_DIR DEST_DIR","\n", sep ="")
  cat("#","\n", sep ="")
  cat("# created on 01.02.2019","\n", sep ="")
  cat("# monica.golumbeanu@unibas.ch","\n", sep ="")
  cat("###########################################","\n", sep ="")
  
  cat("PARAM_TABLE_FILE=$1","\n", sep ="")
  cat("SCAFFOLD_FILE=$2","\n", sep ="")
  cat("BASE_FOLDER=$3","\n", sep ="")
  cat("INPUT_DIR=$4","\n", sep ="")
  cat("DEST_DIR=$5","\n", sep ="")
  
  cat("# Load R", "\n", sep ="")
  cat("ml R/3.6.0-foss-2018b", "\n", sep = "")
    
  cat("# create destination directory","\n", sep ="")
  cat("mkdir -p $DEST_DIR","\n", sep ="")
  
  cat("ID=$(expr ${SLURM_ARRAY_TASK_ID} - 1)","\n", sep ="")
  cat("# echo \"Debug array ID: \" $ID", "\n", sep ="")
  
  # Generate scenarios
  cat("echo \"Generate replacement patterns and scenario xml\"","\n", sep ="")
  
  cat("# Select the parameter names and the scenario  line from the parameter file", "\n", sep ="")
  cat("column_names=$(sed -n \"1p\" < $PARAM_TABLE_FILE)", "\n", sep ="")
  cat("param_line=$(sed -n \"$(expr ${ID} + 2)p\" < $PARAM_TABLE_FILE)", "\n", sep ="")
  
  cat("echo \"Parameter names: \" $column_names", "\n", sep ="")
  cat("echo \"Parameter values: \" $param_line", "\n", sep ="")
  
  cat("# Construct the replacement patterns", "\n", sep ="")
  cat("Rscript ../../../analysisworkflow/1_OM_basic_workflow/create_scenario.R --column_names $column_names --params $param_line --base_folder $BASE_FOLDER --scenario_folder $INPUT_DIR --scaffold_file $SCAFFOLD_FILE","\n", sep ="")
  cat("echo \"Replacement patterns, base and scenario xml created.\"", "\n", sep ="")
  
  # Load OpenMalaria
  cat("# Load OpenMalaria module and change to folder with resource files","\n", sep ="")
  cat("module purge","\n", sep ="")
  cat("ml OpenMalaria/43.0-iomkl-2019.01","\n", sep ="")
  cat("cd /scicore/home/penny/GROUP/M3TPP/OM_schema43","\n", sep ="")
  
  cat("# ml OpenMalaria/38.0-goolf-1.7.20-Python-2.7.11","\n", sep ="")
  cat("# cd /scicore/home/penny/GROUP/M3TPP/OM_schema38","\n", sep ="")
  cat("# module load OpenMalaria/32-RC3-goolf-1.7.20-Python-2.7.11","\n", sep ="")
  cat("# cd /scicore/home/smith/golmon00/OM_schema32","\n", sep ="")
  
  cat("# IMPORTANT: the number of files must equal to the task array length (index starts at 0)","\n", sep ="")
  
  # Run OpenMalaria simulations
  cat("scenario_file=$(echo $param_line | awk '{print $1}')", "\n", sep ="")
  cat("seed=$(echo $param_line | awk '{print $NF}')", "\n", sep ="")
  cat("# echo \"Debug scenario file: \" $scenario_file", "\n", sep ="")
  cat("# echo \"Debug seed: \" $seed", "\n", sep ="")
  
  cat("scenario_file=${INPUT_DIR}${scenario_file}\"_\"${seed}\".xml\"", "\n", sep = "")
  cat("# echo \"Debug scenario file: \" $scenario_file", "\n", sep ="")
  
  cat("# echo \"Running simulation for $scenario_file\"","\n", sep ="")
    
  cat("# Run OpenMalaria on scenario file","\n", sep ="")
  cat("scenario_name=$(basename \"$scenario_file\" \".xml\")","\n", sep ="")
  cat("output1=$DEST_DIR$scenario_name\"_out.txt\"","\n", sep ="")
  cat("output2=$DEST_DIR$scenario_name\"_cts.txt\"","\n", sep ="")
  cat("echo \"Outputs will be written to $output1 and $output2\"","\n", sep ="")
  cat("openMalaria --scenario $scenario_file --output $output1 --ctsout $output2","\n", sep ="")
  cat("echo \"OpenMalaria simulation ended.\"","\n", sep ="")
  cat("parentdir=\"$(dirname \"$INPUT_DIR\")\"","\n", sep ="")
  # cat("# Remove error files","\n", sep ="")
  # cat(paste0("del \"",ERROR_FOLDER,"%A_%a.out\"","\n"))
  
  sink()

  # Open sink
  sink(paste0("/scicore/home/penny/",user,"/M3TPP/Experiments/",exp,"/OM_JOBS/submission","_",exp,"_as.sh"))
  
  cat("#!/bin/bash\n")
  cat("#SBATCH -o /scicore/home/penny/",user,"/M3TPP/Experiments/",exp,"/JOB_OUT/submission%02i.out","\n",sep ="")
    
  cat("#SBATCH --qos=",QOS,"\n", sep ="") #uncomment to change the queue from 6h to 30min (suitable for small jobs)
    
  cat("PARAM_TABLE_FILE=", SIM_FOLDER, "param_tab.txt", "\n\n", sep ="")
  cat("SCAFFOLD_FILE=", SIM_FOLDER, "scaffold.xml", "\n", sep ="")
  cat("BASE_FOLDER=", SIM_FOLDER, "base/", "\n", sep ="")
  cat("SCENARIOS_FOLDER=", SIM_FOLDER, "scenarios/", "\n", sep ="")
  cat("OM_FOLDER=", SIM_FOLDER, "om/", "\n", sep ="")
    
  cat("NUM=$(wc -l < $PARAM_TABLE_FILE)", "\n", sep ="")

  cat(" echo \"Creating $NUM-1 scenarios and running OM simulations... \" ","\n", sep ="")
    
  cat("sbatch -W --array=1-$NUM%380 ./run_OM_as.sh $PARAM_TABLE_FILE $SCAFFOLD_FILE $BASE_FOLDER $SCENARIOS_FOLDER $OM_FOLDER","\n\n", sep ="")
    
  # Close the sink!
  sink()
  
  # either submit run to cluster here or from terminal 
  setwd(paste0("/scicore/home/penny/",user,"/M3TPP/Experiments/",exp,"/OM_JOBS/"))
  sys_command = paste0("sbatch", " submission_",exp,"_as.sh")
  
  # Run  command
  system(sys_command)
  
} 
