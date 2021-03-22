########################################
# script create_folders.R
#
# creates folders for storing simulation outputs
#
########################################

create_folders <- function(exp){
  user <- strsplit(getwd(), "/", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]][5]
  GROUP = "/scicore/home/penny/GROUP/M3TPP/"
  
  dir.create(paste0(GROUP,exp))
  dir.create(paste0("./Experiments/",exp))
  dir.create(paste0("./Experiments/",exp,"/JOB_OUT"))
  dir.create(paste0("./Experiments/",exp,"/OM_JOBS"))
  
  file.copy(paste0("/scicore/home/penny/",user,"/M3TPP/analysisworkflow/1_OM_basic_workflow/1_OM_base_workflow_original.R"), 
            paste0("/scicore/home/penny/",user,"/M3TPP/Experiments/",exp,"/1_OM_base_workflow_",exp,".R"),overwrite=TRUE)
  
  file.copy(paste0("/scicore/home/penny/",user,"/M3TPP/analysisworkflow/2_postprocessing/2_OM_postprocessing_original.R"), 
            paste0("/scicore/home/penny/",user,"/M3TPP/Experiments/",exp,"/2_OM_postprocessing_",exp,".R"),overwrite=FALSE)
  
  file.copy(paste0("/scicore/home/penny/",user,"/M3TPP/analysisworkflow/3_GP_train/3_GP_train_workflow_original.R"), 
            paste0("/scicore/home/penny/",user,"/M3TPP/Experiments/",exp,"/3_GP_train_workflow_",exp,".R"),overwrite=FALSE)
  
  file.copy(paste0("/scicore/home/penny/",user,"/M3TPP/analysisworkflow/4_sensitivity_analysis/4_sensitivityanalsis_workflow_original.R"), 
            paste0("/scicore/home/penny/",user,"/M3TPP/Experiments/",exp,"/4_sensitivityanalysis_workflow_",exp,".R"),overwrite=FALSE)
  
  file.copy(paste0("/scicore/home/penny/",user,"/M3TPP/analysisworkflow/5_optimization/5_optimisation_workflow_original.R"), 
            paste0("/scicore/home/penny/",user,"/M3TPP/Experiments/",exp,"/5_optimisation_workflow_",exp,".R"),overwrite=FALSE)
  
  
  # copying plotting files 
  
  dir.create(paste0("./Experiments/",exp,"/analysis_scripts/"))
  dir.create(paste0("./Experiments/",exp,"/Outputs/"))
  
  file.copy(paste0("/scicore/home/penny/",user,"/M3TPP/analysisworkflow/analysis_scripts/PS_00_plotseasonality_Original.R"), 
            paste0("/scicore/home/penny/",user,"/M3TPP/Experiments/",exp,"/PS_00_plotseasonality_",exp,".R"),overwrite=TRUE)
  
  file.copy(paste0("/scicore/home/penny/",user,"/M3TPP/analysisworkflow/analysis_scripts/PS_02_CompareSimulationResults_Original.R"), 
            paste0("/scicore/home/penny/",user,"/M3TPP/Experiments/",exp,"/PS_02_CompareSimulationResults_",exp,".R"),overwrite=TRUE)
  
}