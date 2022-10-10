############################################################
# import_functions
#
# Helper functions for plotting scripts within M3TPP project workflow
#
# Written by Lydia Burgert
# Adapted by Lydia Braunack-Mayer
############################################################

source(paste0("/scicore/home/penny/",user,"/M3TPP/analysisworkflow/2_postprocessing/postprocessing_resources.R"))

# ----------------------------------------------------------
# import_cont
# ----------------------------------------------------------

import_cont <- function(exp,scenario_name,param_table_file,seeds,follow_up,years_before_interv){

  om_results_folder <- paste0("/scicore/home/penny/GROUP/M3TPP/",exp,"/om/")
  param_table <- read.table(param_table_file, sep= "\t", header = TRUE, as.is = TRUE, stringsAsFactors = FALSE)
  

  
  df_list_prev_210<- list()
  df_list_prev_int<- list()
  df_list_prev_allages<- list()
  df_list_prev_agegroups<- list()
  df_list_incidence_05<- list()
  df_list_incidence_int<- list()
  df_list_incidence_allages<- list()
  df_list_incidence_agegroups<- list()
  df_list_incidence_int_5mo<- list()

scenario<- scenario_name  
  for (k in 1:seeds) { 
  OM_result_file = paste(om_results_folder, scenario, "_",
                         +                                k, "_out.txt", sep="")
  OM_result = read.table(OM_result_file, sep="\t",header=TRUE)
  
  scenario_params <- param_table[which(param_table$Scenario_Name==scenario),][1,]
 
  test <-  calculate_outputs(OM_result, scenario_params, follow_up,years_before_interv,cont=TRUE)
  
  df_list_prev_210[[k]] <-  test$prevalence_210 
  df_list_prev_int[[k]] <-  test$prevalence_int 
  df_list_prev_allages[[k]] <-  test$prevalence_allages 
  df_list_prev_agegroups[[k]] <-  test$prevalence_agegroups 
  df_list_incidence_05[[k]] <-  test$incidence_05 
  df_list_incidence_int[[k]] <-  test$incidence_int 
  df_list_incidence_allages[[k]] <-  test$incidence_allages 
  df_list_incidence_agegroups[[k]] <-  test$incidence_agegroups 
  df_list_incidence_int_5mo[[k]] <-  test$incidence_int_5mo 
  
  }
  
  
  df_prevalence_210 <- rbindlist(df_list_prev_210)[,lapply(.SD,mean), list(time)]  
  df_prevalence_int <- rbindlist(df_list_prev_int)[,lapply(.SD,mean), list(time)]  
  df_prevalence_allages<- rbindlist(df_list_prev_allages)[,lapply(.SD,mean), list(time)]  
  df_prevalence_agegroups <- rbindlist(df_list_prev_agegroups)[,lapply(.SD,mean), list(time)]  
  
  df_incidence_05 <- rbindlist(df_list_incidence_05)[,lapply(.SD,mean), list(time)]  
  df_incidence_int<- rbindlist(df_list_incidence_int)[,lapply(.SD,mean), list(time)]  
  df_incidence_allages<- rbindlist(df_list_incidence_allages)[,lapply(.SD,mean), list(time)]  
  df_incidence_agegroups<- rbindlist(df_list_incidence_agegroups)[,lapply(.SD,mean), list(time)]  
  df_incidence_int_5mo<- rbindlist(df_list_incidence_int_5mo)[,lapply(.SD,mean), list(time)]  
  
  
   
  return(list(
    "prevalence_210"=df_prevalence_210,
    "prevalence_int"=df_prevalence_int,
    "prevalence_allages"=df_prevalence_allages,
    "prevalence_agegroups"=df_prevalence_agegroups,
    "incidence_05"=df_incidence_05,
    "incidence_int"=df_incidence_int,
    "incidence_allages"=df_incidence_allages,
    "incidence_agegroups"=df_incidence_agegroups,
    "incidence_int_5mo"=df_incidence_int_5mo
  ))
} 

# ----------------------------------------------------------
# import_cont_4var
# ----------------------------------------------------------

import_cont_4var <- function(exp,scenario_name,param_table_file,seeds,follow_up,years_before_interv){
  
  om_results_folder <- paste0("/scicore/home/penny/GROUP/M3TPP/",exp,"/om/")
  param_table <- read.table(param_table_file, sep= "\t", header = TRUE, as.is = TRUE, stringsAsFactors = FALSE)
  
  
  
  df_list_prev_210<- list()
  df_list_prev_int<- list()
  df_list_incidence_int<- list()
  df_list_incidence_int_5mo<- list()
  
  scenario<- scenario_name  
  for (k in 1:seeds) { 
    OM_result_file = paste(om_results_folder, scenario, "_",
                           +                                k, "_out.txt", sep="")
    OM_result = read.table(OM_result_file, sep="\t",header=TRUE)
    
    scenario_params <- param_table[which(param_table$Scenario_Name==scenario),][1,]
    
    test <-  calculate_outputs(OM_result, scenario_params, follow_up,years_before_interv,cont=TRUE)
    
    df_list_prev_210[[k]] <-  test$prevalence_210 
    df_list_prev_int[[k]] <-  test$prevalence_int 
    df_list_incidence_int[[k]] <-  test$incidence_int 
    df_list_incidence_int_5mo[[k]] <-  test$incidence_int_5mo 
    
  }
  
  
  df_prevalence_210 <- rbindlist(df_list_prev_210)[,lapply(.SD,mean), list(time)]  
  df_prevalence_int <- rbindlist(df_list_prev_int)[,lapply(.SD,mean), list(time)]  

  df_incidence_int<- rbindlist(df_list_incidence_int)[,lapply(.SD,mean), list(time)]  
  df_incidence_int_5mo<- rbindlist(df_list_incidence_int_5mo)[,lapply(.SD,mean), list(time)]  
  
  
  
  return(list(
    "prevalence_210"=df_prevalence_210,
    "prevalence_int"=df_prevalence_int,
    "incidence_int"=df_incidence_int,
    "incidence_int_5mo"=df_incidence_int_5mo
  ))
} 

# ----------------------------------------------------------
# import_EIRs_cat
# ----------------------------------------------------------

import_EIRs_cat <- function(exp, scenario_id, seeds, timesteps){
  
  # ----------------------------------------------------------
  # This function averages continuous OpenMalaria outputs from multiple random seeds, by setting
  #
  # Inputs
  # exp: string containing the experiment name as it appears on SciCore, e.g. "MyExperiment"
  # scenario_id: character vector containing id(s) for scenario(s) to plot, e.g c("MyExperiment_1", "MyExperiment_2")
  # seeds: integer vector containing selection of seeds to plot, e.g. c(1)
  # timesteps: integer vector containing selection of OpenMalaria timesteps to plot, e.g. 1:73 captures outputs for the first year of simulation
  #
  # Outputs
  # out: list containing 'All', data frame containing all continuous OpenMalaria outputs for the given setting, and
  # 'Average', data frame containing averaged continuous OpenMalaria outputs for the given setting
  #
  # ----------------------------------------------------------
  
  # Load packages
  require(data.table)
  
  # Set up file pathways
  om_results_folder <- paste0("/scicore/home/penny/GROUP/M3TPP/", exp, "/om/")
  
  # Set up list to store function outputs
  out_all <- list()
  
  # Identify and concatenate continuous OpenMalaria outputs for the given setting, seed and selection of timesteps
  for (j in scenario_id) {
    for (i in seeds) {
      om_result <- read.table(paste0(om_results_folder, j, "_", i, "_cts.txt"), sep = "\t", header = TRUE)
      om_result$scenario_id <- j; om_result$seed <- i
      out_all[[paste0("scenario_id", j, "_seed", i)]]<- subset(om_result, timestep %in% timesteps)
    }
  }
  
  # Average continuous OpenMalaria outputs across seeds
  out_all <- rbindlist(out_all)
  out_average <- out_all[, lapply(.SD, mean), list(scenario_id, timestep)]
  
  # Return function outputs
  return(list("All" = out_all, "Average" = out_average))
  
}

# ----------------------------------------------------------
# import_monitoring_outcome
# ----------------------------------------------------------

import_monitoring_outcome <- function(exp, scenario_id, measure, seeds, timesteps){
  
  # ----------------------------------------------------------
  # This function averages OpenMalaria monitoring outcomes from multiple random seeds, by setting and age group
  #
  # Inputs
  # exp: string containing the experiment name as it appears on SciCore, e.g. "MyExperiment"
  # scenario_id: character vector containing id(s) for scenario(s) to plot, e.g c("MyExperiment_1", "MyExperiment_2")
  # measure: integer containing the id for the survey measure to plot (see https://github.com/SwissTPH/openmalaria/wiki/MonitoringOptions)
  # seeds: integer vector containing selection of seeds to plot, e.g. c(1)
  # timesteps: integer vector containing selection of OpenMalaria timesteps to plot, e.g. 1:73 captures outputs for the first year of simulation
  #
  # Outputs
  # out: list containing 'All', data frame containing the selected OpenMalaria monitoring outcome for all simulations for a given setting, and
  # 'Average', data frame containing averaged OpenMalaria monitoring outcome for the given setting
  #
  # ----------------------------------------------------------
  
  # Load packages
  require(data.table)
  
  # Set up file pathways
  om_results_folder <- paste0("/scicore/home/penny/GROUP/M3TPP/", exp, "/om/")
  
  # Set up list to store function outputs
  out_all <- list()
  
  # Identify and concatenate continuous OpenMalaria outputs for the given setting, seed and selection of timesteps
  for (j in scenario_id) {
    for (i in seeds) {
      om_result <- read.table(paste0(om_results_folder, j, "_", i, "_out.txt"), sep = "\t", header = FALSE)
      names(om_result) <- c("timestep", "agegroup", "measures", "outcome")
      om_result$scenario_id <- j; om_result$seed <- i
      out_all[[paste0("scenario_id", j, "_seed", i)]] <- subset(om_result, timestep %in% timesteps & measures %in% measure)
    }
  }
  
  # Average continuous OpenMalaria outputs across seeds, by age group
  out_all <- rbindlist(out_all)
  out_average <- out_all[, lapply(.SD, mean), list(scenario_id, timestep, agegroup)]
  
  # Return function outputs
  return(list("All" = out_all, "Average" = out_average))
  
}

