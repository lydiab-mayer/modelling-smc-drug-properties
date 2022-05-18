####################################################
# Adaptive sampling scheme for a multi-dimensional space
#
# 
# monica.golumbeanu@unibas.ch
# Adapted by lydia.braunack-mayer@swisstph.ch
###################################################

##############
### SET UP ###
##############

# Set working directory
user = strsplit(getwd(), "/", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]][5]
setwd(paste0("/scicore/home/penny/", user, "/M3TPP"))

# Load required packages
library(hetGP)
library(tgp)
library(dplyr)

# Define arguments
args = commandArgs(TRUE)
as_folder = args[1]
gp_file = args[2]
scale = args[3]

# Retained here for testing
as_folder = "/scicore/home/penny/GROUP/M3TPP/iTPP3_bloodstage_mixed/gp/as/"
gp_file = "/scicore/home/penny/GROUP/M3TPP/iTPP3_bloodstage_mixed/gp/trained/inc_red_int_Tot/seeds_iTPP3bloodstagemixed_seas3mo_sp_Mali_2_5_0.04_May_0.020831339_inc_red_int_Tot_cv.RData"
scale = TRUE


########################
### HELPER FUNCTIONS ###
########################

source("~/M3TPP/analysisworkflow/3_GP_train/adaptivesampling/adaptive_design_genOMsimscripts.R")
source("~/M3TPP/analysisworkflow/2_postprocessing/postprocessing_resources.R")
source("~/M3TPP/analysisworkflow/3_GP_train/GP_toolbox.R")

# Function that returns a set of points to be explored with OpenMalaria and added to the training set
# Method: uniform sampling across the parameter space
select_points_uniform = function(GP_model, param_ranges, num_points, n_sample = 100000) {
  
  # Select a random set of points uniformly distributed across the parameter space
  X_random_samples = lhs(n_sample, param_ranges) 
  colnames(X_random_samples) = rownames(param_ranges)
  
  # Predict the output with the GP model on the new points and evaluate the variance   
  pred_obj = predict(x = X_random_samples, object = GP_model)
  new_data = cbind.data.frame(X_random_samples, sqrt(pred_obj$sd2 + pred_obj$nugs), row.names = NULL)
  colnames(new_data) = c(colnames(X_random_samples), "variance")
  
  # Select a set of points with the largest posterior predictive variance
  XX = new_data[order(new_data$variance, decreasing = TRUE), ]
  
  # Return function outputs
  param_tab = XX[1:num_points,]
  variances_vec = c(min(sqrt(pred_obj$sd2 + pred_obj$nugs)), mean(sqrt(pred_obj$sd2 + pred_obj$nugs)), max(sqrt(pred_obj$sd2 + pred_obj$nugs)))
  
  return(list("param_tab" = param_tab, "variances_vec" = variances_vec))
}

# Create parameter table from newly sampled points to be used for building corresponding scenarios 
create_param_table_from_samples = function(sampled_points, split_tab) {
  
  # Remove variance column
  sampled_points$variance = NULL
  
  # load categorical parameter values
  col_id = min(which(names(split_tab) %in% names(sampled_points) == TRUE))
  as_param_tab = split_tab[, c(2:(col_id - 1))]
  as_param_tab = unique(as_param_tab)
  
  # add setting parameters
  as_param_tab = cbind(as_param_tab, sampled_points, row.names = NULL)

  # add seeds and define scenarios names
  SEED = unique(split_tab$seed)
  Scenario_Name = paste("Scenario", 1:nrow(as_param_tab), sep="_")
  
  # construct parameter table
  as_param_tab = cbind(Scenario_Name, as_param_tab, row.names = NULL)
  as_param_tab = merge(as_param_tab, as.data.frame(SEED))
  
  return(as_param_tab)
}

# Function that waits until all jobs have finished (marked in a log file)
wait_for_jobs = function(log_file, n_lines) {
  while((!file.exists(log_file))) {
    # waiting for log file to be created
    Sys.sleep(1)
  }
  while(nrow(read.table(log_file)) < n_lines ) {
    # waiting for scenarios to be created
    Sys.sleep(1)
  }
  return(0)
}


#############################
### RUN ADAPTIVE SAMPLING ###
#############################

## Set up
if (!file.exists(as_folder)) {dir.create(as_folder)}
if(!file.exists(paste0(as_folder, "/postprocessing"))) {dir.create(paste0(as_folder, "/postprocessing"))}
exp_dir = dirname(dirname(as_folder))

## Load the GP model
load(gp_file)
predicted = cv_result$predicted

## Load the parameter table
split_file = gsub("gp/trained/inc_red_int_Tot", "postprocessing", gp_file)
split_file = gsub(paste0("_", predicted, "_cv.RData"), ".txt", split_file)
split_tab = read.table(split_file, header = TRUE)

#if(!file.exists(paste0(as_folder, "/", cv_result$predicted))) {dir.create(paste0(as_folder, "/", cv_result$predicted))}

## Load parameter ranges
load(paste0(exp_dir, "/param_ranges.RData"))

## Get experimant name
exp = basename(exp_dir)

## Initialization

# Set number of new sampled points for updating the GP model
n_samples = 10

# Initalize variances vector
variances = NULL

## Adaptive sampling
for(as_runs in 1:2) {
  
  cat("Initiate adaptive sampling run", as_runs)

  # Obtain new samples and create the parameter table with new data points
  new_points_obj = select_points_uniform(cv_result$GP_model, param_ranges_cont, n_samples)
  param_tab = create_param_table_from_samples(new_points_obj$param_tab, split_tab) 
  
  # Store files
  tab_file = paste0(as_folder, "param_tab.txt")
  
  write.table(param_tab, tab_file, sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)
  variances = rbind(variances, new_points_obj$variances_vec)
  
  # Create the new scenarios and submit OpenMalaria simulations to cluster
  cat("Generate OpenMalaria simulations run", as_runs)
  adaptive_design_genOMsimscripts(exp, QOS = "30min")
  
  # Wait until simulations are complete
  count <- nrow(param_tab)
  while (count < length(list.files(paste0(as_folder, "om/"), pattern = "*_out.txt"))) Sys.sleep(10)
  
  # Run postprocessing, get the new training data and update GP
  cat("Run OM postprocessing run", as_runs)
  pp_params = readRDS(paste0(exp_dir, "/postprocessing/pp_param_values.rds"))
  postprocess.om(dir = paste0(as_folder, "om/"), 
                 param.file = tab_file,
                 date = pp_params$date,
                 fmonth = pp_params$fmonth,
                 months = pp_params$months,
                 year.counterfactual = pp_params$year_counterfactual,
                 year.intervention = pp_params$year_intervention,
                 min.int = pp_params$min_int)
  pp_results = read.table(paste0(as_folder, "/postprocessing/seeds_param_tab.txt"), header = TRUE)
  
  # Prepare data for training emulator
  new_train_data = pp_results[, c(rownames(param_ranges_cont), predicted)]
  n = ncol(new_train_data)
  
  if (scale) {
    for (col in 1:(n - 1)) {
      new_train_data[, col] = (new_train_data[, col] - param_ranges_cont[col, 1]) / (param_ranges_cont[col, 2] - param_ranges_cont[col, 1])
    }
  }
  
  prdata = find_reps(X = as.matrix(new_train_data[, 1:(n-1)]), Z = as.matrix(new_train_data[, n]),
                     rescale = FALSE, normalize = FALSE)
  
  # Update emulator with new data and model fit
  cat("Train emulator run", as_runs)
  train_data = data.frame(rbind(cv_result$train_data, new_train_data))
  cv_result = cv_train_matern(input_data = train_data, 
                              lower = rep(0.001, nrow(param_ranges_cont)),
                              upper = rep(10, nrow(param_ranges_cont)),
                              scale = NULL,
                              test_prop = 0.1,
                              hetGP = TRUE)

  # Clean up simulation files before the next OpenMalaria run 
  system(paste0("rm -r ", as_folder, "err/"))
  system(paste0("rm -r ", as_folder, "om/"))
  system(paste0("rm -r ", as_folder, "base/"))
  system(paste0("rm -r ", as_folder, "scenarios/"))
}

# Write results to file
cat("Adapting sampling completed. Writing results to file")
gp_save = paste0(gsub("_cv.RData", "", gp_file), "_as")
test_GP_plot(cv_result = cv_result, save = paste0(gp_save, ".jpg"))
save(cv_result, file = paste0(gp_save, ".RData"))

