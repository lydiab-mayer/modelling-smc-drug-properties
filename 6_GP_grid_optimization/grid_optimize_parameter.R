##################################################
##################################################
###                                            ###
### STEP 6: GP-BASED GRID SEARCH OPTIMIZATION  ###
###                                            ###
##################################################
##################################################


### -------------------------------------------------------------------------
###
### M3TPP PROJECT:
### Additional script for running optimization procedure
### of key performance characteristics from a grid search method
### using a pre-train GP emulator
### 
### Original script:
### Created 12.02.2021
### narimane.nekkab@swisstph.ch
### lydia.braunack-mayer@swisstph.ch
### josephine.malinga@swisstph.ch
### 
### R version 3.6.0
###
### -------------------------------------------------------------------------

#INPUTS....

exp = "iTPP3_tradeoffs"
GROUP = "/scicore/home/penny/GROUP/M3TPP/"
SIM_FOLDER = paste0(GROUP, exp, "/")
user = strsplit(getwd(), "/", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]][5]
scale = TRUE
# Create grid of parameter values to search
ngrid <- c(10, 10, 10, 10)
# ngrid = c("Coverage" = ..., "Halflife" = ..., "Efficacy" = ...)
pred = "inc_red_int_Avg"

# Set target range size (1 is per 1% jumps, 10 by 10% etc.)
target_range_size = 10

##############
### HEADER ###
##############

# Set seed for replication
set.seed(42)

# Library
library(dplyr)
library(hetGP)


##################
### PARAMETERS ###
##################

# Load parameter ranges
ranges_file = paste0(GROUP, exp, "/param_ranges.RData")
load(ranges_file)
param_ranges_cont

# Get postprocessing seed file name for GP file name
seed_name = list.files(path = paste0(GROUP,exp,"/postprocessing/"), pattern = "seeds", full.names = FALSE)
seed_name = sub('.txt', '', seed_name) 
seed_name


######################
### OUTPUT FOLDERS ###
######################

# Create output main folder and subfolder per predictor
opt_folder = paste0(SIM_FOLDER, "gp/GP_grid_optimization/")
if(!dir.exists(opt_folder)){
  dir.create(opt_folder)
}

if (!dir.exists(paste0(opt_folder, "/", pred, "/"))) {
  dir.create(paste0(opt_folder, "/", pred))
}


#######################################################
###      GRID SEARCH FUNCTION FOR OPTIMIZATION      ###
#######################################################

predicted = i = "inc_red_int_Avg"
j = seed_name[1]
gp_file = paste0("/scicore/home/penny/GROUP/M3TPP/",exp,"/gp/trained/",i,"/",j,"_", i,"_cv.RData")
optim_param
scale
target_range_size
opt_folder
j

# Create grid search optimization function
GP_grid_search_predictions <- function(predicted, gp_file, scale, ngrid, target_range_size, opt_folder, j, param_ranges_cont){
  
  ###############
  ### SCALING ###
  ###############
  
  if(scale == TRUE){
    # Scale parameter ranges to c(0, 1)
    D <- nrow(param_ranges_cont)
    scale_params <- t(param_ranges_cont)
    for (i in 1:D) {
      scale_params[, i] <- (scale_params[, i] - scale_params[1, i]) / (scale_params[2, i] - scale_params[1, i])
    }
  }else{
    scale_params <- t(param_ranges_cont)
  }
  
  ###################
  ### GRID SEARCH ###
  ###################
  
  # Load GP model
  gp_result_name <- load(gp_file)
  gp_result <- cv_result$GP_model
  
  # Create scenarios
  scenarios <- list()
  
  for (i in 1:D) {
    scenarios[[i]] <- seq(scale_params[1, i], scale_params[2, i], length.out = ngrid[i])
  }
  
  scenarios <- expand.grid(scenarios)
  names(scenarios) <- rownames(param_ranges_cont)
  
  ####################
  ### OPTIMIZATION ###
  ####################
  
  # Make predictions using emulator
  preds <- predict(x = as.matrix(scenarios), object = gp_result)
  scenarios$mean <- preds$mean
  scenarios$sd2 <- preds$sd2
  scenarios$nugs <- preds$nugs
  
  # Covert parameter values back to original scale
  if(scale == TRUE){
    for (i in rownames(param_ranges_cont)) {
      scenarios[, i] <- scenarios[, i] * (param_ranges_cont[i, 2] - param_ranges_cont[i, 1]) + param_ranges_cont[i, 1]
    }
  }
  
  # Identify minimum parameter value required to reach target reduction
  scenarios$target_range <- floor(scenarios$mean/target_range_size)*target_range_size
  
  out <- list()
  
  for (i in rownames(param_ranges_cont)) {
    # Optimization
    optimization_results <- scenarios %>%
      group_by_at(c(rownames(param_ranges_cont)[!(rownames(param_ranges_cont) == i)], "target_range")) %>%
      summarise_at(i, min) %>% 
      ungroup()
    
    # Rename optimized variable
    colnames(optimization_results)[colnames(optimization_results) == i] <- paste0("optimal_", i)
    head(optimization_results)
    
    # Write results
    write.table(optimization_results, file=paste0(opt_folder,predicted,"/",j,"_",predicted,"_GP_grid_optimization_",i,".txt")) 
    
    out[[i]] <- optimization_results  
  }
  
  # Save predictions
  write.table(scenarios, file= paste0(opt_folder,predicted,"/",j,"_",predicted,"_GP_grid_scenario_predictions.txt")) 
  
  return(list(scenarios, optimization_results))
}


#############################################################################
###      Performing optimization with grid search for all predictors      ###
#############################################################################

# All predictors

  # All seeds for each predictor
  
for (i in seed_name){
  
  print(i)

  # Update GP file name
  gp_file = paste0("/scicore/home/penny/GROUP/M3TPP/",exp,"/gp/trained/", pred,"/", i, "_", pred, "_cv.RData")

  # Run
  GP_grid_search_predictions(pred, gp_file, optim_param, scale, n_grid, target_range_size, opt_folder, i)
}


