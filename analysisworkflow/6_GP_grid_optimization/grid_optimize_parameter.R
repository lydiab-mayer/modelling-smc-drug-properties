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
### using a pre-trained GP emulator
### 
### Original script:
### Created 29.10.2021
### narimane.nekkab@swisstph.ch
### lydia.braunack-mayer@swisstph.ch
### josephine.malinga@swisstph.ch
### 
### R version 3.6.0
###
### -------------------------------------------------------------------------

##############
### HEADER ###
##############

args <- commandArgs(TRUE)
gp_file <- args[1]
sim_folder <- args[2]
scale <- as.logical(args[3])
ngrid <- as.numeric(args[4])
target_range_size <- as.numeric(args[5])

# # Sample arguments, retained here for testing
# sim_folder <- "/scicore/home/penny/GROUP/M3TPP/iTPP3_ChemoLiver_TreatLiverBlood_3rounds/"
# scale <- TRUE
# ngrid <- "100"
# target_range_size <- 10
# gp_file <- "/scicore/home/penny/GROUP/M3TPP/iTPP3_ChemoLiver_TreatLiverBlood_3rounds/gp/trained/inc_red_int_Tot/seeds_iTPP3ChemoLiverTreatLiverBlood3rounds_seas3mo_Mali_2_5_exp_0.04_May_inc_red_int_Tot_cv.RData"
# ngrid <- as.numeric(ngrid)

print(paste0("gp_file: ", gp_file))
print(paste0("sim_folder: ", sim_folder))
print(paste0("scale: ", scale))
print(paste0("ngrid: ", ngrid))
print(paste0("target_range_size: ", target_range_size))

# Library
library(dplyr)
library(hetGP)
library(tgp)


##################
### PARAMETERS ###
##################

# Load parameter ranges
print("Loading parameter ranges")
load(paste0(sim_folder, "param_ranges.RData"))
param_ranges_cont
# param_ranges <- param_ranges_cont # uncomment for testing


#######################################################
###      GRID SEARCH FUNCTION FOR OPTIMIZATION      ###
#######################################################

# Create grid search optimization function
GP_grid_search_predictions <- function(gp_file, scale, ngrid, target_range_size, param_ranges){
  
  ###############
  ### SCALING ###
  ###############
  
  if(scale == TRUE){
    # Scale parameter ranges to c(0, 1)
    D <- nrow(param_ranges)
    scale_params <- t(param_ranges)
    for (i in 1:D) {
      scale_params[, i] <- (scale_params[, i] - scale_params[1, i]) / (scale_params[2, i] - scale_params[1, i])
    }
  }else{
    scale_params <- t(param_ranges)
  }
  
  ###################
  ### GRID SEARCH ###
  ###################
  
  # Load GP model
  gp_result_name <- load(gp_file)
  gp_result <- cv_result$GP_model
  
  # Generate sample scenarios
  
  # Uncomment to sample uniformly across the entire parameter space
  scenarios <- lhs(ngrid, t(scale_params))
  scenarios <- as.data.frame(scenarios)
  names(scenarios) <- rownames(param_ranges)
  # end uncomment
  
  # # Uncomment to sample uniformly across drug parameters only
  # if (scale == TRUE) {
  #   params <- expand.grid("Coverage1" = seq(0.05, 0.25, 0.05)/0.3,
  #                           "Coverage2" = seq(0.05, 0.25, 0.05)/0.3,
  #                           "Slope" = (6 - 0.5)/(9 - 0.5))
  # } else {
  #   params <- expand.grid("Coverage1" = seq(0.75, 0.95, 0.05), "Coverage2" = seq(0.75, 0.95, 0.05), "Slope" = 6)
  # }
  # ncov <- nrow(params)
  # params <- params[rep(seq_len(ncov), ngrid), ]
  # 
  # samp <- lhs(ngrid, t(scale_params[, c("Halflife", "MaxKillingRate")]))
  # samp <- samp[rep(seq_len(ngrid), ncov), ]
  # colnames(samp) <- c("Halflife", "MaxKillingRate")
  # 
  # scenarios <- cbind(params, samp)
  # scenarios <- scenarios[, rownames(param_ranges)]
  # rownames(scenarios) <- 1:nrow(scenarios)
  # # end uncomment
  
  # # Uncomment to sample uniformly across all parameters but one (slope)
  # if (scale == TRUE) {
  #   params <- expand.grid("Slope" = (6 - 0.5)/(9 - 0.5))
  # } else {
  #   params <- expand.grid("Slope" = 6)
  # }
  # ncov <- nrow(params)
  # params <- params[rep(seq_len(ncov), ngrid), ]
  # 
  # samp <- lhs(ngrid, t(scale_params[, c("Coverage1", "Coverage2", "Halflife", "MaxKillingRate")]))
  # samp <- samp[rep(seq_len(ngrid), ncov), ]
  # colnames(samp) <- c("Coverage1", "Coverage2", "Halflife", "MaxKillingRate")
  # 
  # scenarios <- cbind("Slope" = params, samp)
  # scenarios <- as.data.frame(scenarios[, rownames(param_ranges)])
  # rownames(scenarios) <- 1:nrow(scenarios)
  # # end uncomment
  
  ####################
  ### OPTIMIZATION ###
  ####################
  
  scenarios$nugs <- scenarios$sd2 <- scenarios$mean <- NA

  # Make predictions using emulator - split into loop to reduce memory requirements
  preds <- predict(x = as.matrix(scenarios[, rownames(param_ranges)]), object = gp_result)
  scenarios$mean <- preds$mean
  scenarios$sd2 <- preds$sd2
  scenarios$nugs <- preds$nugs
  
  # Covert parameter values back to original scale
  if(scale == TRUE){
    for (i in rownames(param_ranges)) {
      scenarios[, i] <- scenarios[, i] * (param_ranges[i, 2] - param_ranges[i, 1]) + param_ranges[i, 1]
    }
  }
  
  # Identify minimum parameter value required to reach target reduction
  scenarios$target_range <- floor(scenarios$mean/target_range_size)*target_range_size
  
  out <- list()
  
  for (i in rownames(param_ranges_cont)) {
    # Optimization
    optimization_results <- scenarios %>%
      group_by_at(c(rownames(param_ranges)[!(rownames(param_ranges) == i)], "target_range")) %>%
      summarise_at(i, min) %>% 
      ungroup()
    
    # Rename optimized variable
    colnames(optimization_results)[colnames(optimization_results) == i] <- paste0("optimal_", i)
    head(optimization_results)

    out[[i]] <- optimization_results  
  }
  
  # Return outputs
  out$scenarios <- scenarios
  return(out)
  
}


#############################################################################
###      Performing optimization with grid search for all predictors      ###
#############################################################################

opt_file <- sub(".*agg_", "", sub(".*seeds_", "", gp_file))
opt_file <- sub("*_cv.RData", "", opt_file)
print(paste0("Performing grid search and optimization for: ", opt_file))

# Run
result <- GP_grid_search_predictions(gp_file = gp_file,
                                     scale = scale,
                                     ngrid = ngrid,
                                     target_range_size = target_range_size,
                                     param_ranges = param_ranges_cont)

# Write to file
pred <- basename(dirname(gp_file))
saveRDS(result, file = paste0(sim_folder, "gp/GP_grid_optimization/", pred, "/opt_", opt_file, ".rds"))


