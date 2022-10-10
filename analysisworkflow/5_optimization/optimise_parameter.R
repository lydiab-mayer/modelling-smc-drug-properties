#####################################
#####################################
###                               ###
###     STEP 5: OPTIMISATION      ###
###                               ###
#####################################
#####################################


### -------------------------------------------------------------------------
###
### M3TPP PROJECT:
### Supporting script for performing optimisation step of workflow 
### 
### Original script:
### Created ???
### Monica Golumbeanu 
###
### Adapted script:
### Saved September 2021
### lydia.braunack-mayer@swisstph.ch
###
### R version 3.6.0
###
### -------------------------------------------------------------------------


##############
### HEADER ###
##############

# Load packages
library(hetGP)
library(Rsolnp)

# Packages required by previous versions of code - stored here against future need
#library(nloptr)
#library(metaheuristicOpt)
#library(tools)
#library(rapportools)
#library(stringr)
#library(ggplot2)
#library(ggpubr)

# Load required custom functions
source("../../../analysisworkflow/5_optimization/optimisation_resources.R")

# Set arguments
args <- commandArgs(TRUE)
gp_file <- args[1]
ranges_file <- args[2]
results_folder <- args[3]
opt_setup_file <- args[4]
opt_row <- 1
n_gridpoints <- as.numeric(args[5])
scale <- as.logical(args[6])

# Sample arguments - stored here to facilitate testing
# gp_file <- "/scicore/home/penny/GROUP/M3TPP/E0_LAIExampleLBM/gp/trained/inc_red_int/seeds_E0LAIExampleLBM_wideseasonal_Mali_10_4.9167_exp_0.1_inc_red_int_cv.RData"
# ranges_file <- "/scicore/home/penny/GROUP/M3TPP/E0_LAIExampleLBM/param_ranges.RData"
# results_folder <- "/scicore/home/penny/GROUP/M3TPP/E0_LAIExampleLBM/gp/optimisation/"
# opt_setup_file <- "/scicore/home/penny/GROUP/M3TPP/E0_LAIExampleLBM/gp/optimisation/inc_red_int/Halflife_10_opt_setup.txt"
# opt_row <- 1
# n_gridpoints <- as.numeric("3")
# scale <- as.logical("TRUE")


##################
###   SET UP   ###
##################

# Retrieve the optimisation specifications
print(gp_file)
print(opt_setup_file)
opt_setup = read.table(opt_setup_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)[opt_row,]

# Load GP model and parameter ranges
settings <- strsplit(gp_file, "/", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]][11]
gp_result_name = load(gp_file)
gp_result =  cv_result$GP_model
rm(gp_result_name)
load(ranges_file)

# Initialize optimiser settings
param_ranges <- param_ranges_cont
names_params <- rownames(param_ranges)
opt_var_name <- opt_setup$Param_opt
param_vec <- opt_setup
param_vec$CM_name <- NULL
param_vec$Param_opt <- NULL
n_restarts <- 10
max_it <- 150
n_sim <- 500

# If inputs have been scaled to c(0, 1) when training the emulator, update parameter ranges to match
if (scale) {
  param_ranges[, 1] <- 0 
  param_ranges[, 2] <- 1 
}

LB <- param_ranges[opt_var_name, 1]
UB <- param_ranges[opt_var_name, 2]

# Set up optimisation problem and file to store results
variable <- opt_var_name
cutoff <- opt_setup$reductions

grid <- mapply(seq, param_ranges[rownames(param_ranges) != variable, 1], param_ranges[rownames(param_ranges) != variable, 2], length.out = n_gridpoints)
scenarios <- expand.grid(grid[, 1], grid[, 2])

scenarios$optim <- NA
optim_name <- paste("optimal_", variable, sep = "")
colnames(scenarios) <- c(rownames(param_ranges)[rownames(param_ranges) != variable], optim_name)

set_out <- strsplit(settings, ".RData", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]][1]
outfile_scenarios <- paste(results_folder, set_out, "_", opt_var_name, "_cutoff", cutoff, "_opt.txt", sep = "")

print("Optimisation setup:")
print(scenarios)


##################
## OPTIMISATION ##
##################

# For each combination of parameters, run optimisation problem
for (i in 1:nrow(scenarios)) {
  
  # Initialise problem
  scenario <- scenarios[i, ]
  order <- names_params
  ans_mean <- ans_sd_plus <- ans_sd_minus <- ans_sd2_plus <- ans_sd2_minus <- NULL
  ans_mean$pars <- ans_sd_plus$pars <- ans_sd_minus$pars <- ans_sd2_plus$pars <- ans_sd2_minus$pars <- NA
  print("Run optimisation for:")
  print(scenario)
  
  # Calculate maximum attainable prevalence reduction
  max_param_vals <- scenario[, 1:ncol(scenario) - 1]
  max_param_vals[, opt_var_name] <- param_ranges[opt_var_name, 2]
  
  max_param_vals <- max_param_vals[, order]
  max_red <- predict(x = as.matrix(max_param_vals), gp_result)$mean

  # Run optimisation algorithm only if target reduction is below maximum achievable reduction     
  if(max_red >= cutoff)  {
    result <- tryCatch({
      
      # Optimize for target
      print("Optimising target")
      ans_mean <- gosolnp(pars = NULL, fixed = NULL, fun = get_param, ineqfun = get_p_red,
                          ineqLB = cutoff, ineqUB = c(100), LB = LB, UB = UB, 
                          distr = rep(1, length(LB)), distr.opt = list(), n.restarts = n_restarts, 
                          control = list(maxit = max_it), n.sim = n_sim, GP_model = gp_result, 
                          param_vec = as.vector(max_param_vals), param_name = opt_var_name)

      # Optimize for 1 standard deviation above target 
      print("Optimising one standard deviation above target")
      ans_sd_plus <- gosolnp(pars  = NULL, fixed = NULL, fun = get_param, ineqfun = get_p_red_sd_plus,
                             ineqLB = cutoff, ineqUB = c(100), LB = LB, UB = UB, 
                             distr = rep(1, length(LB)), distr.opt = list(), n.restarts = n_restarts, 
                             control = list(maxit = max_it), n.sim = n_sim, GP_model = gp_result, 
                             param_vec = as.vector(max_param_vals), param_name = opt_var_name)
      
      # Optimize for 1 standard deviation below target 
      print("Optimising one standard deviation below target")
      ans_sd_minus <- gosolnp(pars  = NULL, fixed = NULL, fun = get_param, ineqfun = get_p_red_sd_minus,
                              ineqLB = cutoff, ineqUB = c(100), LB = LB, UB = UB, 
                              distr = rep(1, length(LB)), distr.opt = list(), n.restarts = n_restarts, 
                              control = list(maxit = max_it), n.sim = n_sim, GP_model = gp_result, 
                              param_vec = as.vector(max_param_vals), param_name = opt_var_name)
      
      # Optimize for 2 standard deviations above target 
      print("Optimising two standards deviation above target")
      ans_sd2_plus <- gosolnp(pars  = NULL, fixed = NULL, fun = get_param, ineqfun = get_p_red_sd2_plus,
                              ineqLB = cutoff, ineqUB = c(100), LB = LB, UB = UB, 
                              distr = rep(1, length(LB)), distr.opt = list(), n.restarts = n_restarts, 
                              control = list(maxit = max_it), n.sim = n_sim, GP_model = gp_result, 
                              param_vec = as.vector(max_param_vals), param_name = opt_var_name)
      
      # Optimize for 2 standard deviations below target 
      print("Optimising two standard deviations below target")
      ans_sd2_minus <- gosolnp(pars  = NULL, fixed = NULL, fun = get_param, ineqfun = get_p_red_sd2_minus,
                               ineqLB = cutoff, ineqUB = c(100), LB = LB, UB = UB, 
                               distr = rep(1, length(LB)), distr.opt = list(), n.restarts = n_restarts, 
                               control = list(maxit = max_it), n.sim = n_sim, GP_model = gp_result, 
                               param_vec = as.vector(max_param_vals), param_name = opt_var_name)
      
      }, warning = function(w) {
        
      }, error = function(e) {
        print(paste("An error was handled:", e))
        }, finally = {
            
          })
    } else {
      print("Maximum attainable prevalence reduction is lower than desired level.")
      }
     
    # Store outputs 
    scenarios[i, which(colnames(scenarios) == optim_name)] <- ans_mean$pars
    scenarios[i, "sd_plus"] <- ans_sd_plus$pars
    scenarios[i, "sd_minus"] <- ans_sd_minus$pars
    scenarios[i, "sd2_plus"] <- ans_sd2_plus$pars
    scenarios[i, "sd2_minus"] <- ans_sd2_minus$pars
    
}

# If inputs have been scaled to c(0, 1), transform back to their original scale
if (scale) {
  
  # Transform parameters that have not been optimised
  for (i in names_params[names_params != opt_var_name]) {
      scenarios[, i] <- scenarios[, i] * (param_ranges_cont[i, 2] - param_ranges_cont[i, 1]) + param_ranges_cont[i, 1]
  }
  
  # Transform optimised parameter
  scenarios[, optim_name] <- scenarios[, optim_name] * (param_ranges_cont[opt_var_name, 2] - param_ranges_cont[opt_var_name, 1]) + param_ranges_cont[opt_var_name, 1]
  scenarios[, "sd_plus"] <- scenarios[, "sd_plus"] * (param_ranges_cont[opt_var_name, 2] - param_ranges_cont[opt_var_name, 1]) + param_ranges_cont[opt_var_name, 1]
  scenarios[, "sd_minus"] <- scenarios[, "sd_minus"] * (param_ranges_cont[opt_var_name, 2] - param_ranges_cont[opt_var_name, 1]) + param_ranges_cont[opt_var_name, 1]
  scenarios[, "sd2_plus"] <- scenarios[, "sd2_plus"] * (param_ranges_cont[opt_var_name, 2] - param_ranges_cont[opt_var_name, 1]) + param_ranges_cont[opt_var_name, 1]
  scenarios[, "sd2_minus"] <- scenarios[, "sd2_minus"] * (param_ranges_cont[opt_var_name, 2] - param_ranges_cont[opt_var_name, 1]) + param_ranges_cont[opt_var_name, 1]

}


##################
##    OUTPUTS   ##
##################

# Write outputs to table
write.table(scenarios, file = outfile_scenarios)

