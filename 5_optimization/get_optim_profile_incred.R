###########################
# Solve the following optimization problem: given fixed (default = average) efficacy, 
# halflife and access to treatment, at what minimum coverage should the intervention
# be implemented to achieve a given prevalence reduction level at a gven transmission level
# 
# created 08.01.2019
# monica.golumbeanu@unibas.ch
###########################
library(nloptr)
library(hetGP)
library(Rsolnp)
library(metaheuristicOpt)
library(tools)
library(rapportools)
library(stringr)



# Function returning the predicted mean prevalence reduction for a given input
# Function returning the parameter to be optimized, 
get_param = function(x, GP_model, param_vec, param_name) {
  return(x)
}


get_p_red = function(x, GP_model, param_vec, param_name) {
  param_vec[param_name] = x
  incred = predict(x = c(param_vec), GP_model)$mean
  if(incred<0) {
    incred = 0
  }
  return(incred)
} 

# Function returning the predicted mean prevalence reduction - sd for a given input
get_p_red_sd_minus = function(x, GP_model, param_vec, param_name) {
  param_vec[,param_name] = x
  prediction_res = predict(x = param_vec, GP_model)
  incred = prediction_res$mean - sqrt(prediction_res$sd2 + prediction_res$nugs)
  if(incred<0) {
    incred = 0
  }
  return(incred)
}

# Function returning the predicted mean prevalence reduction + sd for a given input
get_p_red_sd_plus = function(x, GP_model, param_vec, param_name) {
  param_vec[param_name,] = x
  prediction_res = predict(x = t(param_vec), GP_model)
  incred = prediction_res$mean + sqrt(prediction_res$sd2 + prediction_res$nugs)
  if(incred<0) {
    incred = 0
  }
  return(incred)
}

# Function returning the predicted mean prevalence reduction - 2sd for a given input
get_p_red_sd2_minus = function(x, GP_model, param_vec, param_name) {
  param_vec[param_name,] = x
  prediction_res = predict(x = t(param_vec), GP_model)
  incred = prediction_res$mean - sqrt(2*prediction_res$sd2 + prediction_res$nugs)
  if(incred<0) {
    incred = 0
  }
  return(incred)
}

# Function returning the predicted mean prevalence reduction + 2sd for a given input
get_p_red_sd2_plus = function(x, GP_model, param_vec, param_name) {
  param_vec[param_name,] = x
  prediction_res = predict(x = t(param_vec), GP_model)
  incred = prediction_res$mean + sqrt(2*prediction_res$sd2 + prediction_res$nugs)
  if(incred<0) {
    incred = 0
  }
  return(incred)
}



args = commandArgs(TRUE)
gp_file = args[1]
ranges_file = args[2]
results_folder = args[3]
opt_setup_file = args[4]
opt_row = args[5]


# For testing
#gp_file = "/scicore/home/smith/GROUP/smc_lai/E2_LAI/gp_5/trained/incred/seeds_E2_LAI_Mali_4.9167_exp_0.1_incred_cv.RData"
#ranges_file = "/scicore/home/smith/GROUP/smc_lai/E2_LAI/param_ranges.RData"
#results_folder = "/scicore/home/smith/GROUP/smc_lai/E2_LAI/gp_5/optimisation/incred/"
#opt_setup_file = "/scicore/home/smith/GROUP/smc_lai/E2_LAI/gp_5/optimisation/incred/opt_setup_file.txt"
#opt_row = 1

# Retrieve the optimization specifications
opt_setup = read.table(opt_setup_file, header=TRUE, sep="\t", stringsAsFactors = FALSE)[opt_row,]

# Load GP model and parameter ranges
gp_result_name = load(gp_file)
gp_result = get(gp_result_name)
rm(gp_result_name)
load(ranges_file)

# Define incidence and EIR ranges
EIRValuesRange = c(1,5,10,20,50,100,150,200)
inc_ranges = seq(0.1,1,0.1)

# For testing only:
# EIRValuesRange = seq(1, 3, by = 1)
# inc_ranges = c(10, 20)

# Initialize optimizer settings
opt_var_name = "Coverage"
param_vec = c(0,0,as.numeric(opt_setup) )
names(param_vec) <- c("EIR" ,     "Coverage" ,"Halflife", "Efficacy")
LB = param_ranges[as.character(opt_var_name), 1]
UB = param_ranges[as.character(opt_var_name), 2]

opt_df = point_df = NULL

for (i in 1:length(EIRValuesRange)) {
  
 
  param_vec["EIR"] = EIRValuesRange[i]
  for (k in 1:length(inc_ranges)) {
    print(paste(EIRValuesRange[i], inc_ranges[k]))
    # Find minimum parameter value for given EIR and prevalence reduction range
    
    # Calculate maximum attainable prevalence reduction
    max_param_vals = param_vec
    max_param_vals[2] = 1
    max_incred = predict(x = max_param_vals, gp_result$GP_model)$mean
    
    # Run optimization algorithm only of k is below max_incred
    ans_mean = ans_sd_plus = ans_sd_minus = ans_sd2_plus = ans_sd2_minus = NULL
    ans_mean$pars = ans_sd_plus$pars = ans_sd_minus$pars = ans_sd2_plus$pars = ans_sd2_minus$pars = NA
    if(max_incred >= inc_ranges[k] ) {
      result = tryCatch({
        ans_mean = gosolnp(pars  = NULL, fixed = NULL, fun = get_param, ineqfun = get_p_red, 
                           ineqLB = c(inc_ranges[k]), ineqUB = c(1), 
                           LB = LB, UB = UB, distr = rep(1, length(LB)), distr.opt = list(), 
                           n.restarts = 5, control = list(maxit = 100), n.sim = 500,  
                           GP_model = gp_result$GP_model, param_vec = c(param_vec), param_name = opt_var_name)
      
      
        
          
      #   ans_sd_plus = gosolnp(pars  = NULL, fixed = NULL, fun = get_param, ineqfun = get_p_red_sd_plus, 
      #                         ineqLB = c(inc_ranges[k]), ineqUB = c(100), #c(inc_ranges[k+1]), 
      #                         LB = LB, UB = UB, distr = rep(1, length(LB)), distr.opt = list(), 
      #                         n.restarts = 2, control = list(maxit = 100), n.sim = 500,  
      #                         GP_model = gp_result$GP_model, param_vec = t(param_vec), param_name = opt_var_name)
      #   
      #   ans_sd_minus = gosolnp(pars  = NULL, fixed = NULL, fun = get_param, ineqfun = get_p_red_sd_minus, 
      #                          ineqLB = c(inc_ranges[k]), ineqUB = c(100),
      #                          LB = LB, UB = UB, distr = rep(1, length(LB)), distr.opt = list(), 
      #                          n.restarts = 2, control = list(maxit = 100), n.sim = 500,  
      #                          GP_model = gp_result$GP_model, param_vec = t(param_vec), param_name = opt_var_name)
      #   ans_sd2_plus = gosolnp(pars  = NULL, fixed = NULL, fun = get_param, ineqfun = get_p_red_sd2_plus, 
      #                          ineqLB = c(inc_ranges[k]), ineqUB = c(100), #c(inc_ranges[k+1]), 
      #                          LB = LB, UB = UB, distr = rep(1, length(LB)), distr.opt = list(), 
      #                          n.restarts = 2, control = list(maxit = 100), n.sim = 500,  
      #                          GP_model = gp_result$GP_model, param_vec = t(param_vec), param_name = opt_var_name)
      #   
      #   ans_sd2_minus = gosolnp(pars  = NULL, fixed = NULL, fun = get_param, ineqfun = get_p_red_sd2_minus, 
      #                           ineqLB = c(inc_ranges[k]), ineqUB = c(100),
      #                           LB = LB, UB = UB, distr = rep(1, length(LB)), distr.opt = list(), 
      #                           n.restarts = 2, control = list(maxit = 100), n.sim = 500,  
      #                           GP_model = gp_result$GP_model, param_vec = t(param_vec), param_name = opt_var_name)
      }, warning = function(w) {
         
      }, error = function(e) {
        print(paste("An error was handled:", e))
     }, finally = {
         
       })
        point_df$opt_param = ans_mean$pars
        
    } else {
      print("Maximum attainable prevalence reduction is lower than desired level.")
      point_df$opt_param = 0
      
    }
    # Update the results data frame
    point_df$EIR = EIRValuesRange[i]
    point_df$incred = inc_ranges[k]
    
  
    opt_df = rbind.data.frame(opt_df, point_df)
  }
}
# Save the optimisation results
opt_obj = NULL
opt_obj$table = opt_df
opt_obj$opt_setup = opt_setup

opt_obj$seasonality = gp_result$seasonality
opt_obj$biting_pattern = gp_result$biting_pattern
opt_obj$opt_param = opt_var_name
opt_obj$Access = gp_result$Access
opt_obj$LAI_dec = gp_result$Lai_dec

model_name = file_path_sans_ext(basename(gp_file))
model_name = str_remove(model_name, "cv")
opt_file = paste(results_folder, model_name, opt_var_name,"_",opt_row, ".RData", sep="")
save(opt_obj, file = opt_file)
