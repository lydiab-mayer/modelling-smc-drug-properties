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



#this is to check if optimisation is possible
max_noninf= function(x) {
  params <- scenario
  param_vec = matrix(unlist(c(params[1], x, params[c(2,3,4)])), nrow = 1)
  
  colnames(param_vec) <-  c("EIR", "Coverage", "Coverage_SMC_Res", "Halflife", "Efficacy")
  
  ul = predict(x = as.matrix(param_vec), GP_trained_UL)$mean
  ci = predict(x = as.matrix(param_vec), GP_trained_CI)$mean
  
  return( - (ul - ci))
}


return_coverage = function(x) {
  return(x)
}


# Function returning the predicted mean prevalence reduction for a given input
# Function returning the parameter to be optimized, 
get_param = function(x, GP_model1,GP_model2, param_vec, param_name) {
  return(x)
}


get_noninf = function(x, GP_model1,GP_model2, param_vec, param_name) {
  param_vec[param_name] = x	
  ul = predict(x = param_vec, GP_trained_UL)$mean
  ci = predict(x = param_vec, GP_trained_CI)$mean
  return(ul - ci)
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
gp_file1 = args[1]
gp_file2 = args[2]
ranges_file = args[3]
results_folder = args[4]
opt_setup_file = args[5]
opt_row = args[6]


# # For testing
#  gp_file1 = "/scicore/home/penny/burlyd00/smc_lai/E3E4_comp/gp_5/trained/UL/seeds_E3E4_comp_Mali_4.9167_exp_0.1_UL_cv.RData"
#  gp_file2 = "/scicore/home/penny/burlyd00/smc_lai/E3E4_comp/gp_5/trained/CI_high_HR/seeds_E3E4_comp_Mali_4.9167_exp_0.1_CI_high_HR_cv.RData"
# # 
#  ranges_file = "/scicore/home/penny/burlyd00/smc_lai/E3E4_comp/param_ranges.RData"
#  results_folder = "/scicore/home/penny/burlyd00/smc_lai/E3E4_comp/gp_5/optimisation/non_inf/"
#  opt_setup_file = "/scicore/home/penny/burlyd00/smc_lai/E3E4_comp/gp_5/optimisation/non_inf/opt_setup_file_lhs.txt"
#  opt_row = 1

# Retrieve the optimization specifications
opt_setup_file = "/scicore/home/penny/burlyd00/smc_lai/E3E4_comp/gp_5/optimisation/non_inf/opt_setup_file_lhs.txt"
opt_setup = read.table(opt_setup_file, header=TRUE, sep="\t", stringsAsFactors = FALSE)[opt_row,]

EIRValuesRange = c(5,10,20,50,100,200)

# Load GP model and parameter ranges
gp_result_name2 = load(gp_file2)
gp2 = get(gp_result_name2)
GP_trained_CI =gp2$GP_model
rm(gp_result_name2)

gp_result_name1 = load(gp_file1)
gp1 = get(gp_result_name1)
GP_trained_UL <- gp1$GP_model
rm(gp_result_name1)

load(ranges_file)


# Initialize optimizer settings
opt_var_name = "Coverage"
param_vec = c(as.numeric(opt_setup) )
names(param_vec) <- c("EIR","Coverage","Coverage_SMC_Res" ,"Halflife", "Efficacy")
LB = param_ranges[as.character(opt_var_name), 1]
UB = param_ranges[as.character(opt_var_name), 2]

opt_df = point_df = NULL

for (i in 1:length(EIRValuesRange)) {
  
  
  param_vec["EIR"] = EIRValuesRange[i]

     # Find minimum parameter value for given EIR and prevalence reduction range
    
    # Calculate maximum attainable prevalence reduction
    max_param_vals = param_vec
    max_param_vals[2] = 1

    ul = predict(x = max_param_vals, GP_trained_UL)$mean
    ci = predict(x = max_param_vals, GP_trained_CI)$mean

    # Run optimization algorithm only of k is below max_incred
    ans_mean = ans_sd_plus = ans_sd_minus = ans_sd2_plus = ans_sd2_minus = NULL
    ans_mean$pars = ans_sd_plus$pars = ans_sd_minus$pars = ans_sd2_plus$pars = ans_sd2_minus$pars = NA
    if(ul>ci) {
      result = tryCatch({
          
        ans_mean = gosolnp(pars  = NULL, fixed = NULL, fun = get_param, ineqfun = get_noninf, 
                           ineqLB = 0, ineqUB = c(1), 
                           LB = LB, UB = UB, distr = rep(1, 1), distr.opt = list(), 
                           n.restarts = 5, control = list(maxit = 100), n.sim = 500,  
                           GP_model1 = GP_trained_CI,GP_model2=GP_trained_UL, param_vec = c(param_vec), param_name = opt_var_name)
        

      #   ans_sd_plus = gosolnp(pars  = NULL, fixed = NULL, fun = get_param, ineqfun = get_p_red_sd_plus, 
      #                         ineqLB = c(SMC_Coverage_ranges[k]), ineqUB = c(100), #c(SMC_Coverage_ranges[k+1]), 
      #                         LB = LB, UB = UB, distr = rep(1, length(LB)), distr.opt = list(), 
      #                         n.restarts = 2, control = list(maxit = 100), n.sim = 500,  
      #                         GP_model = gp_result$GP_model, param_vec = t(param_vec), param_name = opt_var_name)
      #   
      #   ans_sd_minus = gosolnp(pars  = NULL, fixed = NULL, fun = get_param, ineqfun = get_p_red_sd_minus, 
      #                          ineqLB = c(SMC_Coverage_ranges[k]), ineqUB = c(100),
      #                          LB = LB, UB = UB, distr = rep(1, length(LB)), distr.opt = list(), 
      #                          n.restarts = 2, control = list(maxit = 100), n.sim = 500,  
      #                          GP_model = gp_result$GP_model, param_vec = t(param_vec), param_name = opt_var_name)
      #   ans_sd2_plus = gosolnp(pars  = NULL, fixed = NULL, fun = get_param, ineqfun = get_p_red_sd2_plus, 
      #                          ineqLB = c(SMC_Coverage_ranges[k]), ineqUB = c(100), #c(SMC_Coverage_ranges[k+1]), 
      #                          LB = LB, UB = UB, distr = rep(1, length(LB)), distr.opt = list(), 
      #                          n.restarts = 2, control = list(maxit = 100), n.sim = 500,  
      #                          GP_model = gp_result$GP_model, param_vec = t(param_vec), param_name = opt_var_name)
      #   
      #   ans_sd2_minus = gosolnp(pars  = NULL, fixed = NULL, fun = get_param, ineqfun = get_p_red_sd2_minus, 
      #                           ineqLB = c(SMC_Coverage_ranges[k]), ineqUB = c(100),
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
    opt_df = rbind.data.frame(opt_df, point_df)
}


# Save the optimisation results
opt_obj = NULL
opt_obj$table = opt_df
opt_obj$opt_setup = param_vec

opt_obj$seasonality = gp1$seasonality
opt_obj$biting_pattern = gp1$biting_pattern
opt_obj$opt_param = opt_var_name
opt_obj$Access = gp1$Access
opt_obj$LAI_dec = gp1$Lai_dec

model_name = file_path_sans_ext(basename(gp_file1))
model_name = str_remove(model_name, "cv")
opt_file = paste(results_folder, model_name, opt_var_name,"_",opt_row, ".RData", sep="")
save(opt_obj, file = opt_file)
