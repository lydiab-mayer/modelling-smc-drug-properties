###########################
# Solve the following optimization problem: given fixed efficacy, 
# halflife, EIR and SMC coverage, what is the minimum LAI coveragae that 
# achieves a 0.8 probability of non-inferiority
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


gp_file = "/scicore/home/smith/GROUP/smc_lai/E3_E4_comparison/gp_5/trained/agg_E3_E4_comparison_seasonal_Low_indoor_exp_cv.RData"
results_folder = "/scicore/home/smith/laagmi01/smc_lai/analysis_workflow/analysis_scripts/outfiles"

# Load GP model 
load(gp_file)
GP_trained <- cv_result$GP_model



# Define prevalence and EIR ranges
EIRValuesRange = seq(1, 25, by = 5)
SMC_Coverage = seq(0.4,1, by= 0.1)
Halflife = seq(30, 150, by = 20)
Efficacy = seq(0.7,1, by = 0.1)



scenarios <- expand.grid(EIRValuesRange,SMC_Coverage,Halflife,Efficacy)
names(scenarios) <- c("EIR","Coverage_SMC","Halflife","Efficacy")
scenarios$max_prob_non_inf <- 0
scenarios$optimal_lai_coverage <- 0


for (i in 1:nrow(scenarios)) {
  
  scenario <- scenarios[i, ]
  
  # Function describing the relationship between the LAI coverage
  # and the probability of non-inferiority for a fixed
  # EIR, Halflife, Efficacy and SMC coverage
  
  get_prob_noninf= function(x) {
    params <- scenario
    GP_model <- GP_trained
    param_vec = matrix(unlist(c(params[1], x, params[c(2,3,4)])), nrow = 1)
    
    colnames(param_vec) <-  c("EIR", "Coverage", "Coverage_SMC", "Halflife", "Efficacy")
    prob_noninf = predict(x = as.matrix(param_vec), GP_model)$mean
    
    return(prob_noninf)
  }
  
  
  max_prob_noninf= function(x) {
    params <- scenario
    GP_model <- GP_trained
    param_vec = matrix(unlist(c(params[1], x, params[c(2,3,4)])), nrow = 1)
    
    colnames(param_vec) <-  c("EIR", "Coverage", "Coverage_SMC", "Halflife", "Efficacy")
    prob_noninf = predict(x = as.matrix(param_vec), GP_model)$mean
    
    return(-prob_noninf)
  }
  
  return_coverage = function(x) {
    return(x)
  }
  
  #check if the maximum of the probability of non inferiority is greater than 0.8
  ans_max <- gosolnp(pars  = NULL, fixed = NULL, fun = max_prob_noninf, 
                     LB = 0, UB = 1, distr = rep(1, 1), distr.opt = list(), 
                     n.restarts = 2, control = list(maxit = 100), n.sim = 100)
  
  scenarios[i,]$max_prob_non_inf <- -mean(ans_max$values)
  
  if (-mean(ans_max$values) > 0.8){
    
    ans_mean = gosolnp(pars  = NULL, fixed = NULL, fun = return_coverage, ineqfun = get_prob_noninf, 
                       ineqLB = 0.6 , ineqUB = 1, 
                       LB = 0, UB = 1, distr = rep(1, 1), distr.opt = list(), 
                       n.restarts = 2, control = list(maxit = 100), n.sim = 100)
    
    scenarios[i,]$optimal_lai_coverage <- ans_mean$pars
    
  }
  
}






#Save the optimisation results
opt_res <- list()
opt_res$seasonality = GP_trained$seasonality
opt_res$biting_pattern = GP_trained$biting_pattern
opt_res$LAI_decay = GP_trained$LAI_decay

opt_res$scenarios = scenarios

opt_file = paste(results_folder, "/test_ml", ".RData", sep="")
save(opt_res, file = opt_file)


