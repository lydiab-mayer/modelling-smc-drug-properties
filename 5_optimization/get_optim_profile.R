###########################
# Solve the following optimization problem: given fixed (default = average) efficacy, 
# halflife and access to treatment, at what minimum coverage should the intervention
# be implemented to achieve a given prevalence reduction level at a gven transmission level
# 
# created 08.01.2019
# monica.golumbeanu@unibas.ch
###########################

# This document is still in testing, remove testing set-up before running on cluster! 
# testing specified:
# limited parameter ranges for minimisation
# input files
# decreased n.sim, normally 500

# to be added:
# usage in workflow with changed inputs from bash
# saving results, GP needs to have LAI decay and CM saved in list 
################
library(nloptr)
library(hetGP)
library(Rsolnp)
library(metaheuristicOpt)
library(tools)
library(rapportools)
library(stringr)

<<<<<<< HEAD
source('supp/optimisation_resources.R')
=======
source('./supp/optimisation_resources.R')
>>>>>>> 0729b2856e92fc7694c93c8f460b4b7416f2f0b6
args = commandArgs(TRUE)
gp_file = args[1]
ranges_file = args[2]
results_folder = args[3]


# For testing
gp_file = "/scicore/home/smith/GROUP/smc_lai/E3_E4_comparison/gp_5/trained/agg_E3_E4_comparison_seasonal_Low_indoor_exp_cv.RData"
ranges_file = "/scicore/home/smith/GROUP/smc_lai/E3_E4_comparison/param_ranges.RData"
results_folder = "/scicore/home/smith/GROUP/smc_lai/E3_E4_comparison/gp_5/optimisation/"
#opt_setup_file = "/scicore/home/smith/GROUP/smc_lai/E3_E4_comparison/opt_setup.txt"


# Load GP model and parameter ranges
gp_result_name = load(gp_file)
gp_result = get(gp_result_name)
rm(gp_result_name)
load(ranges_file)
param_ranges <- param_ranges[-which(param_ranges[,1]== param_ranges[,2]),]

# Define prevalence and EIR ranges
#EIRValuesRange = seq(1, 25, by = 1)
#SMC_Coverage = seq(0.4,1, by=1)
#Halflife = seq(30, 150, by = 10)
#Efficacy = seq(0.7,1, by=5)


# For testing only:
EIRValuesRange = c(1, 25)
 SMC_Coverage = c(0.6, 0.8)
 Halflife = c(30, 150)
 Efficacy = c(0.7,1)
Coverage= c(1)

scenarios <- expand.grid(EIRValuesRange,Coverage,SMC_Coverage,Halflife,Efficacy)
names(scenarios) <- c("EIR","Coverage","Coverage_SMC","Halflife","Efficacy")



opt_var_name="Coverage"

LB=param_ranges["Coverage", 1] # if variable to be optimised is read in by file: #LB = param_ranges[opt_var_name, 1]

UB=param_ranges["Coverage", 2]


scenario_res <- list()
for (i in 1:nrow(scenarios )) {

   # Find minimum parameter value for given EIR and prevalence reduction range
    scenario <- scenarios[i, ]
    # Calculate maximum attainable prevalence reduction
    max_param_vals = scenario
    max_prob_noninferiority = predict(x = as.matrix(max_param_vals), gp_result$GP_model)$mean
    
    
    # Run optimization algorithm only if scenario is able to show non-inferiority with best parameter setting
    
    ans_mean = ans_sd_plus = ans_sd_minus = ans_sd2_plus = ans_sd2_minus = NULL
    ans_mean$pars = ans_sd_plus$pars = ans_sd_minus$pars = ans_sd2_plus$pars = ans_sd2_minus$pars = NA
   
     if(max_prob_noninferiority >= 0.8 ) {

     #test
         ans = gosolnp(pars  = NULL, fixed = NULL, fun = get_param, eqfun = get_prob_noninf, eqB = c(0.8), 
                       LB = LB, UB = UB, distr = rep(1, length(LB)), distr.opt = list(), 
                       n.restarts = 2, control = list(maxit = 100), n.sim = 400, 
                     GP_model = gp_result$GP_model, param_vec = scenario, param_name = opt_var_name)


         param_vec[which(names(param_vec)==param_name)] = 0.70742
         prob_noninf = predict(x = as.matrix(param_vec), GP_model)$mean
         prob_noninf
            
        ans_sd_plus = gosolnp(pars  = NULL, fixed = NULL, fun = get_prob_noninf_sd_plus, eqfun = get_prob_noninf_sd_plus,  eqB = c(0.9),
                              LB = LB, UB = UB, distr = rep(1, length(LB)), distr.opt = list(), 
                              n.restarts = 2, control = list(maxit = 100), n.sim = 100,  
                              GP_model = gp_result$GP_model, param_vec = scenario, param_name = opt_var_name)
        
        ans_sd_minus = gosolnp(pars  = NULL, fixed = NULL, fun = get_prob_noninf_sd_minus, eqfun = get_prob_noninf_sd_minus,  eqB = c(0.9),
                               LB = LB, UB = UB, distr = rep(1, length(LB)), distr.opt = list(), 
                               n.restarts = 2, control = list(maxit = 100), n.sim = 100,  
                               GP_model = gp_result$GP_model, param_vec = scenario, param_name = opt_var_name)
       
         ans_sd2_plus = gosolnp(pars  = NULL, fixed = NULL, fun = get_prob_noninf_sd2_plus, eqfun = get_prob_noninf_sd2_plus,  eqB = c(0.9),
                     LB = LB, UB = UB, distr = rep(1, length(LB)), distr.opt = list(), 
                     n.restarts = 2, control = list(maxit = 100), n.sim = 100,  
                     GP_model = gp_result$GP_model, param_vec = scenario, param_name = opt_var_name)
        
        ans_sd2_minus = gosolnp(pars  = NULL, fixed = NULL, fun = get_prob_noninf_sd2_minus, eqfun = get_prob_noninf_sd2_minus, eqB = c(0.9), 
                                LB = LB, UB = UB, distr = rep(1, length(LB)), distr.opt = list(), 
                                n.restarts = 2, control = list(maxit = 100), n.sim = 100,  
                                GP_model = gp_result$GP_model, param_vec = scenario, param_name = opt_var_name)
      

     }
    
    scenario[opt_var_name] <-NULL
    # Update the results data frame
    point_res <- c(scenario)
    point_res$opt_param = ans_mean$pars
    point_res$opt_param_sd_plus = ans_sd_plus$pars
    point_res$opt_param_sd_minus = ans_sd_minus$pars
    point_res$opt_param_sd2_plus = ans_sd2_plus$pars
    point_res$opt_param_sd2_minus = ans_sd2_minus$pars
    
scenario_res[[i]] <- unlist(point_res)
      }

opt_res <- bind_rows(scenario_res)


#Save the optimisation results
opt_res$seasonality = gp_result$seasonality
opt_res$biting_pattern = gp_result$biting_pattern
opt_res$LAI_decay = gp_result$LAI_decay

opt_obj$opt_param = opt_var_name
opt_obj$CM_level = opt_setup$CM_name
model_name = file_path_sans_ext(basename(gp_file))
model_name = str_remove(model_name, "cv_as")
opt_file = paste(results_folder, model_name, opt_setup$CM_name, "_", opt_setup$Param_opt, ".RData", sep="")
save(opt_obj, file = opt_file)

# Testing plot:
# ggplot(opt_df, aes(EIR, prev_red, fill= opt_param)) + geom_tile() + scale_fill_gradient(low = "#c994c7", high = "#756bb1")

