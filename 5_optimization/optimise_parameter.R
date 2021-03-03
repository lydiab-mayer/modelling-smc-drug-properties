library(nloptr)
library(hetGP)
library(Rsolnp)
library(metaheuristicOpt)
library(tools)
library(rapportools)
library(stringr)
library(ggplot2)
library(ggpubr)



args = commandArgs(TRUE)
gp_file = args[1]
ranges_file = args[2]
results_folder = args[3]
opt_setup_file = args[4]
opt_row = args[5]

print(gp_file)
print(ranges_file)

print(opt_setup_file)

print(opt_row)

gp_file = "/scicore/home/penny/GROUP/M3TPP/E0_MAB/gp/trained/prevred_int_y10/seeds_E0MAB_Mali_4.9167_exp_0.1_10_prevred_int_y10_cv.RData"
ranges_file = "/scicore/home/penny/GROUP/M3TPP/E0_MAB/param_ranges.RData"
opt_setup_file = "/scicore/home/penny/GROUP/M3TPP/E0_MAB/gp/optimisation/opt_setup.txt"
results_folder = "/scicore/home/penny/GROUP/M3TPP/E0_MAB/optimisation/"
opt_row=1

# Retrieve the optimization specifications
print(opt_setup_file)
opt_setup = read.table(opt_setup_file, header=TRUE, sep="\t", stringsAsFactors = FALSE)[opt_row,]


# Load GP model and parameter ranges
settings <- strsplit(gp_file, "/", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]][11]


gp_result_name = load(gp_file)
gp_result =  cv_result$GP_model
rm(gp_result_name)
load(ranges_file)
param_ranges <- param_ranges_cont
names_params <- rownames(param_ranges)
# Initialize optimizer settings
opt_var_name = opt_setup$Param_opt
param_vec = opt_setup
param_vec$CM_name = NULL
param_vec$Param_opt = NULL
LB = param_ranges[opt_var_name, 1]
UB = param_ranges[opt_var_name, 2]


# Define optimisation problem
# Here: What is the minimal coverage under which a prevalence reduction of at least 0.1 can be achieved?
variable <- opt_var_name
#target <- gp_result$predicted
cutoff <- opt_setup$reductions*100
prev_ranges
n_gridpoints <- 10
param_ranges <- param_ranges_cont[rownames(param_ranges_cont) != variable,]

grid <- mapply(seq, param_ranges[,1], param_ranges[,2], length.out = n_gridpoints)
scenarios <- expand.grid(grid[,1], grid[,2])

scenarios$optim <- 0
optim_name <- paste("optimal_",variable, sep = "")
colnames(scenarios) <- c(rownames(param_ranges), optim_name)

set_out <- strsplit(settings, ".RData", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]][1]
outfile_scenarios <- paste(results_folder,set_out,"_",opt_var_name, "_opt.RData", sep="")
  

for (i in 1:nrow(scenarios)) {
        
    scenario <- scenarios[i, ]
    
    # Function returning the parameter to be optimized, 
    opt_df = point_df = NULL
    
        print(paste(EIRValuesRange[EIR_lvl], prev_ranges[prev_red_lvl]))
        # Find minimum parameter value for given EIR and prevalence reduction range
        
        # Calculate maximum attainable prevalence reduction
        max_param_vals = t(param_ranges[,2])
          max_param_vals[head(names(scenario),-1),] = head(as.numeric(scenario),-1)
        max_prev_red = predict(x = max_param_vals, gp_result$GP_model)$mean
        
        # Run optimization algorithm only of prev_red_lvl is below max_prev_red
        ans_mean = ans_sd_plus = ans_sd_minus = ans_sd2_plus = ans_sd2_minus = NULL
        ans_mean$pars = ans_sd_plus$pars = ans_sd_minus$pars = ans_sd2_plus$pars = ans_sd2_minus$pars = NA
        
        if(max_prev_red >= prev_ranges[prev_red_lvl] | EIR_lvl < 10) {
          result = tryCatch({
            ans_mean = gosolnp(pars  = NULL, fixed = NULL, fun = get_param, ineqfun = get_p_red, 
                               ineqLB = c(prev_ranges[prev_red_lvl]), ineqUB = c(100), 
                               LB = LB, UB = UB, distr = rep(1, length(LB)), distr.opt = list(), 
                               n.restarts = 10, control = list(maxit = 100), n.sim = 1000,  
                               GP_model = gp_result$GP_model, param_vec = t(param_vec), param_name = opt_var_name)
            
            ans_sd_plus = gosolnp(pars  = NULL, fixed = NULL, fun = get_param, ineqfun = get_p_red_sd_plus, 
                                  ineqLB = c(prev_ranges[prev_red_lvl]), ineqUB = c(100), #c(prev_ranges[prev_red_lvl+1]), 
                                  LB = LB, UB = UB, distr = rep(1, length(LB)), distr.opt = list(), 
                                  n.restarts = 2, control = list(maxit = 100), n.sim = 500,  
                                  GP_model = gp_result$GP_model, param_vec = t(param_vec), param_name = opt_var_name)
            
            ans_sd_minus = gosolnp(pars  = NULL, fixed = NULL, fun = get_param, ineqfun = get_p_red_sd_minus, 
                                   ineqLB = c(prev_ranges[prev_red_lvl]), ineqUB = c(100),
                                   LB = LB, UB = UB, distr = rep(1, length(LB)), distr.opt = list(), 
                                   n.restarts = 2, control = list(maxit = 100), n.sim = 500,  
                                   GP_model = gp_result$GP_model, param_vec = t(param_vec), param_name = opt_var_name)
            ans_sd2_plus = gosolnp(pars  = NULL, fixed = NULL, fun = get_param, ineqfun = get_p_red_sd2_plus, 
                                   ineqLB = c(prev_ranges[prev_red_lvl]), ineqUB = c(100), #c(prev_ranges[prev_red_lvl+1]), 
                                   LB = LB, UB = UB, distr = rep(1, length(LB)), distr.opt = list(), 
                                   n.restarts = 2, control = list(maxit = 100), n.sim = 500,  
                                   GP_model = gp_result$GP_model, param_vec = t(param_vec), param_name = opt_var_name)
            
            ans_sd2_minus = gosolnp(pars  = NULL, fixed = NULL, fun = get_param, ineqfun = get_p_red_sd2_minus, 
                                    ineqLB = c(prev_ranges[prev_red_lvl]), ineqUB = c(100),
                                    LB = LB, UB = UB, distr = rep(1, length(LB)), distr.opt = list(), 
                                    n.restarts = 2, control = list(maxit = 100), n.sim = 500,  
                                    GP_model = gp_result$GP_model, param_vec = t(param_vec), param_name = opt_var_name)
          }, warning = function(w) {
            
          }, error = function(e) {
            print(paste("An error was handled:", e))
          }, finally = {
            
          })
        } else {
          print("Maximum attainable prevalence reduction is lower than desired level.")
        }
        
    scenarios[i,]$optim_name <- ans_mean$pars
}
  
write.table(scenarios, file = outfile_scenarios)

