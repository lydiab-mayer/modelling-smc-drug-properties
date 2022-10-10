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


###
# For testing
 gp_file = "/scicore/home/penny/GROUP/M3TPP/E0_MAB/gp/trained/prevred_int_y10/seeds_E0MAB_Mali_4.9167_exp_0.1_10_prevred_int_y10_cv.RData"
 ranges_file = "/scicore/home/penny/GROUP/M3TPP/E0_MAB/param_ranges.RData"
 results_folder = "/scicore/home/penny/GROUP/M3TPP/E0_MAB/optimisation/"
####


# Load GP model and parameter ranges
gp_result_name = load(gp_file)
gp_result = get(gp_result_name)
rm(gp_result_name)
load(ranges_file)

outfile_scenarios <- paste(results_folder,"scenarios_coverage.txt",sep = "")

# Define optimisation problem
# Here: What is the minimal coverage under which a prevalence reduction of at least 0.1 can be achieved?
variable <- "Coverage"
target <- gp_result$predicted
cutoff <- 10
n_gridpoints <- 10
param_ranges <- param_ranges_cont[rownames(param_ranges_cont) != variable,]

grid <- mapply(seq, param_ranges[,1], param_ranges[,2], length.out = n_gridpoints)
scenarios <- expand.grid(grid[,1], grid[,2])

scenarios$optim <- 0
optim_name <- paste("optimal_",variable, sep = "")
colnames(scenarios) <- c(rownames(param_ranges), optim_name)


for (i in 1:nrow(scenarios)) {
        
    scenario <- scenarios[i, ]
    
    get_prevred = function(x) {
      params <- scenario
      param_vec = matrix(unlist(c(x, params[c(1,2)])), nrow = 1)
      
      colnames(param_vec) <-  c("Coverage",  "Halflife", "Efficacy")
      
      prev_red = predict(x = as.matrix(param_vec), gp_result$GP_model$GP_model)$mean

      return(prev_red - cutoff)
      
    }
    
    return_coverage = function(x) {
      return(x)
    }
    
    ans_mean = gosolnp(pars  = NULL, fixed = NULL, fun = return_coverage, ineqfun = get_prevred, 
                         LB = 0.4, UB = 1, distr = rep(1, 1), distr.opt = list(), 
                         n.restarts = 3, control = list(maxit = 100), n.sim = 200, 
                       x = scenario
    )
      
    scenarios[i,]$optim_name <- ans_mean$pars
}
  
write.table(scenarios, file = outfile_scenarios)

