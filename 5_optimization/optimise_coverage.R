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
k = as.numeric(args[1])

gp_file = args[1]
ranges_file = args[2]
results_folder = args[3]

# For testing
gp_file = "/scicore/home/penny/GROUP/M3TPP/E0_MAB/gp/trained/prevred_int_y10/seeds_E0MAB_Mali_4.9167_exp_0.1_10_prevred_int_y10_cv.RData"
ranges_file = "/scicore/home/penny/GROUP/M3TPP/E0_MAB/param_ranges.RData"
results_folder = "/scicore/home/penny/GROUP/M3TPP/E0_MAB/optimisation/"


# Load GP model and parameter ranges
gp_result_name = load(gp_file)
gp_result = get(gp_result_name)
rm(gp_result_name)
load(ranges_file)
param_ranges_cont <- param_ranges_cont[-which(param_ranges_cont[,1]== param_ranges_cont[,2]),]

seasonalities <- c("Mali","Sen")
decays <- c("exp","hill","wei")
ages <- c(4.9167)
acesses <- c(0.1,0.5)
EIRs <- c(3,4,5,8,9,20,28,47,150)

EIRs_low_acces <- c(3,4,8,28)
EIRs_high_acces <- c(5,9,20,47)

rows <- seq(1,364,25)


scenarios <- expand.grid(seasonalities,decays,ages,acesses,EIRs,rows)
colnames(scenarios) <- c("Seasonality","Decay","Age","Acess","EIR","row")

index <- which(scenarios$Acess == 0.1 & scenarios$EIR %in% EIRs_high_acces)
scenarios <- scenarios[-index,]

index <- which(scenarios$Acess == 0.5 & scenarios$EIR %in% EIRs_low_acces)
scenarios <- scenarios[-index,]



seasonality <- scenarios[k,1]
decay <- scenarios[k,2]
age <- scenarios[k,3]
acess <- scenarios[k,4]
EIR <- scenarios[k,5]
row_number <- scenarios[k,6]



scen_name <- paste(seasonality,"_",age,"_",decay,"_",acess,"_",EIR, sep = "")
print(scen_name)

filename_UL <- paste('/scicore/home/penny/GROUP/smc_lai/E3E4disc_comparison/gp_trained/gp_trained_UL_', scen_name, '.RData', sep = "")
filename_CI <- paste('/scicore/home/penny/GROUP/smc_lai/E3E4disc_comparison/gp_trained/gp_trained_CI_', scen_name, '.RData', sep = "")
filename_pppy_SMC <- paste('/scicore/home/penny/GROUP/smc_lai/E3E4disc_comparison/gp_trained/gp_trained_pppy_SMC_', scen_name, '.RData', sep = "")
filename_pppy_LAI <- paste('/scicore/home/penny/GROUP/smc_lai/E3E4disc_comparison/gp_trained/gp_trained_pppy_LAI_', scen_name, '.RData', sep = "")

outfile_scenarios <- paste('/scicore/home/penny/GROUP/smc_lai/E3E4disc_comparison/outfiles_optimisation_all/scenarios_',scen_name,"_",row_number,".txt",sep = "")


if (file.exists(filename_UL)){
  print("files exist")
 
  #load GPS

  load(filename_UL)
  load(filename_CI)

  
  SMC_Coverage = seq(0.4,1, by = 0.1)
  Halflife = seq(30, 150, by = 10)
  Efficacy = seq(0.7,1, by = 0.1)
  
  scenarios <- expand.grid(SMC_Coverage,Halflife,Efficacy)
  names(scenarios) <- c("Coverage_SMC","Halflife","Efficacy")
  scenarios$max_diff_non_inf <- 0
  scenarios$optimal_lai_coverage <- 0
  
  
  #change this to selecting appropriate rownames

  for (i in row_number:min(c(row_number + 25,nrow(scenarios)))) {
        
    # Find minimum parameter value for given EIR and prevalence reduction range
    scenario <- scenarios[i, ]
    
    #this is for the optimisation
    get_noninf= function(x) {
      params <- scenario
      param_vec = matrix(unlist(c(x, params[c(1,2,3)])), nrow = 1)
      
      colnames(param_vec) <-  c("Coverage", "Coverage_SMC", "Halflife", "Efficacy")
      
      ul = predict(x = as.matrix(param_vec), GP_trained_UL)$mean
      ci = predict(x = as.matrix(param_vec), GP_trained_CI)$mean
      
      return(ul - ci)
      
    }
    
    #this is to check if optimisation is possible
    max_noninf= function(x) {
      params <- scenario
      param_vec = matrix(unlist(c( x, params[c(1,2,3)])), nrow = 1)
      
      colnames(param_vec) <-  c("Coverage", "Coverage_SMC", "Halflife", "Efficacy")
      
      ul = predict(x = as.matrix(param_vec), GP_trained_UL)$mean
      ci = predict(x = as.matrix(param_vec), GP_trained_CI)$mean
      
      return( - (ul - ci))
    }
    
    
    return_coverage = function(x) {
      return(x)
    }
    
    #check if non inferiority is possible
    
    
    ans_max <- gosolnp(pars  = NULL, fixed = NULL, fun = max_noninf, 
                       LB = -1, UB = 1, distr = rep(1, 1), distr.opt = list(), 
                       n.restarts = 3, control = list(maxit = 100), n.sim = 200)
    
    scenarios[i,]$max_diff_non_inf <- -mean(ans_max$values)
    
    if (-mean(ans_max$values) > 0){
      
      ans_mean = gosolnp(pars  = NULL, fixed = NULL, fun = return_coverage, ineqfun = get_noninf, 
                         ineqLB = 0 , ineqUB = 1, 
                         LB = 0.4, UB = 1, distr = rep(1, 1), distr.opt = list(), 
                         n.restarts = 3, control = list(maxit = 100), n.sim = 200)
      
      scenarios[i,]$optimal_lai_coverage <- ans_mean$pars
      
    } else {
      scenarios[i,]$optimal_lai_coverage <- -1
    }
    
  }
  



write.table(scenarios, file = outfile_scenarios)

}
