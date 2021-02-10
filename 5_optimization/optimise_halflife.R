library(nloptr)
library(hetGP)
library(Rsolnp)
library(metaheuristicOpt)
library(tools)
library(rapportools)
library(stringr)
library(ggplot2)
library(ggpubr)



seasonalities <- c("Mali","Sen")
decays <- c("exp","hill","wei")
ages <- c(4.9167)
acesses <- c(0.1)
EIRs_low_acces <- c(3,4,8,28,150)

cpppy_low_access <- c(0.41,
                      0.71,
                      1.4,
                      2.4,
                      3.3)

EIRs <- c(EIRs_low_acces)




settings <- expand.grid(EIRs,seasonalities,decays,ages,acesses)
colnames(settings) <- c("EIR","Seasonality","Decay","Age","Acess")


scens <- c(1,2,3,4,5,11,12,13,14,15)

for(k in seq(1,length(scens)) ) {

seasonality <- settings[scens[k],2]
decay <- settings[scens[k],3]
age <- settings[scens[k],4]
acess <- settings[scens[k],5]
EIR <- settings[scens[k],1]



scen_name <- paste(seasonality,"_",age,"_",decay,"_",acess,"_",EIR, sep = "")
print(scen_name)

filename_UL <- paste('/scicore/home/penny/GROUP/smc_lai/E3E4disc_comparison/gp_trained/gp_trained_UL_', scen_name, '.RData', sep = "")
filename_CI <- paste('/scicore/home/penny/GROUP/smc_lai/E3E4disc_comparison/gp_trained/gp_trained_CI_', scen_name, '.RData', sep = "")
filename_pppy_SMC <- paste('/scicore/home/penny/GROUP/smc_lai/E3E4disc_comparison/gp_trained/gp_trained_pppy_SMC_', scen_name, '.RData', sep = "")
filename_pppy_LAI <- paste('/scicore/home/penny/GROUP/smc_lai/E3E4disc_comparison/gp_trained/gp_trained_pppy_LAI_', scen_name, '.RData', sep = "")

outfile_scenarios <- paste('/scicore/home/penny/GROUP/smc_lai/E3E4_lb/outfiles_optimisation_HL/scenarios_',scen_name,"test.txt",sep = "")


if (file.exists(filename_UL)){
  print("files exist")
 
  #load GPS

  load(filename_pppy_LAI)
  load(filename_pppy_SMC)
  
  
  load(filename_UL)
  load(filename_CI)

  SMC_Coverage = 0.6
  Coverage = seq(0.4,1, by = 0.2)
  Efficacy = seq(0.7,1, by = 0.1)
  
  scenarios <- expand.grid(Coverage,SMC_Coverage,Efficacy)
  names(scenarios) <- c("Coverage","Coverage_SMC","Efficacy")
  scenarios$max_diff_non_inf <- 0
  scenarios$optimal_lai_HL <- 0
  scenarios$inc_red <- 0
  scenarios$add_inc_red <- 0
  
  
  #change this to selecting appropriate rownames

  for (i in seq(1,nrow(scenarios)) ) {
        
    # Find minimum parameter value for given EIR and prevalence reduction range
    scenario <- scenarios[i, ]
    
    #this is for the optimisation
    get_noninf= function(x) {
      params <- scenario
      param_vec = matrix(unlist(c(params[c(1,2)],x, params[c(3)])), nrow = 1)
      
      colnames(param_vec) <-  c("Coverage", "Coverage_SMC", "Halflife", "Efficacy")
      
      ul = predict(x = as.matrix(param_vec), GP_trained_UL)$mean
      ci = predict(x = as.matrix(param_vec), GP_trained_CI)$mean
      
      return(ul - ci)
      
    }
    
    
    get_noninf_sd_plus= function(x) {
      params <- scenario
      param_vec = matrix(unlist(c(params[c(1,2)],x, params[c(3)])), nrow = 1)
      
      colnames(param_vec) <-  c("Coverage", "Coverage_SMC", "Halflife", "Efficacy")
      
      ul = predict(x = as.matrix(param_vec), GP_trained_UL)
      ulpred = ul$mean + sqrt(2*ul$sd2 + ul$nugs)
      ci = predict(x = as.matrix(param_vec), GP_trained_CI)
      cipred = ci$mean + sqrt(2*ci$sd2 + ci$nugs)
      
      return(ulpred - cipred)
      
    }
    
    get_noninf_sd_minus= function(x) {
      params <- scenario
      param_vec = matrix(unlist(c(params[c(1,2)],x, params[c(3)])), nrow = 1)
      
      colnames(param_vec) <-  c("Coverage", "Coverage_SMC", "Halflife", "Efficacy")
      
      ul = predict(x = as.matrix(param_vec), GP_trained_UL)
      ulpred = ul$mean - sqrt(2*ul$sd2 + ul$nugs)
      ci = predict(x = as.matrix(param_vec), GP_trained_CI)
      cipred = ci$mean - sqrt(2*ci$sd2 + ci$nugs)
      
      return(ulpred - cipred)
      
    }
    #this is to check if optimisation is possible
    max_noninf= function(x) {
      params <- scenario
      param_vec = matrix(unlist(c( params[c(1,2)], 150, params[c(3)])), nrow = 1)
      
      colnames(param_vec) <-  c("Coverage", "Coverage_SMC", "Halflife", "Efficacy")
      
      ul = predict(x = as.matrix(param_vec), GP_trained_UL)$mean
      ci = predict(x = as.matrix(param_vec), GP_trained_CI)$mean
      
      return( - (ul - ci))
    }
    
    
    return_coverage = function(x) {
      return(x)
    }
    
    #check if non inferiority is possible and efficacy over 40% 
    
    
    mincpppy= function(x) {
      params <- scenario
      param_vec = matrix(unlist(c( params[c(1,2)], x, params[c(3)])), nrow = 1)
      
      colnames(param_vec) <-  c("Coverage", "Coverage_SMC", "Halflife", "Efficacy")
      
      mincpppy = predict(x = as.matrix(param_vec), GP_trained_LAI)$mean

      
      return( mincpppy)
    }
    mincpppySMC= function(x) {
      params <- scenario
      param_vec = matrix(0.6)
      
      colnames(param_vec) <-  c( "Coverage_SMC")
      
      mincpppy = predict(x = as.matrix(param_vec), GP_trained_SMC)$mean
      
      
      return( mincpppy)
    }
    mincpppySMC <-  mincpppySMC(0.6)
    
basecpppy <- cpppy_low_access[which(EIRs_low_acces==EIR)]
    if (- max_noninf(150) > 0 ){
      
      ans_mean = gosolnp(pars  = NULL, fixed = NULL, fun = return_coverage, ineqfun = get_noninf, 
                         ineqLB = 0 , ineqUB = 1, 
                         LB = 30, UB = 150, distr = rep(1, 1), distr.opt = list(), 
                         n.restarts = 2, control = list(maxit = 100), n.sim = 200)
      
    incred <- 1-(  mincpppy(ans_mean$pars)/basecpppy)
    
    add_incred <- 1-(  mincpppy(ans_mean$pars)/mincpppySMC)
      # ans_sd_plus = gosolnp(pars  = NULL, fixed = NULL, fun = return_coverage, ineqfun = get_noninf_sd_plus, 
      #                    ineqLB = 0 , ineqUB = 1, 
      #                    LB = 30, UB = 150, distr = rep(1, 1), distr.opt = list(), 
      #                    n.restarts = 1, control = list(maxit = 100), n.sim = 500)
      # 
      # ans_sd_minus = gosolnp(pars  = NULL, fixed = NULL, fun = return_coverage, ineqfun = get_noninf_sd_minus, 
      #                    ineqLB = 0 , ineqUB = 1, 
      #                    LB = 30, UB = 150, distr = rep(1, 1), distr.opt = list(), 
      #                    n.restarts = 3, control = list(maxit = 100), n.sim = 200)
      # 
      
      scenarios[i,]$optimal_lai_HL <- ans_mean$pars
      scenarios[i,]$inc_red <- incred
      scenarios[i,]$add_inc_red <- add_incred
    #  scenarios[i,]$optimal_lai_HL_lowsd <- ans_sd_minus$pars
     # scenarios[i,]$optimal_lai_HL_highsd <- ans_sd_plus$pars
      
    } else {
      scenarios[i,]$optimal_lai_HL <- -1
      scenarios[i,]$inc_red <- -1
      scenarios[i,]$add_incred <- -1
      
    }
    
  }
  



write.table(scenarios, file = outfile_scenarios)

}

}
