#### train gps
library(plyr)
library(ggplot2)
library(coda)
library(ggpubr)
library(pracma)
library(gridExtra)
library(reshape2)
library(hetGP)
library(Rsolnp)

args = commandArgs(TRUE)
seasonality_in = args[1]
decay_in = args[2]
acess_in = args[3]
age_in = args[4]
EIR_in = args[5]


if (seasonality_in == 1){
  seasonality <- "Mali"
}

if (seasonality_in == 2){
  seasonality <- "Sen"
}

if (decay_in == 1){
  decay <- "exp"
}

if (decay_in == 2){
  decay <- "hill"
}

if (decay_in == 3){
  decay <- "wei"
}


if (acess_in == 1){
  acess <- 0.1
}

if (acess_in == 2){
  acess <- 0.5
}


if (age_in == 1){
  age <- 4.9167
}

if (EIR_in == 1 & acess == 0.1){
  EIR <- 3
} 

if (EIR_in == 1 & acess == 0.5) {
  EIR <- 5
}

if (EIR_in == 2 & acess == 0.1){
  EIR <- 4
} 
if (EIR_in == 2 & acess == 0.5){
  EIR <- 9
}

if (EIR_in == 3 & acess == 0.1){
  EIR <- 8
} 

if (EIR_in == 3 & acess == 0.5){
  EIR <- 20
}


if (EIR_in == 4 & acess == 0.1){
  EIR <- 28
} 

if (EIR_in ==4 & acess == 0.5){
  EIR <- 47
}

if (EIR_in == 5){
  EIR <- 150
} 






scen_name <- paste(seasonality,"_",age,"_",decay,"_",acess,"_",EIR, sep = "")
print(scen_name)

filename_SMC <- paste('/scicore/home/penny/GROUP/smc_lai/E3_SMCSMCdisc/postprocessing_5/seeds_E3SMCSMCdisc_', scen_name, '.txt', sep = "")
filename_LAI <- paste('/scicore/home/penny/GROUP/smc_lai/E4_SMCLAIdisc/postprocessing_5/seeds_E4SMCLAIdisc_', scen_name, '.txt', sep = "")
filename_comparison <- paste("/scicore/home/penny/GROUP/smc_lai/E3E4disc_comparison/postprocessing_5/seeds_E3E4disc_comparison_", scen_name,".txt", sep = "")

print(filename_SMC)
print(filename_LAI)
print(filename_comparison)
  
if (file.exists(filename_SMC) & file.exists(filename_LAI) & file.exists(filename_comparison)){

  print("files exist")
  seeds_SMC <- read.table(filename_SMC, header = T)
  seeds_LAI <- read.table(filename_LAI, header = T)
  seeds_comparison <- read.table(filename_comparison, header = T)
  
  outfile_UL <- paste('/scicore/home/penny/GROUP/smc_lai/E3E4disc_comparison/gp_trained/gp_trained_UL_', scen_name, '.RData', sep = "")
  outfile_CI <- paste('/scicore/home/penny/GROUP/smc_lai/E3E4disc_comparison/gp_trained/gp_trained_CI_', scen_name, '.RData', sep = "")
  outfile_pppy_SMC <- paste('/scicore/home/penny/GROUP/smc_lai/E3E4disc_comparison/gp_trained/gp_trained_pppy_SMC_', scen_name, '.RData', sep = "")
  outfile_pppy_LAI <- paste('/scicore/home/penny/GROUP/smc_lai/E3E4disc_comparison/gp_trained/gp_trained_pppy_LAI_', scen_name, '.RData', sep = "")
  
  #train GP on upper limit
  n_seeds <- length(unique(seeds_SMC$seed))
  n_points <- round(nrow(seeds_comparison)/n_seeds*0.75)
  index_train = sample(nrow(seeds_comparison)/n_seeds, n_points)
  index_test = setdiff(1:nrow(seeds_comparison)/n_seeds, index_train)
  
  train_data_comparison <-seeds_comparison[seeds_comparison$scen_id %in% index_train, ]
  test_data_comparison <- seeds_comparison[seeds_comparison$scen_id %in% index_test, ]
  
  param_col <- c("Coverage", "Coverage_SMC", "Halflife", "Efficacy")
  response_col <- "UL"
  
  prdata = find_reps(X = as.matrix(train_data_comparison[, param_col]),
                   Z = as.matrix(train_data_comparison[, response_col]),
                   rescale = FALSE, normalize = FALSE)
  
  GP_trained_UL = mleHetGP(X = list(X0 = as.matrix(prdata$X0),
                                  Z0 = as.matrix(prdata$Z0), 
                                  mult = prdata$mult),
                         Z = prdata$Z, 
                         lower = rep(0.0001, length(param_col)), 
                         upper = rep(10, length(param_col) ),
                         covtype = "Matern5_2")
  
  save(GP_trained_UL, file = outfile_UL)
  
  #train GP on CI_high_HR
  param_col <- c("Coverage", "Coverage_SMC", "Halflife", "Efficacy")
  response_col <- "CI_high_HR"
  
  prdata = find_reps(X = as.matrix(train_data_comparison[, param_col]),
                   Z = as.matrix(train_data_comparison[, response_col]),
                   rescale = FALSE, normalize = FALSE)
  
  GP_trained_CI = mleHetGP(X = list(X0 = as.matrix(prdata$X0),
                                  Z0 = as.matrix(prdata$Z0), 
                                  mult = prdata$mult),
                         Z = prdata$Z, 
                         lower = rep(0.0001, length(param_col)), 
                         upper = rep(10, length(param_col) ),
                         covtype = "Matern5_2")
  
  save(GP_trained_CI, file = outfile_CI)



  #train GP on pppy SMC
  
  n_seeds_SMC <- length(unique(seeds_SMC$seed))
  n_points_SMC = round(nrow(seeds_SMC)/n_seeds_SMC*0.75)
  seeds_SMC$id <- rep(seq(1,nrow(seeds_SMC)/n_seeds_SMC),each=n_seeds_SMC)
  
  index_train_SMC = sample(nrow(seeds_SMC)/n_seeds_SMC, n_points_SMC)
  index_test_SMC = setdiff(1:nrow(seeds_SMC)/n_seeds_SMC, index_train_SMC)
  
  train_data_SMC <-seeds_SMC[seeds_SMC$id %in% index_train_SMC,]
  test_data_SMC <- seeds_SMC[seeds_SMC$id %in% index_test_SMC,]
  
  param_col <- c("Coverage_SMC")
  response_col <- c("pppy_y10_all")
  
  prdata = find_reps(X = as.matrix(train_data_SMC[, param_col]),
                   Z = as.matrix(train_data_SMC[, response_col]),
                   rescale = FALSE, normalize = FALSE)
  
  GP_trained_SMC = mleHetGP(X = list(X0 = as.matrix(prdata$X0),
                                   Z0 = as.matrix(prdata$Z0),
                                   mult = prdata$mult),
                          Z = prdata$Z,
                          lower = rep(0.0001, length(param_col)),
                          upper = rep(10, length(param_col) ),
                          covtype = "Matern5_2")

  save(GP_trained_SMC, file = outfile_pppy_SMC)


  #train GP pppy LAI
  n_seeds_LAI <- length(unique(seeds_LAI$seed))
  n_points_LAI = round(nrow(seeds_LAI)/n_seeds_LAI*0.75)
  seeds_LAI$id <- rep(seq(1,nrow(seeds_LAI)/n_seeds_LAI),each=n_seeds_LAI)
  
  index_train_LAI = sample(nrow(seeds_LAI)/n_seeds_LAI, n_points_LAI)
  index_test_LAI = setdiff(1:nrow(seeds_LAI)/n_seeds_LAI, index_train_LAI)
  
  train_data_LAI <-seeds_LAI[seeds_LAI$id %in% index_train_LAI,]
  test_data_LAI <- seeds_LAI[seeds_LAI$id %in% index_test_LAI,]
  
  param_col <- c("Coverage_SMC", "Efficacy","Halflife","Coverage")
  response_col <- "pppy_y10_all"
  
  prdata = find_reps(X = as.matrix(train_data_LAI[, param_col]),
                   Z = as.matrix(train_data_LAI[, response_col]),
                   rescale = FALSE, normalize = FALSE)
  
  GP_trained_LAI = mleHetGP(X = list(X0 = as.matrix(prdata$X0),
                                   Z0 = as.matrix(prdata$Z0),
                                   mult = prdata$mult),
                          Z = prdata$Z,
                          lower = rep(0.0001, length(param_col)),
                          upper = rep(10, length(param_col) ),
                          covtype = "Matern5_2")

save(GP_trained_LAI, file = outfile_pppy_LAI)


}
