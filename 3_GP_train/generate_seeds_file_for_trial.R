rm(list = ls())


source('~/smc_lai/analysis_workflow/3_GP_train/GP_toolbox.R')
source('~/smc_lai/analysis_workflow/analysis_scripts/supp/calcdifferencetrial.R')

load("/scicore/home/smith/GROUP/smc_lai/E5_2_CT_LAI/param_ranges.RData")
# generate a table to run the GP on the difference between two columns in seeds tables from different experiments
EIR = c(5,9,20,47,150)
SMC_HL <- c(10,20,31.295)
LAI_dec <- c("hill","exp","wei")
seasonaility <- c("Mali","Sen")
scenarios <- expand.grid(EIR,SMC_HL,LAI_dec,seasonaility)
names(scenarios) <- c("EIR","SMC_HL","LAI_dec","seasonality")

calcdifferencetrial(scenarios)


