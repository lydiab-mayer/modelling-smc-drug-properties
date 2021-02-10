library(hetGP)
library(ggplot2)
library(ggpubr)
library(reshape2)
library(tgp)


# plot non-inferiority results

rm(list = ls())

mainDir  <- "~/smc_lai/analysis_workflow/analysis_scripts"
setwd(mainDir)
source('~/smc_lai/analysis_workflow/3_GP_train/GP_toolbox.R')
source('~/smc_lai/analysis_workflow/analysis_scripts/supp/E5_resources.R')
dir.create(file.path(paste0(mainDir,"/Outputs/E5_Compare/")), showWarnings = FALSE)


#EIR = c(5,20,70,100,200)
SMC_HL <- c(10,20,31.295)
LAI_dec <- c("hill","exp")
seasonaility <- c("Mali","Sen")
settings <- expand.grid(SMC_HL,LAI_dec,seasonaility)
names(settings) <- c("SMC_HL","LAI_dec","seasonality")

k=3

for (k in 1:nrow(settings)) {
  setting <- settings[k, ]
  train_GP_E5_incred(setting)
}



for (k in 1:nrow(settings)) {
  setting <- settings[k, ]
  pred_GP_E5_incred(setting)
}


  
  
 