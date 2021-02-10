# plot non-inferiority results
library(ggpubr)
library(dplyr)
library(hetGP)
library(ggplot2)
library(ggpubr)
library(reshape2)
library(tgp)

rm(list = ls())

mainDir  <- "~/smc_lai/analysis_workflow/analysis_scripts"
setwd(mainDir)
source('~/smc_lai/analysis_workflow/analysis_scripts/supp/import_functions.R')





Access <- c(0.1)
LAI_dec <- c("exp")
seasonaility <- c("Sen","Mali")
SMC_HL <- c(20,31.295)
settings <- expand.grid(Access,LAI_dec,seasonaility,SMC_HL)

names(settings) <- c("Access","LAI_dec","seasonality","SMC_HL")

points_high <-c(5,9,20,47,150)


listset <- list()
for (i in 1:nrow(settings)) {
  setting <- settings[i,]
  
  listEIR <- list()
  for(k in 1:length(points_high)) {
    
    EIR <- points_high[k]
  agg <- read.table(paste("/scicore/home/penny/GROUP/smc_lai/E5_1_CT_SMC/postprocessing/seeds_E5_1_CT_SMC_",setting[, "seasonality"],"_4.9167_",
                          setting[,"SMC_HL"],"_",
                            setting[, "LAI_dec"],"_",EIR,'.txt', sep = ""), 
                      header = T, as.is = TRUE, stringsAsFactors = FALSE)
  
  df<- data.frame(cbind(incred=mean(agg[,c("incred")]),setting,EIR))
  listEIR[[k]]<- df
  }
  
  dfset <- as.data.frame(do.call(rbind, listEIR))
  
  listset[[i]] <- dfset
}
  
dfsettings <- as.data.frame(do.call(rbind, listset))
