
rm(list = ls())

setwd("~/smc_lai/analysis_workflow/analysis_scripts")

library(nloptr)
library(hetGP)
library(Rsolnp)
library(metaheuristicOpt)
library(tools)
library(rapportools)
library(stringr)
library(ggplot2)


seasonalities <- c("Mali", "Sen")
decays <- c("exp", "hill") 


#acesses <- c("0.1","0.5")
acesses <- c("0.1")
#ages <- c("4.9167", "9.9167")
ages <- c("4.9167")
rows <- seq(26,9100,25) 

settings <- expand.grid(seasonalities, decays, acesses, ages)


for (k in 1:nrow(settings)){
  seasonality <- settings[k,1]
  decay <- settings[k,2]
  acess <- settings[k,3]
  age <- settings[k,4]
  
  scen_name <- paste(seasonality,"_",age,"_",decay,"_",acess, sep = "")
  
  row <- rows[1]
  scenarios <- read.table(paste('/scicore/home/smith/GROUP/smc_lai/plot_optimisation/outfiles/scenarios_', scen_name, '_', row,'.txt', sep = ""), header = T)

  
for (j in 1:length(rows)){
  row <- rows[j]
  a <- read.table(paste('/scicore/home/smith/GROUP/smc_lai/plot_optimisation/outfiles/scenarios_', scen_name, '_', row,'.txt', sep = ""), header = T)
  scenarios[(row - 25):row,] <- a[(row - 25):row,]

}

  row <- 9100
  a <- read.table(paste('/scicore/home/smith/GROUP/smc_lai/plot_optimisation/outfiles/scenarios_', scen_name, '_', row,'.txt', sep = ""), header = T)
  scenarios[(row - 24):row,] <- a[(row - 24):row,]


write.table(scenarios, paste('/scicore/home/smith/burlyd00/smc_lai/E3E4_comp/gp_5/optimisation/processed_ML/scenarios_',scen_name,".txt",sep = ""))

}



