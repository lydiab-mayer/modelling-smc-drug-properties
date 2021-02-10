

rm(list = ls())

setwd("~/workspace/smc_lai/analysis_june")

library(nloptr)
library(hetGP)
library(Rsolnp)
library(metaheuristicOpt)
library(tools)
library(rapportools)
library(stringr)
library(ggplot2)
library(ggpubr)


seasonalities <- c("Mali", "Sen")
decays <- c("exp", "hill") 

#acesses <- c("0.1","0.5")
acesses <- c("0.1")
#ages <- c("4.9167", "9.9167")
ages <- c("4.9167")

settings <- expand.grid(seasonalities, decays, acesses, ages)

plots <- list()
for (l in 1:nrow(settings)){
  print(l)
  
  seasonality <- settings[l,1]
  decay <- settings[l,2]
  acess <- settings[l,3]
  age <- settings[l,4]
  
  scen_name <- paste(seasonality,"_",age,"_",decay,"_",acess, sep = "")
  
  scenarios <- read.table(paste('scenarios_processed/scenarios_', scen_name, '.txt', sep = ""), header = T)
  
  
  scenarios_plot <- scenarios
  index <- which(scenarios$EIR == 11 | scenarios$EIR == 51)
  scenarios_plot <- scenarios_plot[index,]
  index <- which(scenarios_plot$Coverage_SMC_Res == 0.4 | scenarios_plot$Coverage_SMC_Res == 0.8)
  scenarios_plot <- scenarios_plot[index,]
  

  EIR.labs <- c("EIR = 11", "EIR = 51")
  names(EIR.labs) <- c(11,51)
  
  SMC.labs <- c("SMC coverage = 0.4", "SMC coverage = 0.8")
  names(SMC.labs) <- c(0.4,0.8)
  
  
  plots[[l]] <- ggplot(data = scenarios_plot) +
    geom_tile(aes(x = Halflife, y = Efficacy, fill = optimal_lai_coverage)) +
    facet_grid( EIR ~ Coverage_SMC_Res, labeller = labeller(EIR = EIR.labs, Coverage_SMC_Res = SMC.labs)) + 
    scale_fill_gradientn(limits = c(0.4,1),
                         colours=c("green4",  "grey"),
                         na.value = "white") +
    labs(title = paste("Seasonality: ", seasonality, ", Decay: ", decay, sep =""), fill= "LAI coverage")
  
  
 
  
}

p1 <- plots[[1]]
p2 <- plots[[2]]
p3 <- plots[[3]]
p4 <- plots[[4]]

p <- ggarrange(p1,p2,p3,p4, nrow = 2, ncol = 2)



ggsave(paste('outfiles/summary_plot_coverage.pdf', sep = ""), p, width = 15, height = 10)

