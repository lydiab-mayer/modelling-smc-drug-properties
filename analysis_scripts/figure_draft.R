
rm(list = ls())
setwd("~/smc_lai/analysis_workflow/analysis_scripts")
library(plyr)
library(ggplot2)
library(coda)
library(ggpubr)
library(pracma)
library(gridExtra)
library(reshape2)
library(hetGP)
library(Rsolnp)

seasonalities <- c("Mali","Sen")
decays <- c("exp","hill","wei")
ages <- c(4.9167)
acesses <- c(0.1,0.5)
EIRs <- c(3,4,5,8,9,20,28,47,150)

EIRs_low_acces <- c(3,4,8,28)
EIRs_high_acces <- c(5,9,20,47)


scenarios <- expand.grid(seasonalities,decays,ages,acesses,EIRs)
colnames(scenarios) <- c("Seasonality","Decay","Age","Acess","EIR")

index <- which(scenarios$Acess == 0.1 & scenarios$EIR %in% EIRs_high_acces)
scenarios <- scenarios[-index,]

index <- which(scenarios$Acess == 0.5 & scenarios$EIR %in% EIRs_low_acces)
scenarios <- scenarios[-index,]

outfile <- "outfiles_final_analysis/figure_draft.pdf"
optimisation_table <- c()
for (k in 1:nrow(scenarios)){
  
  seasonality <- scenarios[k,1]
  decay <- scenarios[k,2]
  age <- scenarios[k,3]
  acess <- scenarios[k,4]
  EIR <- scenarios[k,5]
  
  scen_name <- paste(seasonality,"_",age,"_",decay,"_",acess,"_",EIR, sep = "")
  
 
  filename_optimisation <- paste('/scicore/home/penny/GROUP/smc_lai/E3E4disc_comparison/outfiles_optimisation_summary/scenarios_all_',scen_name,".txt",sep = "")
  
  if (file.exists(filename_optimisation)){
    print(scen_name)
  a <- read.table(filename_optimisation, header = T)
  a$seasonality <- rep(seasonality,nrow(a))
  a$decay <- rep(decay,nrow(a))
  a$acess <- rep(acess,nrow(a))
  if (acess == 0.1){
    EIR_low <- which(c(EIRs_low_acces,150) == EIR)
    a$EIR <- rep(EIR_low,nrow(a))
  }

  if (acess == 0.5){
    EIR_high <- which(c(EIRs_high_acces,150) == EIR)
    a$EIR <- rep(EIR_high,nrow(a))
  }

  
  
  optimisation_table <- rbind(optimisation_table,a)
  }
  

 
  
  
 

  

} 

index <- which(optimisation_table$optimal_lai_coverage == -1 )
optimisation_table[index,]$optimal_lai_coverage <- rep(100,length(index))
index <- which(optimisation_table$optimal_lai_coverage == 0 )
optimisation_table[index,]$optimal_lai_coverage <- rep(100,length(index))


index_sorted <- order(optimisation_table$optimal_lai_coverage, decreasing = FALSE)
optimisation_table_sorted <- optimisation_table[index_sorted,]

top_1000 <- optimisation_table_sorted[1:1000,c(1,2,3,5,6,7,8)]
  
  
smc_coverages_df <- data.frame(c(0.4,0.5,0.6,0.7,0.8,0.9,1))
colnames(smc_coverages_df) <- c("SMC_Coverage")
smc_coverages_df$frequency <- 0
for (k in 1:nrow(smc_coverages_df)){
  smc_coverages_df[k,]$frequency <- length(which(top_1000$Coverage_SMC == smc_coverages_df[k,]$SMC_Coverage))/
    length(which(optimisation_table$Coverage_SMC == smc_coverages_df[k,]$SMC_Coverage))
}
smc_coverages_df$variable <- "SMC_Coverage"

halflife_df <- data.frame(seq(30, 150, by = 10))
colnames(halflife_df) <- c("Halflife")
halflife_df$frequency <- 0
for (k in 1:nrow(halflife_df)){
  halflife_df[k,]$frequency <- length(which(top_1000$Halflife == halflife_df[k,]$Halflife))/
    length(which(optimisation_table$Halflife == halflife_df[k,]$Halflife))
}

halflife_df$variable <- "Halflife"


efficacy_df <- data.frame(seq(0.7,1, by = 0.1))
colnames(efficacy_df) <- c("Efficacy")
efficacy_df$frequency <- 0

efficacy_df[1,]$frequency <- length(which(top_1000$Efficacy == 0.7))/length(which(optimisation_table$Efficacy == 0.7))
efficacy_df[2,]$frequency <- length(which(top_1000$Efficacy == 0.8))/length(which(optimisation_table$Efficacy == 0.8))
efficacy_df[3,]$frequency <- length(which(top_1000$Efficacy == 0.9))/length(which(optimisation_table$Efficacy == 0.9))
efficacy_df[4,]$frequency <- length(which(top_1000$Efficacy == 1))/length(which(optimisation_table$Efficacy == 1))

efficacy_df$variable <- "Efficacy"


seasonality_df <- data.frame(c("Mali","Sen"))
colnames(seasonality_df) <- c("seasonality")
seasonality_df$frequency <- 0
for (k in 1:nrow(seasonality_df)){
  seasonality_df[k,]$frequency <- length(which(top_1000$seasonality == seasonality_df[k,]$seasonality))/
    length(which(optimisation_table$seasonality == seasonality_df[k,]$seasonality))
}

seasonality_df$variable <- "Seasonality"


decay_df <- data.frame(c("exp","hill","wei"))
colnames(decay_df) <- c("decay")
decay_df$frequency <- 0
for (k in 1:nrow(decay_df)){
  decay_df[k,]$frequency <- length(which(top_1000$decay == decay_df[k,]$decay))/
    length(which(optimisation_table$decay == decay_df[k,]$decay))
}

decay_df$variable <- "Decay"


acess_df <- data.frame(c(0.1,0.5))
colnames(acess_df) <- c("acess")
acess_df$frequency <- 0
for (k in 1:nrow(acess_df)){
  acess_df[k,]$frequency <- length(which(top_1000$acess == acess_df[k,]$acess))/
    length(which(optimisation_table$acess == acess_df[k,]$acess))
}

acess_df$variable <- "Access"

colnames(smc_coverages_df) <- c("value","frequency","variable")
colnames(halflife_df) <- c("value","frequency","variable")
colnames(efficacy_df) <- c("value","frequency","variable")
colnames(seasonality_df) <- c("value","frequency","variable")
colnames(decay_df) <- c("value","frequency","variable")
colnames(acess_df) <- c("value","frequency","variable")



plot_df <- rbind(smc_coverages_df, halflife_df, efficacy_df, seasonality_df, decay_df, acess_df)
plot_df[plot_df$variable == "Halflife",]$value <- c("030","040","050","060","070","080","090","100","110","120","130","140","150")

plot_df$reference <- 0
index <- which(plot_df$variable == "SMC_Coverage")
plot_df[index,]$reference <- sum(plot_df[index,]$frequency)/length(index)

index <- which(plot_df$variable == "Halflife")
plot_df[index,]$reference <- sum(plot_df[index,]$frequency)/length(index)

index <- which(plot_df$variable == "Efficacy")
plot_df[index,]$reference <- sum(plot_df[index,]$frequency)/length(index)

index <- which(plot_df$variable == "Seasonality")
plot_df[index,]$reference <- sum(plot_df[index,]$frequency)/length(index)

index <- which(plot_df$variable == "Decay")
plot_df[index,]$reference <- sum(plot_df[index,]$frequency)/length(index)

index <- which(plot_df$variable == "Access")
plot_df[index,]$reference <- sum(plot_df[index,]$frequency)/length(index)


p <- ggplot(data = plot_df) +
  geom_bar(aes(x = value, y = frequency),position="dodge",stat="identity", fill = "lightblue4")+
  geom_hline(aes(yintercept = reference), colour = "red", linetype = "dashed")+
  facet_wrap("variable",scales="free") +
  xlab("value")

  

    
  
ggsave(filename = outfile, p, width = 10, height = 7)


