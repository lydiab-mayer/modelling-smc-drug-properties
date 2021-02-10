
rm(list = ls())
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


for (k in 1:nrow(scenarios)){
  
  seasonality <- scenarios[k,1]
  decay <- scenarios[k,2]
  age <- scenarios[k,3]
  acess <- scenarios[k,4]
  EIR <- scenarios[k,5]
  
  scen_name <- paste(seasonality,"_",age,"_",decay,"_",acess,"_",EIR, sep = "")

  filename_UL <- paste('/scicore/home/penny/GROUP/smc_lai/E3E4disc_comparison/gp_trained/gp_trained_UL_', scen_name, '.RData', sep = "")
  filename_CI <- paste('/scicore/home/penny/GROUP/smc_lai/E3E4disc_comparison/gp_trained/gp_trained_CI_', scen_name, '.RData', sep = "")
  filename_pppy_SMC <- paste('/scicore/home/penny/GROUP/smc_lai/E3E4disc_comparison/gp_trained/gp_trained_pppy_SMC_', scen_name, '.RData', sep = "")
  filename_pppy_LAI <- paste('/scicore/home/penny/GROUP/smc_lai/E3E4disc_comparison/gp_trained/gp_trained_pppy_LAI_', scen_name, '.RData', sep = "")
  
  filename_optimisation <- paste('/scicore/home/penny/GROUP/smc_lai/E3E4disc_comparison/outfiles_optimisation_summary/scenarios_all_',scen_name,".txt",sep = "")
  
  outfile <- paste("outfiles_final_analysis/figure_4_", scen_name,".pdf", sep = "")
  

  
  if (file.exists(filename_UL) & file.exists(filename_optimisation)){
    print(scen_name)
    
    load(filename_UL)
    load(filename_CI)
    load(filename_pppy_SMC)
    load(filename_pppy_LAI)
    
    optim_results <- read.table(filename_optimisation, header = T)
    
    #plot left (heatmap)
    SMCcoverage_heatmap <- 0.6
    
    LAI_1 <- c(90, 0.8)
    LAI_2 <- c(100, 0.9)
    LAI_3 <- c(70, 0.7)
    LAI_4 <- c(150, 1)
    
    n_lai_props <- 4
    
    my_colours <- c("tan4","turquoise4","yellowgreen","gold1")
    index_1 <- which(optim_results$Coverage_SMC == SMCcoverage_heatmap & optim_results$Halflife == LAI_1[1] & optim_results$Efficacy == LAI_1[2])
    index_2 <- which(optim_results$Coverage_SMC == SMCcoverage_heatmap & optim_results$Halflife == LAI_2[1] & optim_results$Efficacy == LAI_2[2])
    index_3 <- which(optim_results$Coverage_SMC == SMCcoverage_heatmap & optim_results$Halflife == LAI_3[1] & optim_results$Efficacy == LAI_3[2])
    index_4 <- which(optim_results$Coverage_SMC == SMCcoverage_heatmap & optim_results$Halflife == LAI_4[1] & optim_results$Efficacy == LAI_4[2])


    p11 <- ggplot(data = optim_results[optim_results$Coverage_SMC == SMCcoverage_heatmap,]) +
      geom_tile(aes(x = Halflife, y = Efficacy, fill = optimal_lai_coverage)) +
      geom_point(data = optim_results[index_1,], aes(x = Halflife, y = Efficacy), colour = my_colours[1] )+
      geom_point(data = optim_results[index_2,], aes(x = Halflife, y = Efficacy), colour = my_colours[2]  )+
      geom_point(data = optim_results[index_3,], aes(x = Halflife, y = Efficacy), colour = my_colours[3]  )+
      geom_point(data = optim_results[index_4,], aes(x = Halflife, y = Efficacy), colour = my_colours[4]  )+
      scale_fill_gradientn(limits = c(0.4,1),
                       colours=c("navyblue",  "darkorange1"),
                       na.value = "grey") +
      labs(x = "Halflife of protective efficacy [d]", y = "Initial protective efficacy [%]", fill= "optimal LAI\ncoverage") +
      theme(legend.position="bottom")

    #plot right top (optimal LAI for different SMC)
    SMCs <- unique(optim_results$Coverage_SMC)
    
    plot_LAI <- data.frame(rep(SMCs, n_lai_props))
    colnames(plot_LAI) <- c("SMC")
    plot_LAI$HL <- rep(c(LAI_1[1],LAI_2[1],LAI_3[1],LAI_4[1]), each = length(SMCs)) 
    plot_LAI$Ef <- rep(c(LAI_1[2],LAI_2[2],LAI_3[2],LAI_4[2]), each = length(SMCs))
    
    plot_LAI$opt_LAI_cov <- 0
    
    for (k in 1:nrow(plot_LAI)){
      index <- which(optim_results$Halflife == plot_LAI[k,]$HL &
                   optim_results$Efficacy == plot_LAI[k,]$Ef &
                   optim_results$Coverage_SMC == plot_LAI[k,]$SMC)
      plot_LAI[k,]$opt_LAI_cov <- optim_results[index, ]$optimal_lai_coverage
      }
    
    plot_LAI$LAI_props <- rep(seq(1,n_lai_props), each = nrow(plot_LAI)/n_lai_props)
    plot_LAI[plot_LAI == -1] <- NA
    
    my_labs <- c()
    my_lab <- paste("HL: ", LAI_1[1], ", EF: ",LAI_1[2], sep = "")
    my_labs <- c(my_labs,my_lab)
    
    my_lab <- paste("HL: ", LAI_2[1], ", EF: ",LAI_2[2], sep = "")
    my_labs <- c(my_labs,my_lab)
    
    my_lab <- paste("HL: ", LAI_3[1], ", EF: ",LAI_3[2], sep = "")
    my_labs <- c(my_labs,my_lab)
    
    my_lab <- paste("HL: ", LAI_4[1], ", EF: ",LAI_4[2], sep = "")
    my_labs <- c(my_labs,my_lab)
    
    p12 <- ggplot(data = plot_LAI) +
      geom_line(aes(x = SMC, y = opt_LAI_cov, colour = as.factor(LAI_props)), size = 1) +
      geom_vline(aes(xintercept = 0.6), colour = "grey50", linetype = "dashed") +
      scale_colour_manual(values = my_colours, labels = my_labs)+
      labs( y = "opt. LAI coverage", x = "SMC coverage", colour = "LAI properties")+
      theme(legend.position="none")+
      theme(axis.title = element_text(size=9)) 
    

    
    #plot right middle (relative difference cases per person per year with 20% above optimal coverage)
    plot_LAI_cases <- plot_LAI
    plot_LAI_cases$cases_SMC <- 0
    plot_LAI_cases$cases_LAI <- 0
    
    params <- data.frame(cbind(plot_LAI_cases$SMC, plot_LAI_cases$Ef,
                           plot_LAI_cases$HL, plot_LAI_cases$opt_LAI_cov))
    
    colnames(params) <- c("Coverage_SMC", "Efficacy","Halflife","Coverage")
    
    for (k in 1:nrow(params)){
      params[k,]$Coverage <- min(c(1,params[k,]$Coverage + 0.2*params[k,]$Coverage))
      }
    
    index <- which(is.na(params$Coverage) ==F)
    
    cases_SMC <- predict(x = as.matrix(params[index,c(1)]), 
                     object = GP_trained_SMC)
    
    cases_LAI <- predict(x = as.matrix(params[index,c(1,2,3,4)]), 
                     object = GP_trained_LAI)
    plot_LAI_cases$cases_SMC <- NA
    plot_LAI_cases$cases_LAI <- NA
    
    plot_LAI_cases[index,]$cases_SMC <- cases_SMC$mean
    plot_LAI_cases[index,]$cases_LAI <- cases_LAI$mean

    plot_LAI_cases$cases_relative_difference <- (plot_LAI_cases$cases_SMC - plot_LAI_cases$cases_LAI)/plot_LAI_cases$cases_SMC 
    
    p13 <- ggplot(data = plot_LAI_cases) +
      geom_line(aes(x = SMC, y = cases_relative_difference, colour = as.factor(LAI_props)), size = 1) +
      geom_hline(aes(yintercept = 0),  colour = "grey50")+
      geom_vline(aes(xintercept = 0.6), colour = "grey50", linetype = "dashed") +
      scale_colour_manual(values = my_colours, labels = my_labs)+
      labs(y = "rel. diff. in cpppy",
       x = "SMC coverage", colour = "LAI properties")+
      theme(legend.position="none")+
      theme(axis.title = element_text(size=9)) 
    
    #plot right bottom
    LAI_coverages <- seq(0.4,1,0.05)
    plot_LAI_cases_absolute <- data.frame(rep(LAI_coverages, n_lai_props))
    colnames(plot_LAI_cases_absolute) <- c("LAI_coverage")
    
    plot_LAI_cases_absolute$HL <- rep(c(LAI_1[1],LAI_2[1],LAI_3[1],LAI_4[1]), each = length(LAI_coverages)) 
    plot_LAI_cases_absolute$Ef <- rep(c(LAI_1[2],LAI_2[2],LAI_3[2],LAI_4[2]), each = length(LAI_coverages))
    plot_LAI_cases_absolute$SMC <- 0.6
    
    plot_LAI_cases_absolute$LAI_props <- rep(seq(1,n_lai_props), each = nrow(plot_LAI_cases_absolute)/n_lai_props)  
    
    #get cases per person per year for different coverage levels
    params <- data.frame(cbind(plot_LAI_cases_absolute$SMC, 
                           plot_LAI_cases_absolute$Ef, plot_LAI_cases_absolute$HL, 
                           plot_LAI_cases_absolute$LAI_coverage))

colnames(params) <- c("Coverage_SMC", "Efficacy","Halflife","Coverage")
cases_LAI <- predict(x = as.matrix(params[,c(1,2,3,4)]), object = GP_trained_LAI)

plot_LAI_cases_absolute$cases_LAI <- cases_LAI$mean

#get cases per person per year at optimal coverage level for non inferiority
plot_LAI_cases_absolute$opt_LAI_cov <- 0

for (k in 1:nrow(plot_LAI_cases_absolute)){
  index <- which(optim_results$Halflife == plot_LAI_cases_absolute[k,]$HL &
                   optim_results$Efficacy == plot_LAI_cases_absolute[k,]$Ef &
                   optim_results$Coverage_SMC == plot_LAI_cases_absolute[k,]$SMC)
  plot_LAI_cases_absolute[k,]$opt_LAI_cov <- optim_results[index, ]$optimal_lai_coverage
}

#get cases per person per year in SMC case (60%)
plot_LAI_cases_absolute$opt_LAI_cov <- 0

for (k in 1:nrow(plot_LAI_cases_absolute)){
  index <- which(optim_results$Halflife == plot_LAI_cases_absolute[k,]$HL &
                   optim_results$Efficacy == plot_LAI_cases_absolute[k,]$Ef &
                   optim_results$Coverage_SMC == plot_LAI_cases_absolute[k,]$SMC)
  plot_LAI_cases_absolute[k,]$opt_LAI_cov <- optim_results[index, ]$optimal_lai_coverage
}

params <- data.frame(cbind(plot_LAI_cases_absolute$SMC, 
                           plot_LAI_cases_absolute$Ef, plot_LAI_cases_absolute$HL, 
                           plot_LAI_cases_absolute$opt_LAI_cov))

colnames(params) <- c("Coverage_SMC", "Efficacy","Halflife","Coverage")
cases_LAI <- predict(x = as.matrix(params[,c(1,2,3,4)]), object = GP_trained_LAI)

plot_LAI_cases_absolute$cases_LAI_optimal <- cases_LAI$mean





p14 <- ggplot(data = plot_LAI_cases_absolute) +
  geom_line(aes(x = LAI_coverage, y = cases_LAI, colour = as.factor(LAI_props)), size = 1) +
  geom_vline(aes(xintercept = opt_LAI_cov, colour = as.factor(LAI_props)), linetype = "dashed",show_guide = FALSE) +
  scale_colour_manual(values = my_colours, labels = my_labs)+
  labs(y = "cpppy",
       x = "LAI coverage", colour = "LAI\nproperties") +
  xlim(c(0.4,1))+
  theme(legend.position="bottom")+
  theme(axis.title = element_text(size=9)) 

p11 <- ggarrange(p11, labels = c("A"))
p <- ggarrange(p12,p13,p14, ncol = 1, nrow = 3,  heights = c(0.9,0.9,1), labels = c("B","C","D"))
p <- ggarrange(p11,p, ncol = 2, nrow = 1, widths = c(0.9,1))

ggsave(filename = outfile, p, width = 12, height = 7)

}
}


