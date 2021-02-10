
rm(list = ls())
library(plyr)
library(ggplot2)
library(coda)
library(ggpubr)
library(pracma)
library(gridExtra)
library(reshape2)
library(hetGP)
library(data.table)
library(Rsolnp)

seasonality <- c("Mali","Sen") #Sen
LAI_dec <- c("exp","wei","hill")
IntAge <- c(4.9167)
Access <- c(0.1,0.5)


scenarios <- expand.grid(seasonality,LAI_dec,IntAge,Access)
colnames(scenarios) <- c("seasonality","LAI_dec","IntAge","Access")


points_low <-c(3,4,8,28,150)
points_high <- c(5,9,20,47,150)

setting <- scenarios[i,]

i=5
i=1


plot_tradeoffs_E4 <- function(setting){
  if(setting$Access==0.1) {EIRs <- points_low}else{EIRs <- points_high}
  dfsheat <- list()
  for (k in 1:length(EIRs)){
    
    seasonality <- setting[,1]
    decay <- setting[,2]
    age <- setting[,3]
    acess <- setting[,4]
    EIR <- EIRs[k]
    scen_name <- paste(seasonality,"_",age,"_",decay,"_",acess,"_",EIR, sep = "")
    
    filename_UL <- paste('/scicore/home/penny/GROUP/smc_lai/E3E4disc_comparison/gp_trained/gp_trained_UL_', scen_name, '.RData', sep = "")
    filename_CI <- paste('/scicore/home/penny/GROUP/smc_lai/E3E4disc_comparison/gp_trained/gp_trained_CI_', scen_name, '.RData', sep = "")
    filename_pppy_SMC <- paste('/scicore/home/penny/GROUP/smc_lai/E3E4disc_comparison/gp_trained/gp_trained_pppy_SMC_', scen_name, '.RData', sep = "")
    filename_pppy_LAI <- paste('/scicore/home/penny/GROUP/smc_lai/E3E4disc_comparison/gp_trained/gp_trained_pppy_LAI_', scen_name, '.RData', sep = "")
    
    filename_optimisation <- paste('/scicore/home/penny/GROUP/smc_lai/E3E4disc_comparison/outfiles_optimisation_summary/scenarios_all_',scen_name,".txt",sep = "")
    
    outfile <- paste("./Outputs/E4_tradeoffs/figure_", scen_name,"_tradeoffs.pdf", sep = "")
    
    
    
    if (file.exists(filename_UL) & file.exists(filename_optimisation)){
      print(scen_name)
      
      load(filename_UL)
      load(filename_CI)
      load(filename_pppy_SMC)
      load(filename_pppy_LAI)
      
      optim_results <- read.table(filename_optimisation, header = T)
      

    
  
 

      regimen = ifelse(seasonality=="Mali","long season","short season")
      decay = ifelse(LAI_dec=="exp","exponential LAIs",
                     ifelse( LAI_dec=="hill","sigmoid LAIs","bi-phasic LAIs"))
           #plot right middle (relative difference cases per person per year with 20% above optimal coverage)
      df_plot <- data.frame(seq(0.4,1,0.1))
      colnames(df_plot) <- c("LAI_coverage")
      
      Halflife <- c(30,50,70,90,110,130,150)
      Efficacy <- seq(0.7,1, 0.05)
      Coverage_LAI <- seq(0.4,1,0.1)
      LAI_spec <- expand.grid(Efficacy,Halflife,Coverage_LAI)
      LAI_specs <- data.frame(cbind(rep(0.6, nrow(LAI_spec)), LAI_spec))
      colnames(LAI_specs) <- c("Coverage_SMC", "Efficacy","Halflife","Coverage")
      
      cases_LAIs <- predict(x = as.matrix(LAI_specs[,c(1,2,3,4)]), object = GP_trained_LAI)
      
      LAI_specs$cases_LAIs <- cases_LAIs$mean
      
      #get SMC cases at 60% coverage
      
      params <- data.frame(0.6)
      colnames(params) <- c("Coverage_SMC")
      cases_SMC <- predict(x = as.matrix(params), 
                           object = GP_trained_SMC)

  

      LAI_specs$rel_diff <-  (cases_SMC$mean - LAI_specs$cases_LAIs)/cases_SMC$mean
      
      LAI_specs$inc_red <-  (1.4-LAI_specs$cases_LAIs)/1.4
      
      #plots for relative difference
      
     
      lin_reg <- function(x){
        reg <- lm(rel_diff ~ Coverage, data = x)
        slope <- reg$coefficients[2]
        regs <- data.frame(cbind(x, slope))
        return(regs)
      }
      
      slopes <-    LAI_specs %>%
        group_by(Efficacy, Halflife) %>%
        group_modify(~ lin_reg(.x))
      
     colorsEfficacy<- brewer.pal( n=9, name = "YlOrBr")
     colorsHalflife <- brewer.pal( n=9, name = "RdPu")
      
      ggplot(data = slopes) +
        geom_line(aes(x = Halflife, y = slope, color = as.factor(Efficacy*100)))+
      scale_colour_manual(values =colorsEfficacy[2:9])+
        
        labs(x = "Half-life of \nprotective efficacy [days]", y = "slope",
             color= "Initial protective \nefficacy [%]") +
        guides(colour = guide_legend(reverse = TRUE),
               size = guide_legend(title.position="top"))+
              theme_bw()+
        theme(panel.border = element_blank(),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              axis.line = element_line(colour = "black"),
              axis.text.x =element_text(size=15,colour="black",vjust=0.7),
              axis.text.y =element_text(size=15,colour="black"),
              axis.title=element_text(size=15,colour="black", face="bold"))+
        theme(strip.background = element_rect(colour="white", fill="white",
                                              size=1.5, linetype="solid"))+
        theme(strip.text = element_text(size=15, color="black", face="bold.italic"))+ 
        theme(legend.text = element_text( size=15),
              legend.title= element_text(size=15, face="bold"))+   
        theme(plot.title = element_text(hjust = 0.5,size=20,face="bold"))
      
      
      ggplot(data = slopes) +
        geom_line(aes(x = Efficacy*100, y = slope, color = as.factor(Halflife)))+
        scale_colour_manual(values =colorsHalflife[2:9])+
        
        labs(x = "Initial protective \nefficacy [%]", y = "slope",
             color= "Half-life of protective \nefficacy [days]") +
        guides(colour = guide_legend(reverse = TRUE),
               size = guide_legend(title.position="top"))+
        theme_bw()+
        theme(panel.border = element_blank(),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              axis.line = element_line(colour = "black"),
              axis.text.x =element_text(size=15,colour="black",vjust=0.7),
              axis.text.y =element_text(size=15,colour="black"),
              axis.title=element_text(size=15,colour="black", face="bold"))+
        theme(strip.background = element_rect(colour="white", fill="white",
                                              size=1.5, linetype="solid"))+
        theme(strip.text = element_text(size=15, color="black", face="bold.italic"))+ 
        theme(legend.text = element_text( size=15),
              legend.title= element_text(size=15, face="bold"))+   
        theme(plot.title = element_text(hjust = 0.5,size=20,face="bold"))
     
      # plots for incidence reduction
      lin_reg_incred <- function(x){
        reg <- lm(inc_red ~ Coverage, data = x)
        slope <- reg$coefficients[2]
        regs <- data.frame(cbind(x, slope))
        return(regs)
      }
      
      slopes_incred <-    LAI_specs %>%
        group_by(Efficacy, Halflife) %>%
        group_modify(~ lin_reg_incred(.x))
      
      colorsEfficacy<- brewer.pal( n=9, name = "YlOrBr")
      colorsHalflife <- brewer.pal( n=9, name = "RdPu")
      
      ggplot(data = slopes_incred) +
        geom_line(aes(x = Halflife, y = slope, color = as.factor(Efficacy*100)))+
        scale_colour_manual(values =colorsEfficacy[2:9])+
        
        labs(x = "Half-life of \nprotective efficacy [days]", y = "slope",
             color= "Initial protective \nefficacy [%]") +
        guides(colour = guide_legend(reverse = TRUE),
               size = guide_legend(title.position="top"))+
        theme_bw()+
        theme(panel.border = element_blank(),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              axis.line = element_line(colour = "black"),
              axis.text.x =element_text(size=15,colour="black",vjust=0.7),
              axis.text.y =element_text(size=15,colour="black"),
              axis.title=element_text(size=15,colour="black", face="bold"))+
        theme(strip.background = element_rect(colour="white", fill="white",
                                              size=1.5, linetype="solid"))+
        theme(strip.text = element_text(size=15, color="black", face="bold.italic"))+ 
        theme(legend.text = element_text( size=15),
              legend.title= element_text(size=15, face="bold"))+   
        theme(plot.title = element_text(hjust = 0.5,size=20,face="bold"))
      
      
      ggplot(data = slopes_incred) +
        geom_line(aes(x = Efficacy*100, y = slope, color = as.factor(Halflife)))+
        scale_colour_manual(values =colorsHalflife[2:9])+
        
        labs(x = "Initial protective \nefficacy [%]", y = "slope",
             color= "Half-life of protective \nefficacy [days]") +
        guides(colour = guide_legend(reverse = TRUE),
               size = guide_legend(title.position="top"))+
        theme_bw()+
        theme(panel.border = element_blank(),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              axis.line = element_line(colour = "black"),
              axis.text.x =element_text(size=15,colour="black",vjust=0.7),
              axis.text.y =element_text(size=15,colour="black"),
              axis.title=element_text(size=15,colour="black", face="bold"))+
        theme(strip.background = element_rect(colour="white", fill="white",
                                              size=1.5, linetype="solid"))+
        theme(strip.text = element_text(size=15, color="black", face="bold.italic"))+ 
        theme(legend.text = element_text( size=15),
              legend.title= element_text(size=15, face="bold"))+   
        theme(plot.title = element_text(hjust = 0.5,size=20,face="bold"))
      