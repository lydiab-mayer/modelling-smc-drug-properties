
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

seasonality <- c("Mali","Sen")
LAI_dec <- c("exp","wei","hill")
IntAge <- c(4.9167)
Access <- c(0.1,0.5)


scenarios <- expand.grid(seasonality,LAI_dec,IntAge,Access)
colnames(scenarios) <- c("seasonality","LAI_dec","IntAge","Access")


points_low <-c(3,4,8,28,150)
points_high <- c(5,9,20,47,150)

setting <- scenarios[i,]

plot_optim_E4 <- function(setting){
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
    
    outfile <- paste("./Outputs/E4_Fig4/figure_4_", scen_name,".pdf", sep = "")
    
    
    
    if (file.exists(filename_UL) & file.exists(filename_optimisation)){
      print(scen_name)
      
      load(filename_UL)
      load(filename_CI)
      load(filename_pppy_SMC)
      load(filename_pppy_LAI)
      
      optim_results <- read.table(filename_optimisation, header = T)
      
      #plot left (heatmap)
      SMCcoverage_heatmap <- 0.6
      
      LAI_1 <- c(150, 1)   
      LAI_2 <- c(100, 0.9)
      LAI_3 <- c(90, 0.8)
      LAI_4 <- c(70, 0.7)
      
      
      n_lai_props <- 4
      
      my_colours <- c("gold1","yellowgreen", "turquoise4","darkgreen")
      index_1 <- which(optim_results$Coverage_SMC == SMCcoverage_heatmap & optim_results$Halflife == LAI_1[1] & optim_results$Efficacy == LAI_1[2])
      index_2 <- which(optim_results$Coverage_SMC == SMCcoverage_heatmap & optim_results$Halflife == LAI_2[1] & optim_results$Efficacy == LAI_2[2])
      index_3 <- which(optim_results$Coverage_SMC == SMCcoverage_heatmap & optim_results$Halflife == LAI_3[1] & optim_results$Efficacy == LAI_3[2])
      index_4 <- which(optim_results$Coverage_SMC == SMCcoverage_heatmap & optim_results$Halflife == LAI_4[1] & optim_results$Efficacy == LAI_4[2])
      
      dfheatEIR <- optim_results[optim_results$Coverage_SMC == SMCcoverage_heatmap,]
      dfheatEIR$EIR <- EIR
      
      dfsheat[[k]] <- dfheatEIR
      
      p11 <- ggplot(data = optim_results[optim_results$Coverage_SMC == SMCcoverage_heatmap,]) +
        geom_tile(aes(x = Efficacy*100, y = Halflife, fill = optimal_lai_coverage*100)) +
        geom_point(data = optim_results[index_1,], aes(x = Efficacy*100, y = Halflife), colour = my_colours[1] ,size=5)+
        geom_point(data = optim_results[index_2,], aes(x = Efficacy*100, y = Halflife), colour = my_colours[2] ,size=5 )+
        geom_point(data = optim_results[index_3,], aes(x = Efficacy*100, y = Halflife), colour = my_colours[3] ,size=5 )+
        geom_point(data = optim_results[index_4,], aes(x = Efficacy*100, y = Halflife), colour = my_colours[4] ,size=5)+
        scale_fill_gradientn(limits = c(40,100),
                             colours=c("#6d1c68",  "#f7921e"),
                             na.value = "grey") +
        labs(y = "Halflife of protective efficacy [d]", x = "Initial protective efficacy [%]",
             fill= "min. LAI coverage [%] \nto ensure non-inferiority") +
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
              legend.title= element_text(size=15, face="bold"))+   theme(legend.position="bottom")+
        
        theme(plot.title = element_text(hjust = 0.5,size=20,face="bold"))
      
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
      my_lab <- paste("Halflife: ", LAI_1[1], " d, initial efficacy: ",LAI_1[2]*100," %", sep = "")
      my_labs <- c(my_labs,my_lab)
      
      my_lab <- paste("Halflife: ", LAI_2[1], " d, initial efficacy: ",LAI_2[2]*100," %", sep = "")
      my_labs <- c(my_labs,my_lab)
      
      my_lab <- paste("Halflife: ", LAI_3[1], " d, initial efficacy: ",LAI_3[2]*100," %", sep = "")
      my_labs <- c(my_labs,my_lab)
      
      my_lab <- paste("Halflife: ", LAI_4[1], " d, initial efficacy: ",LAI_4[2]*100," %", sep = "")
      my_labs <- c(my_labs,my_lab)
      
      p12 <- ggplot(data = plot_LAI) +
        geom_line(aes(x = SMC*100, y = opt_LAI_cov*100, colour = as.factor(LAI_props)), size = 1.5) +
        geom_vline(aes(xintercept = 60), colour = "grey50", linetype = "dashed",size=1.5) +
        scale_colour_manual(values = my_colours, labels = my_labs)+xlim(40,100)+ylim(40,100)+
        labs( y = "min. LAI coverage [%]", x = "SMC-SP+AQ coverage [%]", colour = "LAI properties")+
        
        theme(axis.title = element_text(size=9)) +
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
        theme(legend.text = element_text( size=12,face="bold"),
              legend.title= element_text(size=13, face="bold"))+       theme(legend.position="bottom")+
        theme(plot.title = element_text(hjust = 0.5,size=20,face="bold"))+guides(color=guide_legend(ncol=2,nrow=2,byrow=TRUE,title.position="top"))
      
      
      
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
        geom_line(aes(x = SMC*100, y = cases_relative_difference, colour = as.factor(LAI_props)), size = 1.5) +
        geom_hline(aes(yintercept = 0),  colour = "grey50",size=1.5)+
        geom_vline(aes(xintercept = 60), colour = "grey50", linetype = "dashed",size=1.5) +
        scale_colour_manual(values = my_colours, labels = my_labs)+
        labs(y = expression(bold(paste("rel. diff. cpppy"["0.25-5y"]))),
             x = "SMC-SP+AQ coverage [%]", colour = "LAI properties")+
        theme(legend.position="none")+
        theme(axis.title = element_text(size=9)) +
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
        theme(legend.text = element_text( size=12,face="bold"),
              legend.title= element_text(size=13, face="bold"))+       theme(legend.position="bottom")+
        theme(plot.title = element_text(hjust = 0.5,size=20,face="bold"))+guides(color=guide_legend(ncol=2,nrow=2,byrow=TRUE,title.position="top"))
      
      
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
      
      
      regimen = ifelse(seasonality=="Mali","long season","short season")
      decay = ifelse(LAI_dec=="exp","exponential LAIs",
                     ifelse( LAI_dec=="hill","sigmoid LAIs","bi-phasic LAIs"))
      
      
      p14 <- ggplot(data = plot_LAI_cases_absolute) +
        geom_line(aes(x = LAI_coverage*100, y = cases_LAI, colour = as.factor(LAI_props)), size = 1.5) +
        geom_vline(aes(xintercept = opt_LAI_cov*100, colour = as.factor(LAI_props)), size = 1.5, linetype = "dashed") +
        scale_colour_manual(values = my_colours, labels = my_labs)+
        labs(y = expression(bold(paste("cpppy"["0.25-5y"]))),
             x = "LAI coverage [%]", colour = "LAI\nproperties") +
        xlim(c(40,100))+
        theme(legend.position="bottom")+
        theme(axis.title = element_text(size=9)) +
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
        theme(legend.text = element_text( size=12,face="bold"),
              legend.title= element_text(size=13, face="bold"))+      theme(legend.position="bottom")+
        theme(plot.title = element_text(hjust = 0.5,size=20,face="bold"))+guides(color=guide_legend(ncol=2,nrow=2,byrow=TRUE,title.position="top"))
      
      
      
      
      pleft <- ggarrange(p11+ guides(fill = guide_colourbar(barwidth = 13, barheight = 3)), labels = c("a"),
                         font.label = list(size = 40, face = "bold"),hjust=-0.1,vjust=1)
      pright <- ggarrange( plotlist = list(p14,p13,p12) , ncol = 1, nrow = 3,  heights = c(0.9,0.9,0.9), labels = c("b","c","d"),
                           font.label = list(size = 40, face = "bold"),hjust=0.8,vjust=1,common.legend = TRUE,legend="bottom")
      
      pall <- ggarrange(plotlist=list(
        ggparagraph(text=" ", face = "italic", size = 0.5, color = "black"),
        
        ggparagraph(text=" ", face = "italic", size = 0.5, color = "black"),
        pleft,               
        pright
      ),
      ncol = 2, nrow = 2, widths = c(1,1,1,1),heights = c(0.05,1))
      
    pall2 <-   annotate_figure(pall,
                      top = text_grob(        paste0(regimen,"; ",decay) ,

                                          color = "black", face = "bold", size = 35),
                      
      )
      
      
      ggsave(
        filename =  outfile,
        plot = last_plot(),
        width = 13,
        height = 9,
        dpi = 1200
      )  
    }
  }  
  
return(  rbindlist(dfsheat))
}



plotlist <- list()
for (k in 1:nrow(scenarios) ) {
  
  setting <- scenarios[k,]

  test <- plot_optim_E4(setting)
  test$Access <- setting[,"Access"]
  test$LAI_dec <- setting[,"LAI_dec"]
  test$seasonality <- setting[,"seasonality"]
  plotlist[[k]] <- test
}

alldf <- rbindlist(plotlist)
alldf$regimen = ifelse(alldf[, "seasonality"]=="Mali","long season","short season")
alldf$decay = ifelse(alldf[, "LAI_dec"]=="exp","exponential LAIs",
                  ifelse( alldf[, "LAI_dec"]=="hill","sigmoid LAIs","bi-phasic LAIs"))

save(alldf, file = paste0("./Outputs/E4_Fig4/heatmapsupp.RData"))

