
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


i=5
i=1

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
    
    outfile <- paste("./Outputs/E4_Fig4/figure_4_", scen_name,"_new.pdf", sep = "")
    
    
    
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
      optim_results$condition <- ifelse(optim_results$optimal_lai_coverage<0.6 & optim_results$optimal_lai_coverage>0, 1,0)
      
      p11 <- ggplot(data = optim_results[optim_results$Coverage_SMC == SMCcoverage_heatmap,]) +
        geom_tile(aes(x = Efficacy*100, y = Halflife, fill = optimal_lai_coverage*100)) +

        geom_rect(data=optim_results[which(optim_results$condition==1 & optim_results$Coverage_SMC==0.6),] ,
                   size=0.5, fill=NA, colour="lightblue2",
    aes(xmin=Efficacy*100 - 5, xmax=Efficacy*100 + 5, ymin=Halflife - 5, ymax=Halflife +5))+
      
        geom_point(data = optim_results[index_1,], aes(x = Efficacy*100, y = Halflife), colour = my_colours[1] ,size=7)+
        geom_point(data = optim_results[index_2,], aes(x = Efficacy*100, y = Halflife), colour = my_colours[2] ,size=7 )+
        geom_point(data = optim_results[index_3,], aes(x = Efficacy*100, y = Halflife), colour = my_colours[3] ,size=7 )+
        geom_point(data = optim_results[index_4,], aes(x = Efficacy*100, y = Halflife), colour = my_colours[4] ,size=7)+
        scale_fill_gradientn(limits = c(40,100),
                             colours=c("#6d1c68",  "#f7921e"),
                             na.value = "grey") +
        labs(y = "Half-life of protective efficacy [days]", x = "Initial protective efficacy [%]",
             fill= "minimal \nLAI coverage [%] \nto ensure non-inferiority") +
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
      
      # #plot right top (optimal LAI for different SMC)
      # SMCs <- unique(optim_results$Coverage_SMC)
      # 
      # plot_LAI <- data.frame(rep(SMCs, n_lai_props))
      # colnames(plot_LAI) <- c("SMC")
      # plot_LAI$HL <- rep(c(LAI_1[1],LAI_2[1],LAI_3[1],LAI_4[1]), each = length(SMCs)) 
      # plot_LAI$Ef <- rep(c(LAI_1[2],LAI_2[2],LAI_3[2],LAI_4[2]), each = length(SMCs))
      # 
      # plot_LAI$opt_LAI_cov <- 0
      # 
      # for (k in 1:nrow(plot_LAI)){
      #   index <- which(optim_results$Halflife == plot_LAI[k,]$HL &
      #                    optim_results$Efficacy == plot_LAI[k,]$Ef &
      #                    optim_results$Coverage_SMC == plot_LAI[k,]$SMC)
      #   plot_LAI[k,]$opt_LAI_cov <- optim_results[index, ]$optimal_lai_coverage
      # }
      # 
      # plot_LAI$LAI_props <- rep(seq(1,n_lai_props), each = nrow(plot_LAI)/n_lai_props)
      # plot_LAI[plot_LAI == -1] <- NA
      # 
      # my_labs <- c()
      # my_lab <- paste("Half-life: ", LAI_1[1], " d, initial efficacy: ",LAI_1[2]*100," %", sep = "")
      # my_labs <- c(my_labs,my_lab)
      # 
      # my_lab <- paste("Half-life: ", LAI_2[1], " d, initial efficacy: ",LAI_2[2]*100," %", sep = "")
      # my_labs <- c(my_labs,my_lab)
      # 
      # my_lab <- paste("Half-life: ", LAI_3[1], " d, initial efficacy: ",LAI_3[2]*100," %", sep = "")
      # my_labs <- c(my_labs,my_lab)
      # 
      # my_lab <- paste("Half-life: ", LAI_4[1], " d, initial efficacy: ",LAI_4[2]*100," %", sep = "")
      # my_labs <- c(my_labs,my_lab)
      # 
      # p12 <- ggplot(data = plot_LAI) +
      #   geom_line(aes(x = SMC*100, y = opt_LAI_cov*100, colour = as.factor(LAI_props)), size = 1.5) +
      #   scale_colour_manual(values = my_colours, labels = my_labs)+xlim(40,100)+ylim(40,100)+
      #   labs( y = "min. LAI coverage [%]", x = "SMC-SP+AQ coverage [%]", colour = "LAI properties")+
      #   
      #   theme(axis.title = element_text(size=9)) +
      #   theme_bw()+
      #   theme(panel.border = element_blank(),
      #         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      #         axis.line = element_line(colour = "black"),
      #         axis.text.x =element_text(size=15,colour="black",vjust=0.7),
      #         axis.text.y =element_text(size=15,colour="black"),
      #         axis.title=element_text(size=15,colour="black", face="bold"))+
      #   theme(strip.background = element_rect(colour="white", fill="white",
      #                                         size=1.5, linetype="solid"))+
      #   theme(strip.text = element_text(size=15, color="black", face="bold.italic"))+ 
      #   theme(legend.text = element_text( size=12,face="bold"),
      #         legend.title= element_text(size=13, face="bold"))+       theme(legend.position="bottom")+
      #   theme(plot.title = element_text(hjust = 0.5,size=20,face="bold"))+guides(color=guide_legend(ncol=2,nrow=2,byrow=TRUE,title.position="top"))
      # 
      # 
      # 
      # #plot right middle (relative difference cases per person per year with 20% above optimal coverage)
      # plot_LAI_cases <- plot_LAI
      # plot_LAI_cases$cases_SMC <- 0
      # plot_LAI_cases$cases_LAI <- 0
      # 
      # params <- data.frame(cbind(plot_LAI_cases$SMC, plot_LAI_cases$Ef,
      #                            plot_LAI_cases$HL, plot_LAI_cases$opt_LAI_cov))
      # 
      # colnames(params) <- c("Coverage_SMC", "Efficacy","Halflife","Coverage")
      # 
      # for (k in 1:nrow(params)){
      #   params[k,]$Coverage <- min(c(1,params[k,]$Coverage + 0.2*params[k,]$Coverage))
      # }
      # 
      # index <- which(is.na(params$Coverage) ==F)
      # 
      # cases_SMC <- predict(x = as.matrix(params[index,c(1)]), 
      #                      object = GP_trained_SMC)
      # 
      # cases_LAI <- predict(x = as.matrix(params[index,c(1,2,3,4)]), 
      #                      object = GP_trained_LAI)
      # plot_LAI_cases$cases_SMC <- NA
      # plot_LAI_cases$cases_LAI <- NA
      # 
      # plot_LAI_cases[index,]$cases_SMC <- cases_SMC$mean
      # plot_LAI_cases[index,]$cases_LAI <- cases_LAI$mean
      # 
      # plot_LAI_cases$cases_relative_difference <- (plot_LAI_cases$cases_SMC - plot_LAI_cases$cases_LAI)/plot_LAI_cases$cases_SMC 
      # 
      # p13 <- ggplot(data = plot_LAI_cases) +
      #   geom_line(aes(x = SMC*100, y = cases_relative_difference, colour = as.factor(LAI_props)), size = 1.5) +
      #   geom_hline(aes(yintercept = 0),  colour = "grey50",size=1.5)+
      #   scale_colour_manual(values = my_colours, labels = my_labs)+
      #   labs(y = expression(bold(paste("relative diff. cpppy"["0.25-5y"]))),
      #        x = "SMC-SP+AQ coverage [%]", colour = "LAI properties")+
      #   theme(legend.position="none")+
      #   theme(axis.title = element_text(size=9)) +
      #   theme_bw()+
      #   theme(panel.border = element_blank(),
      #         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      #         axis.line = element_line(colour = "black"),
      #         axis.text.x =element_text(size=15,colour="black",vjust=0.7),
      #         axis.text.y =element_text(size=15,colour="black"),
      #         axis.title=element_text(size=15,colour="black", face="bold"))+
      #   theme(strip.background = element_rect(colour="white", fill="white",
      #                                         size=1.5, linetype="solid"))+
      #   theme(strip.text = element_text(size=15, color="black", face="bold.italic"))+ 
      #   theme(legend.text = element_text( size=12,face="bold"),
      #         legend.title= element_text(size=13, face="bold"))+       theme(legend.position="bottom")+
      #   theme(plot.title = element_text(hjust = 0.5,size=20,face="bold"))+guides(color=guide_legend(ncol=2,nrow=2,byrow=TRUE,title.position="top"))
      # 
      # 
      # #plot right bottom
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
      # 
      # 
      # p14 <- ggplot(data = plot_LAI_cases_absolute) +
      #   geom_line(aes(x = LAI_coverage*100, y = cases_LAI, colour = as.factor(LAI_props)), size = 1.5) +
      #   geom_vline(aes(xintercept = opt_LAI_cov*100, colour = as.factor(LAI_props)), size = 1.5, linetype = "dashed") +
      #   scale_colour_manual(values = my_colours, labels = my_labs)+
      #   labs(y = expression(bold(paste("cpppy"["0.25-5y"]))),
      #        x = "LAI coverage [%]", colour = "LAI properties") +
      #   xlim(c(40,100))+
      #   theme(legend.position="bottom")+
      #   theme(axis.title = element_text(size=9)) +
      #   theme_bw()+
      #   theme(panel.border = element_blank(),
      #         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      #         axis.line = element_line(colour = "black"),
      #         axis.text.x =element_text(size=15,colour="black",vjust=0.7),
      #         axis.text.y =element_text(size=15,colour="black"),
      #         axis.title=element_text(size=15,colour="black", face="bold"))+
      #   theme(strip.background = element_rect(colour="white", fill="white",
      #                                         size=1.5, linetype="solid"))+
      #   theme(strip.text = element_text(size=15, color="black", face="bold.italic"))+ 
      #   theme(legend.text = element_text( size=12,face="bold"),
      #         legend.title= element_text(size=13, face="bold"))+      theme(legend.position="bottom")+
      #   theme(plot.title = element_text(hjust = 0.5,size=20,face="bold"))+guides(color=guide_legend(ncol=2,nrow=2,byrow=TRUE,title.position="top"))
      # 
      
      #plot right new
      
      #plot right middle (relative difference cases per person per year with 20% above optimal coverage)
      df_plot <- data.frame(seq(0.4,1,0.1))
      colnames(df_plot) <- c("LAI_coverage")
      df_plot$cases_SMC <- 0
      df_plot$cases_LAI_1 <- 0
      df_plot$cases_LAI_2 <- 0
      df_plot$cases_LAI_3 <- 0
      df_plot$cases_LAI_4 <- 0
      
      
      #get SMC cases at 60% coverage
      
      params <- data.frame(0.6)
      colnames(params) <- c("Coverage_SMC")
      cases_SMC <- predict(x = as.matrix(params), 
                           object = GP_trained_SMC)
      df_plot$cases_SMC <- cases_SMC$mean
      
      
      #get LAI cases properties 1
      params <- data.frame(cbind(rep(0.6,nrow(df_plot)), rep(LAI_1[2],nrow(df_plot)),rep(LAI_1[1],nrow(df_plot)),df_plot$LAI_coverage ))
      colnames(params) <- c("Coverage_SMC", "Efficacy","Halflife","Coverage")
      cases_LAI <- predict(x = as.matrix(params[,c(1,2,3,4)]), object = GP_trained_LAI)
      df_plot$cases_LAI_1 <- cases_LAI$mean
      
      #get LAI cases properties 2
      params <- data.frame(cbind(rep(0.6,nrow(df_plot)), rep(LAI_2[2],nrow(df_plot)),rep(LAI_2[1],nrow(df_plot)),df_plot$LAI_coverage ))
      colnames(params) <- c("Coverage_SMC", "Efficacy","Halflife","Coverage")
      cases_LAI <- predict(x = as.matrix(params[,c(1,2,3,4)]), object = GP_trained_LAI)
      df_plot$cases_LAI_2 <- cases_LAI$mean
      
      #get LAI cases properties 3
      params <- data.frame(cbind(rep(0.6,nrow(df_plot)), rep(LAI_3[2],nrow(df_plot)),rep(LAI_3[1],nrow(df_plot)),df_plot$LAI_coverage ))
      colnames(params) <- c("Coverage_SMC", "Efficacy","Halflife","Coverage")
      cases_LAI <- predict(x = as.matrix(params[,c(1,2,3,4)]), object = GP_trained_LAI)
      df_plot$cases_LAI_3 <- cases_LAI$mean
      
      #get LAI cases properties 4
      params <- data.frame(cbind(rep(0.6,nrow(df_plot)), rep(LAI_4[2],nrow(df_plot)),rep(LAI_4[1],nrow(df_plot)),df_plot$LAI_coverage ))
      colnames(params) <- c("Coverage_SMC", "Efficacy","Halflife","Coverage")
      cases_LAI <- predict(x = as.matrix(params[,c(1,2,3,4)]), object = GP_trained_LAI)
      df_plot$cases_LAI_4 <- cases_LAI$mean
      
      
      
      #get relative difference
      
      df_plot$rel_diff_1 <- (df_plot$cases_SMC - df_plot$cases_LAI_1)/df_plot$cases_SMC
      df_plot$rel_diff_2 <- (df_plot$cases_SMC - df_plot$cases_LAI_2)/df_plot$cases_SMC
      df_plot$rel_diff_3 <- (df_plot$cases_SMC - df_plot$cases_LAI_3)/df_plot$cases_SMC
      df_plot$rel_diff_4 <- (df_plot$cases_SMC - df_plot$cases_LAI_4)/df_plot$cases_SMC
      
      
      
      df <- melt(df_plot[,c(1,7,8,9,10)], id = "LAI_coverage")
      
      
      my_labs <- c()
      my_lab <- paste("Half-life: ", LAI_1[1], " days, initial efficacy: ",LAI_1[2]*100," %", sep = "")
      my_labs <- c(my_labs,my_lab)

      my_lab <- paste("Half-life: ", LAI_2[1], " days, initial efficacy: ",LAI_2[2]*100," %", sep = "")
      my_labs <- c(my_labs,my_lab)

      my_lab <- paste("Half-life: ", LAI_3[1], " days, initial efficacy: ",LAI_3[2]*100," %", sep = "")
      my_labs <- c(my_labs,my_lab)

      my_lab <- paste("Half-life: ", LAI_4[1], " days, initial efficacy: ",LAI_4[2]*100," %", sep = "")
      my_labs <- c(my_labs,my_lab)
      
      optLAIcoverages <- optim_results[c(index_1,index_2,index_3,index_4),5]
      
     # df$optLAIcoverage <- rep(optLAIcoverages,each=7)
      
      df_plot$optLAIcoverages_1 <- rep(optLAIcoverages[1],7)
      df_plot$optLAIcoverages_2 <- rep(optLAIcoverages[2],7)
      df_plot$optLAIcoverages_3 <- rep(optLAIcoverages[3],7)
      df_plot$optLAIcoverages_4 <- rep(optLAIcoverages[4],7)
      
      my_colours <- c("gold1","yellowgreen", "turquoise4","darkgreen")
      
      
      p151 <- ggplot(data = df_plot) +
        # geom_vline( aes(xintercept = optLAIcoverages_1*100), size = 1.5, linetype = "dashed", color=my_colours[4]) +
        
        geom_rect(data=NULL,aes(xmin=40,xmax=optLAIcoverages_1*100,ymin=-0.5,ymax=0),fill="gray")+
        geom_rect(data=NULL,aes(xmin=optLAIcoverages_1*100,xmax=60,ymin=-0.5,ymax=df_plot[3,"rel_diff_1"] ),fill="lightblue2")+
        
        geom_line(aes(x = LAI_coverage*100, y = rel_diff_1, colour="Half-life: 150 days, initial efficacy: 100 %"), size = 1) +
        scale_colour_manual(values = my_colours[1], labels = my_labs[1])+
        geom_vline(aes(xintercept = 60), linetype = "dashed") +
        geom_hline(aes(yintercept = 0), linetype = "dashed") +
        labs(y = expression(bold(paste("relative difference \ncases per person per year"["0.25-5y"]))),
             x = "LAI coverage [%]", colour = "LAI properties") +
        xlim(c(37,100))+ylim(c(-0.5,0.75))+
        theme(legend.position="bottom")+
        # geom_segment( aes(x = 65, y =0.72, xend =  75, yend = 0.72),
        #               colour = 'black',size=0.5,linejoin = c('bevel'),
        #               arrow = arrow(length = unit(0.1, 'inches'))) +
        # geom_segment( aes(x = 55, y =0.72, xend =45, yend = 0.72),
        #               colour = 'black',size=0.5,linejoin = c('bevel'),
        #               arrow = arrow(length = unit(0.1, 'inches')))+
        # 
        # geom_segment( aes(x = 38, y =-0.05, xend =  38, yend = -0.15),
        #               colour = 'black',size=0.5,linejoin = c('bevel'),
        #               arrow = arrow(length = unit(0.1, 'inches'))) +
        # geom_segment( aes(x = 38, y =0.05, xend =38, yend = 0.15),
      #               colour = 'black',size=0.5,linejoin = c('bevel'),
      #               arrow = arrow(length = unit(0.1, 'inches')))+
      #annotate("text", x = 61, y = 0.6 , label = "Higher LAI \ncoverage",size=4,
        #        hjust = 0, fontface =2)+
        #annotate("text", x = 59, y = 0.6 , label = "Higher SMC \ncoverage",size=4,
        #        hjust = 1, fontface =2)+
        #  annotate("text", x = 40, y = 0.15 , label = "More cases \nwith SMC",size=4,
        #           hjust = 0, fontface =2)+
        # annotate("text", x = 40, y = -0.15 , label = "More cases \nwith LAIs",size=4,
        #         hjust = 0, fontface =2)+
        theme_bw()+
        theme(panel.border = element_blank(),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              axis.line = element_line(colour = "black"),
              axis.text.x =element_text(size=13,colour="black",vjust=0.7),
              axis.text.y =element_text(size=13,colour="black"),
              axis.title=element_text(size=13,colour="black", face="bold"))+
        theme(strip.background = element_rect(colour="white", fill="white",
                                              size=1.5, linetype="solid"))+
        theme(strip.text = element_text(size=15, color="black", face="bold.italic"))+
        theme(legend.text = element_text( size=11,face="bold"),
              legend.title= element_text(size=13, face="bold"))+      theme(legend.position="bottom")+
        theme(plot.title = element_text(hjust = 0.5,size=20,face="bold"))+guides(color=guide_legend(ncol=2,nrow=2,byrow=TRUE,title.position="top"))
      
      p152 <- ggplot(data = df_plot) +
        # geom_vline( aes(xintercept = optLAIcoverages_4*100), size = 1.5, linetype = "dashed", color=my_colours[4]) +
        
        geom_rect(data=NULL,aes(xmin=40,xmax=optLAIcoverages_2*100,ymin=-0.5,ymax=0),fill="gray")+
        
        geom_line(aes(x = LAI_coverage*100, y = rel_diff_2, colour="Half-life: 100 days, initial efficacy: 90 %"), size = 1) +
        scale_colour_manual(values = my_colours[2], labels = my_labs[2])+
        geom_vline(aes(xintercept = 60), linetype = "dashed") +
        geom_hline(aes(yintercept = 0), linetype = "dashed") +
        labs(y = expression(bold(paste("relative difference \ncases per person per year"["0.25-5y"]))),
             x = "LAI coverage [%]", colour = "LAI properties") +
        xlim(c(37,100))+ylim(c(-0.5,0.75))+
        theme(legend.position="bottom")+
        # geom_segment( aes(x = 65, y =0.72, xend =  75, yend = 0.72),
        #               colour = 'black',size=0.5,linejoin = c('bevel'),
        #               arrow = arrow(length = unit(0.1, 'inches'))) +
        # geom_segment( aes(x = 55, y =0.72, xend =45, yend = 0.72),
        #               colour = 'black',size=0.5,linejoin = c('bevel'),
        #               arrow = arrow(length = unit(0.1, 'inches')))+
        # 
        # geom_segment( aes(x = 38, y =-0.05, xend =  38, yend = -0.15),
        #               colour = 'black',size=0.5,linejoin = c('bevel'),
        #               arrow = arrow(length = unit(0.1, 'inches'))) +
        # geom_segment( aes(x = 38, y =0.05, xend =38, yend = 0.15),
      #               colour = 'black',size=0.5,linejoin = c('bevel'),
      #               arrow = arrow(length = unit(0.1, 'inches')))+
      #annotate("text", x = 61, y = 0.6 , label = "Higher LAI \ncoverage",size=4,
        #        hjust = 0, fontface =2)+
        #annotate("text", x = 59, y = 0.6 , label = "Higher SMC \ncoverage",size=4,
        #        hjust = 1, fontface =2)+
        #  annotate("text", x = 40, y = 0.15 , label = "More cases \nwith SMC",size=4,
        #           hjust = 0, fontface =2)+
        # annotate("text", x = 40, y = -0.15 , label = "More cases \nwith LAIs",size=4,
        #         hjust = 0, fontface =2)+
        theme_bw()+
        theme(panel.border = element_blank(),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              axis.line = element_line(colour = "black"),
              axis.text.x =element_text(size=13,colour="black",vjust=0.7),
              axis.text.y =element_text(size=13,colour="black"),
              axis.title=element_text(size=13,colour="black", face="bold"))+
        theme(strip.background = element_rect(colour="white", fill="white",
                                              size=1.5, linetype="solid"))+
        theme(strip.text = element_text(size=15, color="black", face="bold.italic"))+
        theme(legend.text = element_text( size=11,face="bold"),
              legend.title= element_text(size=13, face="bold"))+      theme(legend.position="bottom")+
        theme(plot.title = element_text(hjust = 0.5,size=20,face="bold"))+guides(color=guide_legend(ncol=2,nrow=2,byrow=TRUE,title.position="top"))
      
      p153 <- ggplot(data = df_plot) +
        # geom_vline( aes(xintercept = optLAIcoverages_3*100), size = 1.5, linetype = "dashed", color=my_colours[4]) +
        
        geom_rect(data=NULL,aes(xmin=40,xmax=optLAIcoverages_3*100,ymin=-0.5,ymax=0),fill="gray")+
        
        geom_line(aes(x = LAI_coverage*100, y = rel_diff_3, colour="Half-life: 90 days, initial efficacy: 80 %"), size = 1) +
        scale_colour_manual(values = my_colours[3], labels = my_labs[3])+
        geom_vline(aes(xintercept = 60), linetype = "dashed") +
        geom_hline(aes(yintercept = 0), linetype = "dashed") +
        labs(y = expression(bold(paste("relative difference \ncases per person per year"["0.25-5y"]))),
             x = "LAI coverage [%]", colour = "LAI properties") +
        xlim(c(37,100))+ylim(c(-0.5,0.75))+
        theme(legend.position="bottom")+
        # geom_segment( aes(x = 65, y =0.72, xend =  75, yend = 0.72),
        #               colour = 'black',size=0.5,linejoin = c('bevel'),
        #               arrow = arrow(length = unit(0.1, 'inches'))) +
        # geom_segment( aes(x = 55, y =0.72, xend =45, yend = 0.72),
        #               colour = 'black',size=0.5,linejoin = c('bevel'),
        #               arrow = arrow(length = unit(0.1, 'inches')))+
        # 
        # geom_segment( aes(x = 38, y =-0.05, xend =  38, yend = -0.15),
        #               colour = 'black',size=0.5,linejoin = c('bevel'),
        #               arrow = arrow(length = unit(0.1, 'inches'))) +
        # geom_segment( aes(x = 38, y =0.05, xend =38, yend = 0.15),
      #               colour = 'black',size=0.5,linejoin = c('bevel'),
      #               arrow = arrow(length = unit(0.1, 'inches')))+
      #annotate("text", x = 61, y = 0.6 , label = "Higher LAI \ncoverage",size=4,
        #        hjust = 0, fontface =2)+
        #annotate("text", x = 59, y = 0.6 , label = "Higher SMC \ncoverage",size=4,
        #        hjust = 1, fontface =2)+
        #  annotate("text", x = 40, y = 0.15 , label = "More cases \nwith SMC",size=4,
        #           hjust = 0, fontface =2)+
        # annotate("text", x = 40, y = -0.15 , label = "More cases \nwith LAIs",size=4,
        #         hjust = 0, fontface =2)+
        theme_bw()+
        theme(panel.border = element_blank(),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              axis.line = element_line(colour = "black"),
              axis.text.x =element_text(size=13,colour="black",vjust=0.7),
              axis.text.y =element_text(size=13,colour="black"),
              axis.title=element_text(size=13,colour="black", face="bold"))+
        theme(strip.background = element_rect(colour="white", fill="white",
                                              size=1.5, linetype="solid"))+
        theme(strip.text = element_text(size=15, color="black", face="bold.italic"))+
        theme(legend.text = element_text( size=11,face="bold"),
              legend.title= element_text(size=13, face="bold"))+      theme(legend.position="bottom")+
        theme(plot.title = element_text(hjust = 0.5,size=20,face="bold"))+guides(color=guide_legend(ncol=2,nrow=2,byrow=TRUE,title.position="top"))
      
      p154 <- ggplot(data = df_plot) +
       # geom_vline( aes(xintercept = optLAIcoverages_4*100), size = 1.5, linetype = "dashed", color=my_colours[4]) +
        
        geom_rect(data=NULL,aes(xmin=40,xmax=optLAIcoverages_4*100,ymin=-0.5,ymax=0),fill="gray")+
        geom_line(aes(x = LAI_coverage*100, y = rel_diff_4, colour="Half-life: 70 days, initial efficacy: 70 %"), size = 1) +
        scale_colour_manual(values = my_colours[4], labels = my_labs[4])+
        geom_vline(aes(xintercept = 60), linetype = "dashed") +
        geom_hline(aes(yintercept = 0), linetype = "dashed") +
        labs(y = expression(bold(paste("relative difference \ncases per person per year"["0.25-5y"]))),
             x = "LAI coverage [%]", colour = "LAI properties") +
        xlim(c(37,100))+ylim(c(-0.5,0.75))+
        theme(legend.position="bottom")+
        # geom_segment( aes(x = 65, y =0.72, xend =  75, yend = 0.72),
        #               colour = 'black',size=0.5,linejoin = c('bevel'),
        #               arrow = arrow(length = unit(0.1, 'inches'))) +
        # geom_segment( aes(x = 55, y =0.72, xend =45, yend = 0.72),
        #               colour = 'black',size=0.5,linejoin = c('bevel'),
        #               arrow = arrow(length = unit(0.1, 'inches')))+
        # 
        # geom_segment( aes(x = 38, y =-0.05, xend =  38, yend = -0.15),
        #               colour = 'black',size=0.5,linejoin = c('bevel'),
        #               arrow = arrow(length = unit(0.1, 'inches'))) +
        # geom_segment( aes(x = 38, y =0.05, xend =38, yend = 0.15),
        #               colour = 'black',size=0.5,linejoin = c('bevel'),
        #               arrow = arrow(length = unit(0.1, 'inches')))+
        #annotate("text", x = 61, y = 0.6 , label = "Higher LAI \ncoverage",size=4,
         #        hjust = 0, fontface =2)+
        #annotate("text", x = 59, y = 0.6 , label = "Higher SMC \ncoverage",size=4,
         #        hjust = 1, fontface =2)+
      #  annotate("text", x = 40, y = 0.15 , label = "More cases \nwith SMC",size=4,
      #           hjust = 0, fontface =2)+
       # annotate("text", x = 40, y = -0.15 , label = "More cases \nwith LAIs",size=4,
        #         hjust = 0, fontface =2)+
        theme_bw()+
        theme(panel.border = element_blank(),
              panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              axis.line = element_line(colour = "black"),
              axis.text.x =element_text(size=13,colour="black",vjust=0.7),
              axis.text.y =element_text(size=13,colour="black"),
              axis.title=element_text(size=13,colour="black", face="bold"))+
        theme(strip.background = element_rect(colour="white", fill="white",
                                              size=1.5, linetype="solid"))+
        theme(strip.text = element_text(size=15, color="black", face="bold.italic"))+
        theme(legend.text = element_text( size=11,face="bold"),
              legend.title= element_text(size=13, face="bold"))+      theme(legend.position="bottom")+
        theme(plot.title = element_text(hjust = 0.5,size=20,face="bold"))+guides(color=guide_legend(ncol=2,nrow=2,byrow=TRUE,title.position="top"))
      
      
      
      # p15 <- ggplot(data = df) +
      #   geom_line(aes(x = LAI_coverage*100, y = value, colour = as.factor(variable)), size = 1) +
      #   scale_colour_manual(values = my_colours, labels = my_labs)+
      #   geom_vline(aes(xintercept = 60), linetype = "dashed") +
      #   geom_hline(aes(yintercept = 0), linetype = "dashed") +
      #   geom_vline( aes(xintercept = optLAIcoverage*100, colour = as.factor(variable)), size = 1.5, linetype = "dashed") +
      #   labs(y = expression(bold(paste("relative difference cpppy"["0.25-5y"]))),
      #        x = "LAI coverage [%]", colour = "LAI\nproperties") +
      #   xlim(c(37,100))+
      #   theme(legend.position="bottom")+
      #   geom_segment( aes(x = 65, y =0.72, xend =  75, yend = 0.72),
      #                 colour = 'black',size=0.5,linejoin = c('bevel'),
      #                 arrow = arrow(length = unit(0.1, 'inches'))) +
      #   geom_segment( aes(x = 55, y =0.72, xend =45, yend = 0.72),
      #                 colour = 'black',size=0.5,linejoin = c('bevel'),
      #                 arrow = arrow(length = unit(0.1, 'inches')))+
      #   
      #   geom_segment( aes(x = 38, y =-0.05, xend =  38, yend = -0.15),
      #                 colour = 'black',size=0.5,linejoin = c('bevel'),
      #                 arrow = arrow(length = unit(0.1, 'inches'))) +
      #   geom_segment( aes(x = 38, y =0.05, xend =38, yend = 0.15),
      #                 colour = 'black',size=0.5,linejoin = c('bevel'),
      #                 arrow = arrow(length = unit(0.1, 'inches')))+
      #   annotate("text", x = 61, y = 0.6 , label = "LAI coverage \nhigher than \nSMC coverage",size=5,
      #            hjust = 0, fontface =2)+
      #   annotate("text", x = 59, y = 0.6 , label = "LAI coverage \nlower than \nSMC coverage",size=5,
      #            hjust = 1, fontface =2)+
      #   annotate("text", x = 40, y = 0.07 , label = "More cases \nwith SMC-SP+AQ",size=5,
      #            hjust = 0, fontface =2)+
      #   annotate("text", x = 40, y = -0.07 , label = "More cases \nwith LAIs",size=5,
      #            hjust = 0, fontface =2)+
      #   theme_bw()+
      #   theme(panel.border = element_blank(),
      #         panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      #         axis.line = element_line(colour = "black"),
      #         axis.text.x =element_text(size=15,colour="black",vjust=0.7),
      #         axis.text.y =element_text(size=15,colour="black"),
      #         axis.title=element_text(size=15,colour="black", face="bold"))+
      #   theme(strip.background = element_rect(colour="white", fill="white",
      #                                         size=1.5, linetype="solid"))+
      #   theme(strip.text = element_text(size=15, color="black", face="bold.italic"))+
      #   theme(legend.text = element_text( size=12,face="bold"),
      #         legend.title= element_text(size=13, face="bold"))+      theme(legend.position="bottom")+
      #   theme(plot.title = element_text(hjust = 0.5,size=20,face="bold"))+guides(color=guide_legend(ncol=2,nrow=2,byrow=TRUE,title.position="top"))

        
      
      
      pleft <- ggarrange(p11+ guides(fill = guide_colourbar(barwidth = 13, barheight = 3)), labels = c("a"),
                         font.label = list(size = 40, face = "bold"),hjust=-0.1,vjust=1)
      # pright <- ggarrange( plotlist = list(p14,p13,p12) , ncol = 1, nrow = 3,  heights = c(0.9,0.9,0.9), labels = c("b","c","d"),
      #                      font.label = list(size = 40, face = "bold"),hjust=0.8,vjust=1,common.legend = TRUE,legend="bottom")
      
      pright <- ggarrange( plotlist = list(p151,  
                                           ggparagraph(text=" ", face = "italic", size = 0.02, color = "black"),   
                                           p152,p153,
                                           ggparagraph(text=" ", face = "italic", size = 0.02, color = "black"),   
                                           p154) , ncol = 3, nrow = 2, labels = c("b","","c","d","","e"),
                           font.label = list(size = 40, face = "bold"),
                           hjust=0.8,vjust=1, widths = c(1,0.05,1))
      
      pall <- ggarrange(plotlist=list(
        ggparagraph(text=" ", face = "italic", size = 0.05, color = "black"),       
        ggparagraph(text=" ", face = "italic", size = 0.05, color = "black"),       
        ggparagraph(text=" ", face = "italic", size = 0.05, color = "black"),       
        ggparagraph(text=" ", face = "italic", size = 0.05, color = "black"),       
        pleft,      ggparagraph(text=" ", face = "italic", size = 0.05, color = "black"),         
        pright,        ggparagraph(text=" ", face = "italic", size = 0.05, color = "black")

      ),
      ncol = 4, nrow = 2, widths = c(1,0.1,1.3,0.1),heights = c(0.05,1))
      
    #pall2 <-   annotate_figure(pall,
    #                  top = text_grob(        paste0(regimen,"; ",decay) ,

    #                                      color = "black", face = "bold", size = 35),
                      
#      )
      
      
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

toplot <- c(1,5)

plotlist <- list()
for (k in 1:length(toplot) ) {
  
  setting <- scenarios[toplot[k],]

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

