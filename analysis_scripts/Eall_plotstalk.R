# plot non-inferiority results
library(ggpubr)
library(dplyr)
library(ggplot2)
library(cowplot)   
library(patchwork) 

rm(list = ls())

mainDir  <- "~/smc_lai/analysis_workflow/analysis_scripts"
setwd(mainDir)
source('~/smc_lai/analysis_workflow/3_GP_train/GP_toolbox.R')
source('~/smc_lai/analysis_workflow/analysis_scripts/supp/E5_resources.R')



library(hetGP)
library(ggplot2)
library(ggpubr)
library(reshape2)
library(tgp)

SMC_HL <- c(31.295)
LAI_dec <- c("hill")
seasonaility <- c("Sen")
setting<- expand.grid(SMC_HL,LAI_dec,seasonaility)
names(setting) <- c("SMC_HL","LAI_dec","seasonality")

plot <- pred_GP_E5_noninf_talk(setting)

q2 <- ggarrange(plotlist = list(plot,ggparagraph(text=" ", face = "italic", size = 0.5, color = "black")
) ,
heights=c(3,3), 
widths=c(1,0.1
),
nrow=1,ncol=2)

ggsave(
  filename =  paste0("./Outputs/Eall_talk/noninf_normSMC_Sen.jpg"),
  plot = last_plot(),
  
  width = 7,
  height = 6,
  dpi = 600)



SMC_HL <- c(31.295)
LAI_dec <- c("hill")
seasonaility <- c("Mali")
setting<- expand.grid(SMC_HL,LAI_dec,seasonaility)
names(setting) <- c("SMC_HL","LAI_dec","seasonality")

  plot <- pred_GP_E5_noninf_talk(setting)

  q2 <- ggarrange(plotlist = list(plot,ggparagraph(text=" ", face = "italic", size = 0.5, color = "black")
                                ) ,
                  heights=c(3,3), 
                  widths=c(1,0.1
                  ),
                  nrow=1,ncol=2)

ggsave(
  filename =  paste0("./Outputs/Eall_talk/noninf_normSMC.jpg"),
  plot = last_plot(),
  
  width = 7,
  height = 6,
  dpi = 600)




plot2 <- pred_GP_E5_incred_talk(setting)



q2 <- ggarrange(plotlist = list(plot2[[1]],ggparagraph(text=" ", face = "italic", size = 0.5, color = "black")
) ,
heights=c(3,3), 
widths=c(1,0.1
),
nrow=1,ncol=2)

ggsave(
  filename =  paste0("./Outputs/Eall_talk/inc_red_normSMCtest.jpg"),
  plot = q2,
  
  width = 6,
  height = 6,
  dpi = 600,
  
)

leg <- plot2[[2]]

ggsave(
  filename =  paste0("./Outputs/Eall_talk/inc_red_normleg.jpg"),
  plot = leg,
  
  width = 6,
  height = 4,
  dpi = 600,
  
)


seasonality <- c("Mali")
LAI_dec <- c("hill")
IntAge <- c(4.9167)
Access <- c(0.5)


setting <- expand.grid(seasonality,LAI_dec,IntAge,Access)
colnames(setting) <- c("seasonality","LAI_dec","IntAge","Access")




points_high <- c(20)


EIR <- c(20)

    seasonality <- setting[,1]
    decay <- setting[,2]
    age <- setting[,3]
    acess <- setting[,4]
    scen_name <- paste(seasonality,"_",age,"_",decay,"_",acess,"_",EIR, sep = "")
    
    filename_UL <- paste('/scicore/home/penny/GROUP/smc_lai/E3E4disc_comparison/gp_trained/gp_trained_UL_', scen_name, '.RData', sep = "")
    filename_CI <- paste('/scicore/home/penny/GROUP/smc_lai/E3E4disc_comparison/gp_trained/gp_trained_CI_', scen_name, '.RData', sep = "")
    filename_pppy_SMC <- paste('/scicore/home/penny/GROUP/smc_lai/E3E4disc_comparison/gp_trained/gp_trained_pppy_SMC_', scen_name, '.RData', sep = "")
    filename_pppy_LAI <- paste('/scicore/home/penny/GROUP/smc_lai/E3E4disc_comparison/gp_trained/gp_trained_pppy_LAI_', scen_name, '.RData', sep = "")
    
    filename_optimisation <- paste('/scicore/home/penny/GROUP/smc_lai/E3E4disc_comparison/outfiles_optimisation_summary/scenarios_all_',scen_name,".txt",sep = "")
    
    outfile <- paste("./Outputs/E4_Fig4/figure_4_", scen_name,".pdf", sep = "")
    
    
    
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
      #  geom_point(data = optim_results[index_1,], aes(x = Efficacy*100, y = Halflife), colour = my_colours[1] ,size=7)+
      #  geom_point(data = optim_results[index_2,], aes(x = Efficacy*100, y = Halflife), colour = my_colours[2] ,size=7 )+
       # geom_point(data = optim_results[index_3,], aes(x = Efficacy*100, y = Halflife), colour = my_colours[3] ,size=7 )+
       # geom_point(data = optim_results[index_4,], aes(x = Efficacy*100, y = Halflife), colour = my_colours[4] ,size=7)+
        scale_fill_gradientn(limits = c(40,100),
                             colours=c("#6d1c68",  "#f7921e"),
                             na.value = "grey") +
        labs(y = "Half-life of protective efficacy [d]", x = "Initial protective efficacy [%]",
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
        theme(legend.text = element_text( size=10),
              legend.title= element_text(size=15, face="bold"))+   theme(legend.position="bottom")+
        
        theme(plot.title = element_text(hjust = 0.5,size=20,face="bold"))
      
   
      ggsave(
        filename =  paste0("./Outputs/Eall_talk/opt.jpg"),
        plot = p11,
        
        width = 7,
        height = 6,
        dpi = 600,
        
      )

