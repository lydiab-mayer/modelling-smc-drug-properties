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

dir.create(file.path(paste0(mainDir,"/Outputs/E0_Cont/")), showWarnings = FALSE)




Access <- c(0.1,0.5)
LAI_dec <- c("exp")
seasonaility <- c("Sen","Mali")
settings <- expand.grid(Access,LAI_dec,seasonaility)

names(settings) <- c("Access","LAI_dec","seasonality")

points_high <-c(5,9,20,47,150)
  points_low <- c(3,4,8,28,150)
plotlist <- list()
for (i in 1:nrow(settings)) {
  setting <- settings[i,]
  agg <- read.table(paste("/scicore/home/penny/GROUP/smc_lai/E0_Cont/postprocessing_5/agg_E0_Cont_",setting[, "seasonality"],"_4.9167_",
                            setting[, "LAI_dec"],"_",setting[, "Access"],'.txt', sep = ""), 
                      header = T, as.is = TRUE, stringsAsFactors = FALSE)
  
  if(setting$Access==0.5) {
    points <- points_high
    sampled <- agg[which(agg$EIR %in% points),]
    sampledcpp <- sampled$pppy_y4_all
    
  } else {
    points <- points_low
    sampled <- agg[which(agg$EIR %in% points),]
    sampledcpp <- sampled$pppy_y4_all
  }
  
  print(sampled[order(sampled$EIR),c("EIR","pppy_y4_all","iprev_int_y5","prev_210_y5","Access","Seasonality")])
  
  p1 <- ggplot()+geom_line(data=agg, aes(EIR,pppy_y4_all,color="age 0-5"))+
    geom_point(data=sampled, aes(EIR,pppy_y4_all,color="sampled points"))+
    
   # geom_line(data=agg, aes(EIR,pppy_y4_210,color="age 2-10"))+
  #  geom_point(data=sampled, aes(EIR,pppy_y4_210,color="sampled points"))+
    
    labs(x="EIR",y="cumulative cases per person")+ylim(0,3)
  
 
  p2 <- ggplot()+geom_line(data=agg, aes(EIR,iprev_int_y5*100,color="age 0-5")) +
    geom_point(data=sampled, aes(EIR,iprev_int_y5*100,color="sampled points"))+
    
    geom_line(data=agg, aes(EIR,prev_210_y5*100,color="age 2-10"))+labs(x="EIR",y="prevalence")+ylim(0,100)+
  geom_point(data=sampled, aes(EIR,prev_210_y5*100,color="sampled points"))
    
  p <- ggarrange(p1,p2, common.legend = TRUE)
  
  ggsave(
    filename= paste0("./Outputs/E0_Cont/agg_E0_Cont_",setting[, "seasonality"],"_4.9167_",
                            setting[, "LAI_dec"],"_",setting[, "Access"],".jpg"),
    plot = last_plot(),
    
    width = 7,
    height = 5,
    dpi = 600 )
  
  if(setting$Access==0.5) {
    cpp <- sampledcpp
    ATC <- "high access to care"
  } else {
    cpp <- sampledcpp
    ATC <- "low access to care"
     }
  
  
  if(setting$seasonality=="Sen") {
    seas <- "short season"
  } else {
    seas <- "long season"
  }
  
  
  
  p3 <- ggplot()+geom_line(data=agg, aes(iprev_int_y5*100,pppy_y4_all,linetype="prev 0.25-5"),size=1.5) +
   # geom_line(data=agg, aes(iprev_int_y5*100,pppy_y4_210,color="cpp 2-10",linetype="prev 0.25-5"),size=1.5) +
    geom_line(data=agg, aes(prev_210_y5*100,pppy_y4_all,linetype="prev 2-10"),size=1.5) +
  #  geom_line(data=agg, aes(prev_210_y5*100,pppy_y4_210,color="cpp 2-10",linetype="prev 2-10"),size=1.5)+
  labs(x=expression(paste("Mean ", italic("Pf"), "PR", " [%]")),
       y="clinical incidence [0.25-5y] \n (events per person per year)")+
  geom_hline(yintercept = cpp, linetype="dotted", 
               color = "black", size=0.5)+
    scale_linetype_discrete(name = " ",
                          labels = c(expression(paste("Mean ", italic("Pf"), "PR"["0.25-5y"])),
                                                expression(paste("Mean ", italic("Pf"), "PR"["2-10y"]))))+

    theme(panel.border = element_blank(), panel.background = element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
          legend.text=element_text(size=15),
          axis.line = element_line(colour = "black"),
          axis.text=element_text(size=20,colour="black"),
          axis.text.x=element_text(margin = margin(t = 5)),
          axis.title=element_text(size=20,face="bold"), 
          legend.title=element_text(size=20,face="bold"))+ theme(legend.key=element_blank()) +
    theme(legend.background=element_blank()) + theme(legend.text = element_text(size=25,face="bold"))+
    theme(legend.key.width = unit(1.5,"cm"),legend.position = "right")+
    ggtitle(paste0(seas,"\n", ATC))+
   theme(plot.title = element_text(size=22,face="bold"))+xlim(0,80)+ylim(0,3.5)

  ggsave(
    filename= paste0("./Outputs/E0_Cont/prev_cpp_agg_E0_Cont_",setting[, "seasonality"],"_4.9167_",
                     setting[, "LAI_dec"],"_",setting[, "Access"],".jpg"),
    plot = last_plot(),
    
    width = 7,
    height = 5,
    dpi = 600 )
  plotlist[[i]] <- p3
  }


p <- ggarrange(plotlist = plotlist, common.legend = TRUE, legend="bottom")
ggsave(
  filename= paste0("./Outputs/E0_Cont/all_prev_cpp_agg_E0_Cont_high.jpg"),
  plot = last_plot(),
  
  width = 14,
  height = 10,
  dpi = 1200 )
  
for (i in 1:nrow(settings)) {
  
  setting <- settings[i,]

res <- import_EIRs(setting)

res2 <- data.frame(res %>% group_by(EIR) %>% summarise(diffmean = mean(diff) ))

if(setting$Access==0.5) {
  points <- points_high
  sampled <- res2[which(res2$EIR %in% points),]
} else {
  points <- points_low
  sampled <- res2[which(res2$EIR %in% points),]
  
}


p <- ggplot()+geom_point(data=res2, aes(EIR,diffmean)) +labs(x="EIR",y="peak shift")+
  geom_point(data=sampled, aes(EIR,diffmean,colour="sampled")) 


ggsave(
  filename= paste0("./Outputs/E0_Cont/peakshift_E0_Cont_",setting[, "seasonality"],"_4.9167_",
                   setting[, "LAI_dec"],"_",setting[, "Access"],".jpg"),
  plot = last_plot(),
  
  width = 5,
  height = 5,
  dpi = 600)

}