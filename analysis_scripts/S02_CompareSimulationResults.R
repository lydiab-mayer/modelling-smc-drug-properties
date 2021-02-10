#=========================================
#
#	Analysis of OM experiments smc_lai
# -plotting cases averted in same setting example plot
# 
#	author:  Lydia Burgert
#			
# input: 
# output: 
#=========================================
library(dplyr)
library(ggrepel)

rm(list = ls())

source("./supp/import_functions.R")
setwd("~/smc_lai/analysis_workflow/analysis_scripts")
sim_id <- 1

# import Results to compare
SMC1_low <- import_cpppy(Exp="E3_SMCSMCdisc",scen_ids=seq(6651,6655),param_table_file="seeds_E3SMCSMCdisc_Mali_4.9167_hill_0.1_4.txt")
SMC2_low <- import_cpppy(Exp="E3_SMCSMCdisc",scen_ids=seq(526,530),param_table_file="seeds_E3SMCSMCdisc_Mali_4.9167_hill_0.1_4.txt")
LAI1_low <- import_cpppy(Exp="E4_SMCLAIdisc",scen_ids=seq(6651,6655),param_table_file="seeds_E4SMCLAIdisc_Mali_4.9167_hill_0.1_4.txt")
LAI2_low <- import_cpppy(Exp="E4_SMCLAIdisc",scen_ids=seq(526,530),param_table_file="seeds_E4SMCLAIdisc_Mali_4.9167_hill_0.1_4.txt")

SMC1_high <- import_cpppy(Exp="E3_SMCSMCdisc",scen_ids=seq(6651,6655),param_table_file="seeds_E3SMCSMCdisc_Mali_4.9167_hill_0.1_150.txt")
SMC2_high <- import_cpppy(Exp="E3_SMCSMCdisc",scen_ids=seq(526,530),param_table_file="seeds_E3SMCSMCdisc_Mali_4.9167_hill_0.1_150.txt")
LAI1_high <- import_cpppy(Exp="E4_SMCLAIdisc",scen_ids=seq(6651,6655),param_table_file="seeds_E4SMCLAIdisc_Mali_4.9167_hill_0.1_150.txt")
LAI2_high <- import_cpppy(Exp="E4_SMCLAIdisc",scen_ids=seq(526,530),param_table_file="seeds_E4SMCLAIdisc_Mali_4.9167_hill_0.1_150.txt")

cont_high <- import_cpppy(Exp="E0_Cont",scen_ids=c(691, 692, 693, 694, 695),param_table_file="seeds_E0_Cont_Mali_4.9167_exp_0.1.txt")
  cont_low <- import_cpppy(Exp="E0_Cont",scen_ids=c(46 ,47, 48, 49, 50),param_table_file="seeds_E0_Cont_Mali_4.9167_exp_0.1.txt")

  
  colors <- c("LAI1" = "#517594", "SMC1"="#42733C","cont"="black")
types <- c("high" = "solid", "low"="dashed")

p1 <- ggplot()+
  geom_line(data=cont_low[[1]], aes(x=trialtime*5/30,y=cppcum,color="cont"),size=1)+
  
 geom_line(data=SMC1_low[[1]], aes(x=trialtime*5/30,y=cppcum,color="SMC1",linetype="high"),size=1)+
  geom_line(data=SMC2_low[[1]], aes(x=trialtime*5/30,y=cppcum,color="SMC1",linetype="low"),size=1)+
  geom_line(data=LAI1_low[[1]], aes(x=trialtime*5/30,y=cppcum,color="LAI1",linetype="high"),size=1)+
  geom_line(data=LAI2_low[[1]], aes(x=trialtime*5/30,y=cppcum,color="LAI1",linetype="low"),size=1)+
  scale_linetype_manual(values = types) +
  scale_color_manual(values = colors,
                     labels=c("Control","LAI",
                              "SMC-SP+AQ"))+
  labs( x ="Month", y = expression(paste("cpp"["0.25-5y"])),
        color = "Intervention", linetype="Coverage")+
  scale_x_continuous(breaks=c(0,1,2,3,4,5,6,7,8,9,10,11), 
                     labels=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"), 
                     limits = c(0,11))+
  guides(color=guide_legend(nrow=3,ncol=1,byrow=TRUE,title.position = "top"))+
  
  guides(linetype=guide_legend(nrow=2,ncol=1,byrow=TRUE,title.position = "top"))+
  theme(panel.border = element_blank(), panel.background = element_blank(),title = element_text(size=15, face="bold"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.text=element_text(size=13),
        axis.line = element_line(colour = "black"),
        axis.text=element_text(size=13,colour="black"),
        axis.text.x=element_text(margin = margin(t = 5)),
        axis.title=element_text(size=13,face="bold"),  
        legend.title=element_text(size=13))+ theme(legend.position = "bottom") +
  theme(legend.background=element_blank()) +ggtitle(paste("low initial clinical incidence"))+
  theme(legend.key.width = unit(2.5,"cm"))
  



p2 <- ggplot()+
  geom_line(data=cont_high[[1]], aes(x=trialtime*5/30,y=cppcum,color="cont"),size=1)+
  
  geom_line(data=SMC1_high[[1]], aes(x=trialtime*5/30,y=cppcum,color="SMC1",linetype="high"),size=1)+
  geom_line(data=SMC2_high[[1]], aes(x=trialtime*5/30,y=cppcum,color="SMC1",linetype="low"),size=1)+
  geom_line(data=LAI1_high[[1]], aes(x=trialtime*5/30,y=cppcum,color="LAI1",linetype="high"),size=1)+
  geom_line(data=LAI2_high[[1]], aes(x=trialtime*5/30,y=cppcum,color="LAI1",linetype="low"),size=1)+
  scale_linetype_manual(values = types) +
  scale_color_manual(values = colors,
                     labels=c("Control","LAI",
                              "SMC-SP+AQ"))+
  labs( x ="Month", y = expression(paste("cpp"["0.25-5y"])),
        color = "Intervention", linetype="Coverage")+
  scale_x_continuous(breaks=c(0,1,2,3,4,5,6,7,8,9,10,11), 
                     labels=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"), 
                     limits = c(0,11))+guides(color=guide_legend(nrow=3,ncol=1,byrow=TRUE,title.position = "top"))+
  
  guides(linetype=guide_legend(nrow=2,ncol=1,byrow=TRUE,title.position = "top"))+
  theme(panel.border = element_blank(), panel.background = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.text=element_text(size=13),title = element_text(size=15, face="bold"),
        axis.line = element_line(colour = "black"),
        axis.text=element_text(size=13,colour="black"),
        axis.text.x=element_text(margin = margin(t = 5)),
        axis.title=element_text(size=13,face="bold"),  
        legend.title=element_text(size=13))+ theme(legend.position = "bottom") +
  theme(legend.background=element_blank()) +ggtitle(paste("high initial clinical incidence"))+   
  theme(legend.key.width = unit(2.5,"cm"))



q <- ggarrange(plotlist = list(
 p1,ggparagraph(text=" ", face = "italic", size = 0.5, color = "black"),
 p2,ggparagraph(text=" ", face = "italic", size =0.5, color = "black")
)  
, nrow=4  , ncol=1,heights  = c(1,0.1,1,0.1),
labels=c("a", "", "b",""),
font.label = list(size = 30, face = "bold"),hjust=-0.1,vjust=1,common.legend = TRUE,legend="bottom")


ggsave(
  filename =  paste0("./Outputs/Exampleplot.jpg"),
  plot = last_plot(),
  
  width = 6,
  height = 8,
  dpi = 600,
  
)

# plot survival estimates 

p1 <- ggplot()+
  geom_line(data=cont_low[[2]], aes(x=trialtime*5/30,y=KM*100,color="cont"),size=1)+
  
  geom_line(data=SMC1_low[[2]], aes(x=trialtime*5/30,y=KM*100,color="SMC1",linetype="high"),size=1)+
  geom_line(data=SMC2_low[[2]], aes(x=trialtime*5/30,y=KM*100,color="SMC1",linetype="low"),size=1)+
  geom_line(data=LAI1_low[[2]], aes(x=trialtime*5/30,y=KM*100,color="LAI1",linetype="high"),size=1)+
  geom_line(data=LAI2_low[[2]], aes(x=trialtime*5/30,y=KM*100,color="LAI1",linetype="low"),size=1)+
  scale_linetype_manual(values = types) +ylim(0,100)+
  scale_color_manual(values = colors,
                     labels=c("Control","LAI",
                              "SMC-SP+AQ"))+
  labs( x ="Month", y = expression(paste("Survival [%]")),
        color = "Intervention", linetype="Coverage")+
  scale_x_continuous(breaks=c(0,1,2,3,4,5,6,7,8,9,10,11), 
                     labels=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"), 
                     limits = c(0,11))+guides(color=guide_legend(nrow=3,ncol=1,byrow=TRUE,title.position = "top"))+
  
  guides(linetype=guide_legend(nrow=2,ncol=1,byrow=TRUE,title.position = "top"))+
  theme(panel.border = element_blank(), panel.background = element_blank(),title = element_text(size=15, face="bold"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.text=element_text(size=13),
        axis.line = element_line(colour = "black"),
        axis.text=element_text(size=13,colour="black"),
        axis.text.x=element_text(margin = margin(t = 5)),
        axis.title=element_text(size=13,face="bold"),  
        legend.title=element_text(size=13))+ theme(legend.position = "bottom") +
  theme(legend.background=element_blank()) +ggtitle(paste("low initial clinical incidence"))+
  theme(legend.key.width = unit(2.5,"cm"))




p2 <- ggplot()+
  geom_line(data=cont_high[[2]], aes(x=trialtime*5/30,y=KM*100,color="cont"),size=1)+
  
  geom_line(data=SMC1_high[[2]], aes(x=trialtime*5/30,y=KM*100,color="SMC1",linetype="high"),size=1)+
  geom_line(data=SMC2_high[[2]], aes(x=trialtime*5/30,y=KM*100,color="SMC1",linetype="low"),size=1)+
  geom_line(data=LAI1_high[[2]], aes(x=trialtime*5/30,y=KM*100,color="LAI1",linetype="high"),size=1)+
  geom_line(data=LAI2_high[[2]], aes(x=trialtime*5/30,y=KM*100,color="LAI1",linetype="low"),size=1)+
  scale_linetype_manual(values = types) +ylim(0,100)+
  scale_color_manual(values = colors,
                     labels=c("Control","LAI",
                              "SMC-SP+AQ"))+
  labs( x ="Month", y = expression(paste("Survival [%]")),
        color = "Intervention", linetype="Coverage")+
  scale_x_continuous(breaks=c(0,1,2,3,4,5,6,7,8,9,10,11), 
                     labels=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"), 
                     limits = c(0,11))+guides(color=guide_legend(nrow=3,ncol=1,byrow=TRUE,title.position = "top"))+
  
  guides(linetype=guide_legend(nrow=2,ncol=1,byrow=TRUE,title.position = "top"))+
  theme(panel.border = element_blank(), panel.background = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.text=element_text(size=13),title = element_text(size=15, face="bold"),
        axis.line = element_line(colour = "black"),
        axis.text=element_text(size=13,colour="black"),
        axis.text.x=element_text(margin = margin(t = 5)),
        axis.title=element_text(size=13,face="bold"),  
        legend.title=element_text(size=13))+ theme(legend.position = "bottom") +
  theme(legend.background=element_blank()) +ggtitle(paste("high initial clinical incidence"))+   
  theme(legend.key.width = unit(2.5,"cm"))



q <- ggarrange(plotlist = list(
  p1,ggparagraph(text=" ", face = "italic", size = 0.5, color = "black"),
  p2,ggparagraph(text=" ", face = "italic", size =0.5, color = "black")
)  
, nrow=4  , ncol=1,heights  = c(1,0.1,1,0.1),
labels=c("a", "", "b",""),
font.label = list(size = 30, face = "bold"),hjust=-0.1,vjust=1,common.legend = TRUE,legend="bottom")


ggsave(
  filename =  paste0("./Outputs/survivalExample.jpg"),
  plot = last_plot(),
  
  width = 6,
  height = 8,
  dpi = 600,
  
)

qsurv <- ggarrange(plotlist = list(
  p1,ggparagraph(text=" ", face = "italic", size = 0.5, color = "black"),
  p2,ggparagraph(text=" ", face = "italic", size =0.5, color = "black")
)  
, nrow=1  , ncol=4,heights  = c(1,1,1,1),widths=c(1,0.1,1,0.1),
labels=c("c", "", "d",""),
font.label = list(size = 30, face = "bold"),hjust=-0.1,vjust=1,common.legend = TRUE,legend="bottom")



## plot cases
p1 <- ggplot()+
  
  geom_line(data=cont_low[[1]], aes(x=trialtime*5/30,y=sum,color="cont"),size=1)+
  
  geom_line(data=SMC1_low[[1]], aes(x=trialtime*5/30,y=sum,color="SMC1",linetype="high"),size=1)+
  geom_line(data=SMC2_low[[1]], aes(x=trialtime*5/30,y=sum,color="SMC1",linetype="low"),size=1)+
  geom_line(data=LAI1_low[[1]], aes(x=trialtime*5/30,y=sum,color="LAI1",linetype="high"),size=1)+
  geom_line(data=LAI2_low[[1]], aes(x=trialtime*5/30,y=sum,color="LAI1",linetype="low"),size=1)+
  scale_linetype_manual(values = types) +
  scale_color_manual(values = colors,
                     labels=c("Control","LAI",
                              "SMC-SP+AQ"))+
  labs( x ="Month", y = expression(paste("cases"["0.25-5y"])),
        color = "Intervention", linetype="Coverage")+ylim(0,125)+
  scale_x_continuous(breaks=c(0,1,2,3,4,5,6,7,8,9,10,11), 
                     labels=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"), 
                     limits = c(0,11))+guides(color=guide_legend(nrow=3,ncol=1,byrow=TRUE,title.position = "top"))+
  
  guides(linetype=guide_legend(nrow=2,ncol=1,byrow=TRUE,title.position = "top"))+
  theme(panel.border = element_blank(), panel.background = element_blank(),title = element_text(size=15, face="bold"),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.text=element_text(size=13),
        axis.line = element_line(colour = "black"),
        axis.text=element_text(size=13,colour="black"),
        axis.text.x=element_text(margin = margin(t = 5)),
        axis.title=element_text(size=13,face="bold"),  
        legend.title=element_text(size=13))+ theme(legend.position = "bottom") +
  theme(legend.background=element_blank()) +ggtitle(paste("low initial clinical incidence"))+
  theme(legend.key.width = unit(2.5,"cm"))




p2 <- ggplot()+
  geom_line(data=cont_high[[1]], aes(x=trialtime*5/30,y=sum,color="cont"),size=1)+
  
  geom_line(data=SMC1_high[[1]], aes(x=trialtime*5/30,y=sum,color="SMC1",linetype="high"),size=1)+
  geom_line(data=SMC2_high[[1]], aes(x=trialtime*5/30,y=sum,color="SMC1",linetype="low"),size=1)+
  geom_line(data=LAI1_high[[1]], aes(x=trialtime*5/30,y=sum,color="LAI1",linetype="high"),size=1)+
  geom_line(data=LAI2_high[[1]], aes(x=trialtime*5/30,y=sum,color="LAI1",linetype="low"),size=1)+

  scale_linetype_manual(values = types) +
  scale_color_manual(values = colors,
                     labels=c("Control","LAI",
                              "SMC-SP+AQ"))+
  labs( x ="Month", y = expression(paste("cases"["0.25-5y"])),
        color = "Intervention", linetype="Coverage")+
  scale_x_continuous(breaks=c(0,1,2,3,4,5,6,7,8,9,10,11), 
                     labels=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"), 
                     limits = c(0,11))+ylim(0,125)+
  guides(color=guide_legend(nrow=3,ncol=1,byrow=TRUE,title.position = "top"))+
  
  guides(linetype=guide_legend(nrow=2,ncol=1,byrow=TRUE,title.position = "top"))+
  theme(panel.border = element_blank(), panel.background = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        legend.text=element_text(size=13),title = element_text(size=15, face="bold"),
        axis.line = element_line(colour = "black"),
        axis.text=element_text(size=13,colour="black"),
        axis.text.x=element_text(margin = margin(t = 5)),
        axis.title=element_text(size=13,face="bold"),  
        legend.title=element_text(size=13))+ theme(legend.position = "bottom") +
  theme(legend.background=element_blank()) +ggtitle(paste("high initial clinical incidence"))+   
  theme(legend.key.width = unit(2.5,"cm"))



qcases <- ggarrange(plotlist = list(
  p1+ theme(legend.position = "none"),ggparagraph(text=" ", face = "italic", size = 0.5, color = "black"),
  p2+ theme(legend.position = "none"),ggparagraph(text=" ", face = "italic", size =0.5, color = "black")
)  
, nrow=1  , ncol=4,heights  = c(1,1,1,1),widths=c(1,0.1,1,0.1),
labels=c("a", "", "b",""),
font.label = list(size = 30, face = "bold"),hjust=-0.1,vjust=1)



qall <- ggarrange(plotlist = list(
  qcases,qsurv)  , heights=c(1,1.3)
, nrow=2  , ncol=1,common.legend = TRUE,legend="bottom")



ggsave(
  filename =  paste0("./Outputs/all_samescale_high.jpg"),
  plot = last_plot(),
  
  width = 12,
  height = 8,
  dpi = 1200,
  
)
