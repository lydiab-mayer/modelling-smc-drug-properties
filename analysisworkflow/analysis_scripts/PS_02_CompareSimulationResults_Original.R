##############################
# script for plotting single OM simulations and comparing prevalence and incidence in multiple age groups
# 
#
# created 22.02.2021
#lydia.burgert@unibas.ch 
#############################

# Setup
rm(list = ls())
library(tgp)
library(dplyr)
library(ggrepel)

set.seed(42)

# User 
user = strsplit(getwd(), "/", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]][5]

# Working directory and group directory
setwd(paste0("/scicore/home/penny/",user,"/M3TPP"))
GROUP_dr = "/scicore/home/penny/GROUP/M3TPP/"

# Source function scripts
source(paste0("./analysisworkflow/analysis_scripts/import_functions.R"))

# insert experiment name here
exp ="..."

# import parameter table for experiment 
param_table_file <- paste0(GROUP_dr,exp,"/param_tab.txt")
param_table <- read.table(param_table_file, sep= "\t", header = TRUE, as.is = TRUE, stringsAsFactors = FALSE)


# choose which categories to plot 

load(paste0(GROUP_dr,exp,"/param_ranges_cat.RData"))

# extract unique experiment combinations 
settings <- unique(param_table[,names(param_cat)])

# set the setting from the settings table 
set_setting <- 1

# filter for settings param table 
params_setting <- merge( param_table,settings[set_setting,], by.x = names(settings))

# define follow up and year before intervention
follow_up=5
years_before_interv=5

# set simulation to plot 
scenario_name <- c("...") #input scenario name from param_table$Scenario_Name column

sim_1 <- import_cont(exp,
                        scenario_name,
                        param_table_file,
                        seeds=max(params_setting$SEED),
                        follow_up,
                        years_before_interv)


scenario_name <- c("...") #input scenario name from param_table$Scenario_Name column

sim_2 <- import_cont(exp,
                     scenario_name,
                     param_table_file,
                     seeds=max(params_setting$SEED),
                     follow_up,
                     years_before_interv)


# plotting incidence   
colors <- c("sim_1" = "#517594", "sim_2"="#42733C")

plot_name <- "cumulative cases per person all age groups "
p1 <- ggplot()+
 geom_line(data=subset(sim_1[["incidence_allages"]],year==10), aes(x=(time-min(time))*5/30,y=cppcum,color="sim_1"),size=1)+
  geom_line(data=subset(sim_2[["incidence_allages"]],year==10), aes(x=(time-min(time))*5/30,y=cppcum,color="sim_2"),size=1)+
  scale_color_manual(values = colors,
                     labels=c("sim_1",
                              "sim_2"))+
  labs( x ="Month", y = expression(paste("cumulative cases per person"["allages"])),
        color = "Simulation")+
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
  theme(legend.background=element_blank()) +ggtitle(plot_name)+
  theme(legend.key.width = unit(2.5,"cm"))
 
ggsave(
  filename= paste0("./Experiments/",exp,"/Outputs/PS_02_",plot_name,".jpg"),
  plot = last_plot(),
  
  width = 9,
  height = 6,
  dpi = 600 )

plot_name <- "cases per person 0-5"
p2 <- ggplot()+
  geom_line(data=subset(sim_1[["incidence_05"]],year==10), aes(x=(time-min(time))*5/30,y=cpp,color="sim_1"),size=1)+
  geom_line(data=subset(sim_2[["incidence_05"]],year==10), aes(x=(time-min(time))*5/30,y=cpp,color="sim_2"),size=1)+
  scale_color_manual(values = colors,
                     labels=c("sim_1",
                              "sim_2"))+
  labs( x ="Month", y = expression(paste("cases per person"["0.25-5y"])),
        color = "Simulation")+
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
  theme(legend.background=element_blank()) +ggtitle(plot_name)+
  theme(legend.key.width = unit(2.5,"cm"))

ggsave(
  filename= paste0("./Experiments/",exp,"/Outputs/PS_02_",plot_name,".jpg"),
  plot = last_plot(),
  
  width = 9,
  height = 6,
  dpi = 600 )

plot_name = "prevalence last year 2-10 years"

p3 <- ggplot()+
  geom_line(data=subset(sim_1[["prevalence_210"]],year==10), aes(x=(time-min(time))*5/30,y=prev*100,color="sim_1"),size=1)+
  geom_line(data=subset(sim_2[["prevalence_210"]],year==10), aes(x=(time-min(time))*5/30,y=prev*100,color="sim_2"),size=1)+
  scale_color_manual(values = colors,
                     labels=c("sim_1",
                              "sim_2"))+
  labs( x ="Month", y = expression(paste(italic("Pf"), "PR"["2-10y"], " [%]" )),
        color = "Simulation")+
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
  theme(legend.background=element_blank()) +ggtitle(plot_name)+
  theme(legend.key.width = unit(2.5,"cm"))

ggsave(
  filename= paste0("./Experiments/",exp,"/Outputs/PS_02_",plot_name,".jpg"),
  plot = last_plot(),
  
  width = 9,
  height = 6,
  dpi = 600 )

plot_name = "prevalence comparison 2-10"

types <- c("after" = "solid", "before"="dashed")
colors <- c("sim_1" = "#517594", "sim_2"="#42733C","cont"="black")

p4 <- ggplot()+
  geom_line(data=subset(sim_1[["prevalence_210"]],year==10), aes(x=(time-min(time))*5/30,y=prev*100,color="sim_1",linetype="after"),size=1)+
  geom_line(data=subset(sim_2[["prevalence_210"]],year==10), aes(x=(time-min(time))*5/30,y=prev*100,color="sim_2",linetype="after"),size=1)+
  geom_line(data=subset(sim_1[["prevalence_210"]],year==5), aes(x=(time-min(time))*5/30,y=prev*100,color="cont",linetype="before"),size=1)+

   scale_color_manual(values = colors,
                     labels=c("sim_1",
                              "sim_2"))+
  labs( x ="Month", y = expression(paste(italic("Pf"), "PR"["2-10y"], " [%]" )),
        color = "Simulation", linetype="Intervention")+
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
  theme(legend.background=element_blank()) +ggtitle(plot_name)+
  theme(legend.key.width = unit(2.5,"cm"))

ggsave(
  filename= paste0("./Experiments/",exp,"/Outputs/PS_02_",plot_name,".jpg"),
  plot = last_plot(),
  
  width = 9,
  height = 6,
  dpi = 600 )
