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

dir.create(paste0("./Experiments/",exp,"/Outputs"))

GROUP_dr = "/scicore/home/penny/GROUP/M3TPP/"

# Source function scripts
source(paste0("./analysisworkflow/analysis_scripts/supp/import_functions.R"))

# insert experiment name here
exp ="..."

# import parameter table for experiment 
param_table_file <- paste0(GROUP_dr,exp,"/param_tab.txt")
param_table <- read.table(param_table_file, sep= "\t", header = TRUE, as.is = TRUE, stringsAsFactors = FALSE)

# choose which categories to plot 

load(paste0(GROUP_dr,exp,"/param_ranges_cat.RData"))

# extract unique experiment combinations 
settings <- unique(param_table[,names(param_cat)])

### setting 1 #####

# set the setting from the settings table 
set_setting <- 1
# filter for settings param table 
params_setting <- merge( param_table,settings[set_setting,], by.x = names(settings))[1, ]

# rename the object here to your prefered setting identifier
AGO_EIR5 = import_EIRs_cat(exp,params_setting, seeds=max(params_setting$SEED))

print(paste0("input EIR=",sum(AGO_EIR5$input.EIR)))
print(paste0("simulated EIR=",sum(AGO_EIR5$simulated.EIR)))

### setting 2 #####
# set the setting from the settings table 
set_setting <- 2
# filter for settings param table 
params_setting <- merge( param_table,settings[set_setting,], by.x = names(settings))[1, ]

# rename the object here to your prefered setting identifier
AGO_EIR10 = import_EIRs_cat(exp,params_setting, seeds=max(params_setting$SEED))

print(paste0("input EIR=",sum(AGO_EIR10$input.EIR)))
print(paste0("simulated EIR=",sum(AGO_EIR10$simulated.EIR)))

#### plot settings ####

colors <- c("AGO_EIR10" = "#517594", "AGO_EIR5"="#42733C")
types <- c("simulated" = "solid", "input"="dashed")

plot_name <- "input vs simulated EIR"


p = ggplot()+
 geom_line(data=AGO_EIR10, aes(x=timestep*5/30, y=input.EIR,color="AGO_EIR10",linetype="input")) + 
  geom_line(data=AGO_EIR10, aes(x=timestep*5/30, y=simulated.EIR,color="AGO_EIR10",linetype="simulated")) + 
  geom_line(data=AGO_EIR5, aes(x=timestep*5/30, y=input.EIR,color="AGO_EIR5",linetype="input")) + 
  geom_line(data=AGO_EIR5, aes(x=timestep*5/30, y=simulated.EIR,color="AGO_EIR5",linetype="simulated")) + 
   scale_color_manual(values = colors,
                     labels=names(colors),)+
  labs( x ="Month", y = "EIR",
        color = "Setting")+
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
  filename= paste0("./Experiments/",exp,"/Outputs/PS_00_",plot_name,".jpg"),
  plot = last_plot(),
  
  width = 9,
  height = 6,
  dpi = 600 )
