##############################
# script for plotting continuous parameters over output of interest 
# 
# this script needs to be used after the GPs have been trained, it needs a well trained GP to shows sensible results
# output: predicted mean and confidence intervals over the parameter of interest, depends also on the parameter ranges of the other parameters 
#created 22.03.2021
#lydia.burgert@unibas.ch 
#############################

# Setup
rm(list = ls())
library(tgp)
library(dplyr)
library(ggrepel)
library(tidyverse)
library(hetGP)
library(ggpubr)

set.seed(42)

# User 
user = strsplit(getwd(), "/", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]][5]

# Working directory and group directory
setwd(paste0("/scicore/home/penny/",user,"/M3TPP"))
GROUP_dr = "/scicore/home/penny/GROUP/M3TPP/"

# insert experiment name here
exp ="..."

# insert predicted variable here 
predicted = "..."

# search for gps trained on the desired predicted output
files <- list.files(path = paste0(GROUP_dr,exp,"/gp/trained/",predicted,"/"), pattern = "cv.RData")

# laod parameter ranges over which to plot change in predicted variable 
load(paste0(GROUP_dr,exp,"/param_ranges.RData"))

dataplotlist <- list()
for (k in 1:length(files)){

# choose setting     
setting <- files[k]
scen_name <- str_remove(setting, c(".txt"))

 # load the gaussian process model for the setting of interest  
 gp <- paste0(GROUP_dr,exp,"/gp/trained/",predicted,"/", setting, sep = "")
 load(gp)

 # generate a grid for prediction with the gp
 n_gridpoints=50
 grid <- as.data.frame(mapply(seq, param_ranges_cont[,1], param_ranges_cont[,2], length.out = n_gridpoints))
 colnames(grid) <- rownames(param_ranges_cont)
 
 grid2 <- expand.grid(grid[, rownames(param_ranges_cont)])

# predidct at every grid point using the gp
 pred <- predict(x = as.matrix(grid2),  object = cv_result$GP_model)
 pred_df   <- data.frame(cbind(grid2,pred_mean=pred$mean,pred_sd=pred$sd2))

 # set colours for the single plots   
colors = c("#8fccb4", "#ffcda1","#e89bb9")
plotlist <- list()
  for (j in 1:length(rownames(param_ranges_cont))){ 
 
   to_plot <- rownames(param_ranges_cont)[j]
   
   # aggregate the results for plotting 
   mean <- setNames(aggregate(pred_df[, "pred_mean"], by = list(pred_df[,to_plot ] ) , median), c(to_plot, "mean") )
  CI_high <- setNames(aggregate(pred_df[, "pred_mean"], by = list(pred_df[,to_plot ] ) , FUN=quantile, probs=0.975), c(to_plot,"CIhigh" ))
  CI_low <- setNames(aggregate(pred_df[, "pred_mean"], by = list(pred_df[,to_plot ] ) , FUN=quantile, probs=0.025), c(to_plot,"CIlow" ))
   
  aggdata <- merge(mean, CI_high) %>%
     merge(CI_low)
  
  # plot for predictor of interest 
  plotlist[[j]]<- ggplot(data= aggdata)+
   geom_line(aes_string(to_plot,"mean"),size=1,color=colors[j])+
  geom_ribbon(aes_string(
      x = to_plot,
      ymin = "CIlow",
      ymax = "CIhigh"
     ), fill = colors[j],alpha=0.2)+
     labs( y= predicted,x =to_plot) + ggtitle(to_plot)+
   theme(panel.border = element_blank(), panel.background = element_blank(),
         panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
         legend.text=element_text(size=15),
         axis.line = element_line(colour = "black"),
         axis.text=element_text(size=15,colour="black"),
         axis.text.x=element_text(margin = margin(t = 5)),
         axis.title=element_text(size=15,face="bold"), 
         legend.title=element_text(size=15,face="bold"))+ theme(legend.key=element_blank()) +
   theme(legend.background=element_blank()) + theme(legend.text = element_text(size=15,face="bold"),legend.position = "bottom")+
   theme(legend.key=element_blank(),  legend.text = element_text(size=15,face="bold"))+ 
   guides(colour = guide_legend(title.position = "top",ncol=3,nrow=2,byrow=TRUE))+ 
   theme(axis.title.y=element_text(vjust=1))+theme(legend.key.width = unit(1,"cm"))+ 
   guides(colour = guide_legend(override.aes = list(size=7), title.position = "top",ncol=3,nrow=2,byrow=TRUE))

   } 
  
  

p <- ggarrange( plotlist = plotlist ,nrow = 1,
font.label = list(size = 40, face = "bold"),hjust=0,vjust=1 )

# save plot 
ggsave(
   filename= paste0("./Experiments/",exp,"/Outputs/PS_03_",scen_name,".jpg"),
   plot = last_plot(),
   
   
   dpi = 600 )

               
}

