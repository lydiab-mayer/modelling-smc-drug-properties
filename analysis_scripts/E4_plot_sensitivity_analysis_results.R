

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

setwd("~/smc_lai/analysis_workflow/analysis_scripts")

#S_eff: first order indices
#T_eff: total sensitivity indices (includes interaction)

dir.create(file.path(paste0(mainDir,"/Outputs/E4_sens_new/")), showWarnings = FALSE)

Access <- c(0.1,0.5)
LAI_dec <- c("hill","exp","wei")
seasonaility <- c("Mali","Sen")
IntAge <- c(4.9167)

settings <- expand.grid(Access,LAI_dec,seasonaility,IntAge)
names(settings) <- c("Access","LAI_dec","seasonality","IntAge")

settings <- subset(settings,  Access ==0.1 )

points_low <-c(3,4,8,28,150)
points_high <- c(5,9,20,47,50)

### extract data for range plots 

df_data <- list()
    dfs <- list()
   for (k in 1:nrow(settings)){
     
       setting <- settings[k,]
   if(setting$Access==0.1) {EIRs <- points_low}else{EIRs <- points_high}
   
 
   save_id <- paste0("E4_sens_new_",setting[, "seasonality"],"_",setting[, "LAI_dec"],"_",setting[, "Access"])
  
   dataplotlist_total_lhl <- list()
   dataplotlist_first_lhl <- list()
    dataplotlist_pred_lhl <- list()
  for (i in 1:length(EIRs)) {

    EIR <- EIRs[i]
    if(setting[, "seasonality"]=="Sen"){sens <- "highss"} else {sens <- "highls"}
  # file <- paste('/scicore/home/penny/GROUP/smc_lai/E3E4_lb/sensitivity/gp_trained_pppy_LAI_',setting[, "seasonality"],"_",
  #             setting[ , "IntAge"]  ,'_',setting[, "LAI_dec"],"_",setting[, "Access"],"_",EIR,"_sidx_",sens,".RData", sep = "")
  #       load(file)
  #    
  # 
  #  dataplotlist_total_lhl[[i]] <- c(sobol_idx_list$T_eff[-1],EIR,sobol_idx_list$Access,sobol_idx_list$seasonality,sobol_idx_list$LAI_dec)
  #  dataplotlist_first_lhl[[i]] <- c(sobol_idx_list$S_eff[-1],EIR,sobol_idx_list$Access,sobol_idx_list$seasonality,sobol_idx_list$LAI_dec)
  #  
   
   scen_name <- paste0( setting[, "seasonality"],"_","4.9167","_",setting[, "LAI_dec"],"_",setting[, "Access"],"_",EIR)
 
      filename_pppy_LAI <- paste('/scicore/home/penny/GROUP/smc_lai/E3E4disc_comparison/gp_trained/gp_trained_pppy_LAI_', scen_name, '.RData', sep = "")
      filename_pppy_LAI <- paste('/scicore/home/penny/GROUP/smc_lai/E3E4disc_comparison/gp_trained/gp_trained_pppy_LAI_', scen_name, '.RData', sep = "")
      
       if(file.exists(filename_pppy_LAI)){
       load(filename_pppy_LAI)
    }

    
    Efficacy <- seq(0.7,1,0.05)
    Halflife <- seq(30,150,3)
    Coverage <- seq(0.4,1,0.05)
    Xcand <- expand.grid(Halflife,Efficacy,Coverage)
    names(Xcand) <- c("Halflife","Efficacy","Coverage")
    Xcand$Coverage_SMC <- 0.01
   
    param_col <- c("Coverage_SMC", "Efficacy","Halflife","Coverage")
    
    
    prediction_for_plot <- predict(x = as.matrix(Xcand[,param_col]), 
                                   object = GP_trained_LAI)
    pppy   <- data.frame(cbind(Xcand,pred=prediction_for_plot$mean,EIR,setting$Access,setting$seasonality,setting$LAI_dec ))
    pppy$EIR <- EIR
    dataplotlist_pred_lhl[[i]] <-pppy 
  }
 
    
   dataplotlist_total_shl <- list()
   dataplotlist_first_shl <- list()
    dataplotlist_pred_shl <- list()
  for (i in 1:length(EIRs)) {
   EIR <- EIRs[i]
    if(setting[, "seasonality"]=="Sen"){sens <- "lowss"} else {sens <- "lowls"}
   #  file <- paste('/scicore/home/penny/GROUP/smc_lai/E3E4_lb/sensitivity/gp_trained_pppy_LAI_',setting[, "seasonality"],"_",
   #                setting[ , "IntAge"]  ,'_',setting[, "LAI_dec"],"_",setting[, "Access"],"_",EIR,"_sidx_",sens,".RData", sep = "")
   #      load(file)
   #   
   # 
   # dataplotlist_total_shl[[i]] <- c(sobol_idx_list$T_eff[-1],EIR,sobol_idx_list$Access,sobol_idx_list$seasonality,sobol_idx_list$LAI_dec)
   # dataplotlist_first_shl[[i]] <- c(sobol_idx_list$S_eff[-1],EIR,sobol_idx_list$Access,sobol_idx_list$seasonality,sobol_idx_list$LAI_dec)
   # 
   # 
   scen_name <- paste0( setting[, "seasonality"],"_","4.9167","_",setting[, "LAI_dec"],"_",setting[, "Access"],"_",EIR)
    filename_pppy_LAI <- paste('/scicore/home/penny/GROUP/smc_lai/E3E4disc_comparison/gp_trained/gp_trained_pppy_LAI_', scen_name, '.RData', sep = "")
    if(file.exists(filename_pppy_LAI)){
       load(filename_pppy_LAI)
    }

    
    Efficacy <- seq(0.7,1,0.05)
    Halflife <- seq(30,150,3)
    Coverage <- seq(0.4,1,0.05)
    Xcand <- expand.grid(Halflife,Efficacy,Coverage)
    names(Xcand) <- c("Halflife","Efficacy","Coverage")
    Xcand$Coverage_SMC <- 0.01
   
    param_col <- c("Coverage_SMC", "Efficacy","Halflife","Coverage")
    
    
    prediction_for_plot <- predict(x = as.matrix(Xcand[,param_col]), 
                                   object = GP_trained_LAI)
    pppy   <- data.frame(cbind(Xcand,pred=prediction_for_plot$mean,EIR,setting$Access,setting$seasonality,setting$LAI_dec ))
    pppy$EIR <- EIR
    dataplotlist_pred_shl[[i]] <-pppy 
  }
   
 # total_effects_shl <- as.data.frame(do.call(rbind, dataplotlist_total_shl),stringsAsFactors =FALSE)
 # names(total_effects_shl) <- c(  "Efficacy" ,"Halflife","Coverage","EIR","Access","seasonality","LAI_dec")
 # total_effects_shl$Type <- "shl"
 # 
 # total_effects_lhl <- as.data.frame(do.call(rbind, dataplotlist_total_lhl),stringsAsFactors =FALSE)
 # names(total_effects_lhl) <- c(  "Efficacy" ,"Halflife","Coverage","EIR","Access","seasonality","LAI_dec")
 # total_effects_lhl$Type <- "lhl"
 # 
 #  dataplot <- data.frame(rbind(total_effects_shl,total_effects_lhl))
 # 
 # dfs[[k]]<- dataplot
 
 df_pred1 <-  as.data.frame(do.call(rbind, dataplotlist_pred_shl),stringsAsFactors =FALSE)
 df_pred1$Type <- "shl"
 df_pred2 <-  as.data.frame(do.call(rbind, dataplotlist_pred_lhl),stringsAsFactors =FALSE)
 df_pred2$Type <- "lhl"
 
 df_pred <- as.data.frame(rbind(df_pred1,df_pred2),stringsAsFactors =FALSE)
   names(df_pred) <-c(  "Halflife" ,"Efficacy","Coverage","Coverage_SMC","pred","EIR","Access","seasonality","LAI_dec","Type")

 df_data[[k]] <- df_pred
 
   }
  
  dataplot <- as.data.frame(do.call(rbind, df_data))
  
  dataplot$Coverage <- as.numeric(as.character(dataplot$Coverage))
  dataplot$Halflife <- as.numeric(as.character(dataplot$Halflife))
  dataplot$Efficacy <- as.numeric(as.character(dataplot$Efficacy))
  


 
 
 df <- dataplot
 colors <- c(  "Coverage"="#579c97","Halflife"="#cb6ca2" , "Efficacy"="#f7921e")
 
 

df$regimen = ifelse(df[, "seasonality"]=="Mali","long season","short season")
df$decay = ifelse(df[, "LAI_dec"]=="exp","exponential LAIs",
               ifelse( df[, "LAI_dec"]=="hill","sigmoidal LAIs","bi-phasic LAIs"))

### extract data for side plots 

subsettings <- subset(settings, seasonality=="Mali" & LAI_dec=="hill" & Access==0.1)
df_data <- list()
dfs <- list()
for (k in 1:nrow(subsettings)){
   
   setting <- subsettings[k,]
   if(setting$Access==0.1) {EIRs <- points_low}else{EIRs <- points_high}
   
   
   save_id <- paste0("E4_sens_",setting[, "seasonality"],"_",setting[, "LAI_dec"],"_",setting[, "Access"])
   
   dataplotlist_total_lhl <- list()
   dataplotlist_first_lhl <- list()
   dataplotlist_pred_lhl <- list()
   for (i in 1:length(EIRs)) {
      
      EIR <- EIRs[i]
      if(setting[, "seasonality"]=="Sen"){sens <- "highss"} else {sens <- "highls"}
       file <- paste('/scicore/home/penny/GROUP/smc_lai/E3E4_lb/sensitivity/gp_trained_pppy_LAI_',setting[, "seasonality"],"_",
                   setting[ , "IntAge"]  ,'_',setting[, "LAI_dec"],"_",setting[, "Access"],"_",EIR,"_sidx_",sens,".RData", sep = "")
             load(file)
          
       
        dataplotlist_total_lhl[[i]] <- c(sobol_idx_list$T_eff[-1],EIR,sobol_idx_list$Access,sobol_idx_list$seasonality,sobol_idx_list$LAI_dec)
        dataplotlist_first_lhl[[i]] <- c(sobol_idx_list$S_eff[-1],EIR,sobol_idx_list$Access,sobol_idx_list$seasonality,sobol_idx_list$LAI_dec)
        
      
       }
   
   
   dataplotlist_total_shl <- list()
   dataplotlist_first_shl <- list()
   dataplotlist_pred_shl <- list()
   for (i in 1:length(EIRs)) {
      EIR <- EIRs[i]
      if(setting[, "seasonality"]=="Sen"){sens <- "lowss"} else {sens <- "lowls"}
        file <- paste('/scicore/home/penny/GROUP/smc_lai/E3E4_lb/sensitivity/gp_trained_pppy_LAI_',setting[, "seasonality"],"_",
                      setting[ , "IntAge"]  ,'_',setting[, "LAI_dec"],"_",setting[, "Access"],"_",EIR,"_sidx_",sens,".RData", sep = "")
            load(file)
         
       
       dataplotlist_total_shl[[i]] <- c(sobol_idx_list$T_eff[-1],EIR,sobol_idx_list$Access,sobol_idx_list$seasonality,sobol_idx_list$LAI_dec)
       dataplotlist_first_shl[[i]] <- c(sobol_idx_list$S_eff[-1],EIR,sobol_idx_list$Access,sobol_idx_list$seasonality,sobol_idx_list$LAI_dec)
         }
   
    total_effects_shl <- as.data.frame(do.call(rbind, dataplotlist_total_shl),stringsAsFactors =FALSE)
    names(total_effects_shl) <- c(  "Efficacy" ,"Halflife","Coverage","EIR","Access","seasonality","LAI_dec")
    total_effects_shl$Type <- "shl"
    
    total_effects_lhl <- as.data.frame(do.call(rbind, dataplotlist_total_lhl),stringsAsFactors =FALSE)
    names(total_effects_lhl) <- c(  "Efficacy" ,"Halflife","Coverage","EIR","Access","seasonality","LAI_dec")
    total_effects_lhl$Type <- "lhl"
    
     dataplot <- data.frame(rbind(total_effects_shl,total_effects_lhl))
    
    dfs[[k]]<- dataplot
}

dataplot <- as.data.frame(do.call(rbind, dfs))


dataplot$Coverage <- as.numeric(as.character(dataplot$Coverage))
dataplot$Halflife <- as.numeric(as.character(dataplot$Halflife))
dataplot$Efficacy <- as.numeric(as.character(dataplot$Efficacy))

dataplot$Sum <- rowSums(dataplot[, c(1,2,3)])

dataplot$Coverage <- dataplot$Coverage/dataplot$Sum 
dataplot$Halflife <- dataplot$Halflife/dataplot$Sum 
dataplot$Efficacy <- dataplot$Efficacy/dataplot$Sum 

df2 <-   reshape(dataplot,
                direction = "long",
                varying = list(names(dataplot)[1:3]),
                times =  c(names(dataplot)[1:3]),
                
                v.names = "Value",
                idvar = c("EIR","Access","seasonality","LAI_dec","Coverage","Type"),
                timevar = "Sobol")

colors <- c(  "Coverage"="#579c97","Halflife"="#cb6ca2" , "Efficacy"="#f7921e")



df2$regimen = ifelse(df2[, "seasonality"]=="Mali","long season","short season")
df2$decay = ifelse(df2[, "LAI_dec"]=="exp","exponential LAIs",
                  ifelse( df2[, "LAI_dec"]=="hill","sigmoidal LAIs","bi-phasic LAIs"))


plots <- c("shl","lhl")
plotlist <- list()
pblist <- list()
for (j in 1:2) {
   df_plot <- subset(df2, df2$Type==plots[j])
   if(unique(df_plot$Access)==0.5){ 
      labels <-     c(0.4,1,1.6,2.2,2.9)
      title <- "high access to healthcare"
   }else{
      labels <-  c(0.45,  0.72,1.3,  2.3, 3.2) 
      title <- "low access to healthcare"
    }
   
  
  pp <-  ggplot(df_plot[which(df_plot$seasonality=="Mali" & df_plot$LAI_dec =="hill" &df_plot$Access ==0.1 ), ], aes_string(fill="Sobol", y="Value", x="EIR") )  +
      geom_bar( stat="identity", width = 0.8, position= position_fill(reverse = TRUE))+
      scale_fill_manual(values = colors,
                        labels= c("Deployment coverage","Initial protective efficacy","Protective efficacy half-life"))+ labs(fill=" ") +
      scale_x_discrete( limits = c(as.character(unique(df_plot$EIR))),
                        labels= labels ) +
      labs( x=  expression(paste("Initial cases per person per year"["0.25-5y"])),
            y ="Relative importance")+
      theme_bw()+scale_y_continuous(breaks = seq(0,1,0.2), limits = c(0,1)) + labs(color="Dose [mg]") +
      theme(panel.border = element_blank(),
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            axis.line = element_line(colour = "black"),
            axis.text.x =element_text(size=20,colour="black",vjust=0.7),
            axis.text.y =element_text(size=20,colour="black"),
            axis.title=element_text(size=20,colour="black", face="bold"))+
      theme(strip.background = element_rect(colour="white", fill="white",
                                            size=1.5, linetype="solid"))+
      theme(strip.text = element_text(size=20, color="black", face="bold.italic"))+ theme(legend.position="bottom")+
      theme(legend.text = element_text( size=20),
            legend.title= element_text(size=15, face="bold"))+
      theme(plot.title = element_text(size=17,face="bold"))+ guides(fill = guide_legend(reverse = TRUE,nrow=3),override.aes = list(size=7)) 

   ggsave(
      filename =  paste0("./Outputs/E4_sens_new/",title,".jpg"),
      plot = last_plot(),
      width = 6,
      height = 3,
      dpi = 1200
   )  
   
   
   
   
   
   
   pblist[[j]] <- pp
   }

# generate other plots

data_var <- df # as.data.frame(do.call(rbind, df_data))
data_var$Coverage <- as.numeric(as.character(data_var$Coverage))
data_var$Halflife <- as.numeric(as.character(data_var$Halflife))
data_var$Efficacy <- as.numeric(as.character(data_var$Efficacy))
data_var$Coverage_SMC <- as.numeric(as.character(data_var$Efficacy))

data_var$pred <- as.numeric(as.character(data_var$pred))
data_var$Access <- as.numeric(as.character(data_var$Access))


data_var_plot <- subset(data_var, Access==0.1 & LAI_dec=="exp" & seasonality=="Mali")

data_var_plot2 <- subset(data_var, Access==0.1 & LAI_dec=="hill" & seasonality=="Mali")

normalize <- function(x){
   return((x-min(x)) / (max(x)-min(x)))
}

genplotlist <- function(x){
data_var_plot <- x
   #data %>%
   #group_by(EIR) %>%
   #mutate(pred = normalize(pred))

params <- c("Coverage","Efficacy","Halflife")
params_names <- c("Deployment coverage","Initial protective efficacy","Protective efficacy half-life")


aggdata_Cov_mean <- aggregate(data_var_plot[,c(1,2,3,5)], by = list(data_var_plot$EIR,data_var_plot$Coverage ) , median)
aggdata_Cov_mean$CIhigh <- aggregate(data_var_plot[,c(1,2,3,5)], by = list(data_var_plot$EIR,data_var_plot$Coverage ) , FUN=quantile, probs=0.975)[,"pred"]
aggdata_Cov_mean$CIlow <- aggregate(data_var_plot[,c(1,2,3,5)], by = list(data_var_plot$EIR,data_var_plot$Coverage ) , FUN=quantile, probs=0.025)[,"pred"]

aggdata_Eff_mean <- aggregate(data_var_plot[,c(1,2,3,5)], by = list(data_var_plot$EIR,data_var_plot$Efficacy ) , median)
aggdata_Eff_mean$CIhigh <- aggregate(data_var_plot[,c(1,2,3,5)], by = list(data_var_plot$EIR,data_var_plot$Efficacy ) , FUN=quantile, probs=0.975)[,"pred"]
aggdata_Eff_mean$CIlow <- aggregate(data_var_plot[,c(1,2,3,5)], by = list(data_var_plot$EIR,data_var_plot$Efficacy ) , FUN=quantile, probs=0.025)[,"pred"]

aggdata_HL_mean <- aggregate(data_var_plot[,c(1,2,3,5)], by = list(data_var_plot$EIR,data_var_plot$Halflife ) , median)
aggdata_HL_mean$CIhigh <- aggregate(data_var_plot[,c(1,2,3,5)], by = list(data_var_plot$EIR,data_var_plot$Halflife ) , FUN=quantile, probs=0.975)[,"pred"]
aggdata_HL_mean$CIlow <- aggregate(data_var_plot[,c(1,2,3,5)], by = list(data_var_plot$EIR,data_var_plot$Halflife ) , FUN=quantile, probs=0.025)[,"pred"]

labels <-  c(0.45,  0.72,1.3,  2.3, 3.2) 

colors <- c(  "3"="#d1fdb4", "4"="#8fccb4" ,"8"="#579c97","28" ="#2a6d7a","150"= "#0e3f5c" )
             

p1 <- ggplot(data= aggdata_Cov_mean)+
   geom_line(aes(Coverage*100,pred, color = as.character(Group.1)),size=1)+
  geom_ribbon(aes(
      x = Coverage*100,
      ymin = CIlow,
      ymax = CIhigh,
      fill = as.character(Group.1)),alpha=0.2 )+labs( y=  expression(paste("cases per \nperson per year"["0.25-5y"])),
                                                      x ="Deployment coverage [%]",
                                                      colour =  expression(paste("Initial cases per person per year"["0.25-5y"]))) +
   scale_fill_manual(values=colors,
                     labels=labels,
                     breaks=c("3","4","8","28","150"),drop = FALSE,guide=FALSE) +
   scale_color_manual(values=colors,
                      labels=labels,
                      breaks=c("3","4","8","28","150"),drop = FALSE) +
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

colors <- c(  "3"="#fef5eb", "4"="#ffcda1" ,"8"="#f7921e","28" ="#d35001","150"= "#802a07" )

p2 <- ggplot(data= aggdata_Eff_mean)+
   geom_line(aes(Efficacy*100,pred, color = as.character(Group.1)),size=1)+
   geom_ribbon(aes(
      x = Efficacy*100,
      ymin = CIlow,
      ymax = CIhigh,
      fill = as.character(Group.1)),alpha=0.2 )+labs( y=  expression(paste("cases per \nperson per year"["0.25-5y"])),
                                                      x ="Initial protective efficacy [%]",
                                                      colour =  expression(paste("Initial cases per person per year"["0.25-5y"]))) +
   scale_fill_manual(values=colors,
                     labels=labels,
                     breaks=c("3","4","8","28","150"),drop = FALSE,guide=FALSE) +
   scale_color_manual(values=colors,
                      labels=labels,
                      breaks=c("3","4","8","28","150"),drop = FALSE) +
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

colors <- c(  "3"="#f3cad2", "4"="#e89bb9" ,"8"="#cb6ca2","28" ="#9e4387","150"= "#6d1c68" )


p3 <- ggplot(data= aggdata_HL_mean)+
   geom_line(aes(Halflife,pred, color = as.character(Group.1)),size=1)+
   geom_ribbon(aes(
      x = Halflife,
      ymin = CIlow,
      ymax = CIhigh,
      fill = as.character(Group.1)),alpha=0.2 )+labs( y=  expression(paste("cases per \nperson per year"["0.25-5y"])),
                                                      x ="Protective efficacy half-life [days]",
                                                      colour =  expression(paste("Initial cases per person per year"["0.25-5y   "]))) +
   scale_fill_manual(values=colors,
                     labels=labels,
                     breaks=c("3","4","8","28","150"),drop = FALSE,guide=FALSE) +
   scale_color_manual(values=colors,
                      labels=labels,
                      breaks=c("3","4","8","28","150"),drop = FALSE) +
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

plots <- list(p3,p2,p1)

return(plots)
}

data_var_plot <- subset(data_var, Access==0.1 & LAI_dec=="exp" & seasonality=="Mali")
data_var_plot2 <- subset(data_var, Access==0.1 & LAI_dec=="hill" & seasonality=="Mali")
data_var_plot3 <- subset(data_var, Access==0.1 & LAI_dec=="wei" & seasonality=="Mali")

list_exp_ls <- genplotlist(data_var_plot)
list_hill_ls <- genplotlist(data_var_plot2)
list_wei_ls <- genplotlist(data_var_plot3)

data_var_plot <- subset(data_var, Access==0.1 & LAI_dec=="exp" & seasonality=="Sen")
data_var_plot2 <- subset(data_var, Access==0.1 & LAI_dec=="hill" & seasonality=="Sen")
data_var_plot3 <- subset(data_var, Access==0.1 & LAI_dec=="wei" & seasonality=="Sen")

list_exp_ss <- genplotlist(data_var_plot)
list_hill_ss <- genplotlist(data_var_plot2)
list_wei_ss <- genplotlist(data_var_plot3)




# suppplement plot 
p01 <- ggarrange( plotlist = list(
   ggparagraph(text=" ", face = "italic", size = 0.5, color = "black") ,
   list_exp_ls[[1]],
                                 ggparagraph(text=" ", face = "italic", size = 0.5, color = "black") ,
                                 list_exp_ls[[2]],
                                 ggparagraph(text=" ", face = "italic", size = 0.5, color = "black"),
                                 list_exp_ls[[3]],
                                 ggparagraph(text=" ", face = "italic", size = 0.5, color = "black")
) , ncol = 7, nrow = 1, labels = c("","a","","b","","c"),
font.label = list(size = 40, face = "bold"),hjust=0,vjust=1,heights = c(2,2,2,2,2,2,2),widths=c(0.1,0.9,0.1,0.9,0.1,0.9,0.1) )
p01 <- annotate_figure(p01,
                      top = text_grob(  "long season; exponential LAIs" ,
                                        color = "black", face = "bold", size = 35))
p02 <- ggarrange( plotlist = list(
   ggparagraph(text=" ", face = "italic", size = 0.5, color = "black") ,
   list_wei_ls[[1]],
                                 ggparagraph(text=" ", face = "italic", size = 0.5, color = "black") ,
                                 list_wei_ls[[2]],
                                 ggparagraph(text=" ", face = "italic", size = 0.5, color = "black"),
                                 list_wei_ls[[3]],
                                 ggparagraph(text=" ", face = "italic", size = 0.5, color = "black")
) , ncol = 7, nrow = 1, labels = c("","d","","e","","f"),
font.label = list(size = 40, face = "bold"),hjust=0,vjust=1,heights = c(2,2,2,2,2,2,2),widths=c(0.1,0.9,0.1,0.9,0.1,0.9,0.1) )
p02 <- annotate_figure(p02,
                      top = text_grob(  "long season; bi-phasic LAIs" ,
                                        color = "black", face = "bold", size = 35))







p1 <- ggarrange( plotlist = list(   ggparagraph(text=" ", face = "italic", size = 0.5, color = "black") ,
   list_hill_ls[[1]],
                                  ggparagraph(text=" ", face = "italic", size = 0.5, color = "black") ,
                                 list_hill_ls[[2]],
                                  ggparagraph(text=" ", face = "italic", size = 0.5, color = "black"),
                                 list_hill_ls[[3]],
                                  ggparagraph(text=" ", face = "italic", size = 0.5, color = "black")
) , ncol = 7, nrow = 1, labels = c("","g","","h","","i"),
font.label = list(size = 40, face = "bold"),hjust=0,vjust=1,heights = c(2,2,2,2,2,2,2),widths=c(0.1,0.9,0.1,0.9,0.1,0.9,0.1) )
p1 <- annotate_figure(p1,
                        top = text_grob(  "long season; sigmoidal LAIs" ,
                                          color = "black", face = "bold", size = 35))

p2 <- ggarrange( plotlist = list(   ggparagraph(text=" ", face = "italic", size = 0.5, color = "black") ,
list_exp_ss[[1]],
                                    ggparagraph(text=" ", face = "italic", size = 0.5, color = "black") ,
                                 list_exp_ss[[2]],
                                    ggparagraph(text=" ", face = "italic", size = 0.5, color = "black"),
                                 list_exp_ss[[3]],
                                    ggparagraph(text=" ", face = "italic", size = 0.5, color = "black")
) , ncol = 7, nrow = 1, labels = c("","a","","b","","c"),
font.label = list(size = 40, face = "bold"),hjust=0,vjust=1,heights = c(2,2,2,2,2,2,2),widths=c(0.1,0.9,0.1,0.9,0.1,0.9,0.1) )

p2 <- annotate_figure(p2,
                          top = text_grob(  "short season; exponential LAIs" ,
                                            color = "black", face = "bold", size = 35))


p3 <- ggarrange( plotlist = list(
   ggparagraph(text=" ", face = "italic", size = 0.5, color = "black") ,
   list_wei_ss[[1]],
                                 ggparagraph(text=" ", face = "italic", size = 0.5, color = "black") ,
                                 list_wei_ss[[2]],
                                 ggparagraph(text=" ", face = "italic", size = 0.5, color = "black"),
                                 list_wei_ss[[3]],
                                 ggparagraph(text=" ", face = "italic", size = 0.5, color = "black")
) , ncol = 7, nrow = 1, labels = c("","d","","e","","f"),
font.label = list(size = 40, face = "bold"),hjust=0,vjust=1,heights = c(2,2,2,2,2,2,2),widths=c(0.1,0.9,0.1,0.9,0.1,0.9,0.1) )

p3 <- annotate_figure(p3,
                      top = text_grob(  "short season; bi-phasic LAIs" ,
                                        color = "black", face = "bold", size = 35))


p4 <- ggarrange( plotlist = list(ggparagraph(text=" ", face = "italic", size = 0.5, color = "black") ,
                                  list_hill_ss[[1]],
                                 ggparagraph(text=" ", face = "italic", size = 0.5, color = "black") ,
                                 list_hill_ss[[2]],
                                 ggparagraph(text=" ", face = "italic", size = 0.5, color = "black"),
                                 list_hill_ss[[3]],
                                 ggparagraph(text=" ", face = "italic", size = 0.5, color = "black")
) , ncol = 7, nrow = 1, labels = c("","g","","h","","i"),
font.label = list(size = 40, face = "bold"),hjust=0,vjust=1,heights = c(2,2,2,2,2,2,2),widths=c(0.1,0.9,0.1,0.9,0.1,0.9,0.1) )

p4 <- annotate_figure(p4,
                      top = text_grob(  "short season; sigmoidal LAIs" ,
                                        color = "black", face = "bold", size = 35))

pall <- ggarrange(plotlist=list(p01,p02,p1)      ,
                  ncol = 1, nrow = 3 )





ggsave(
   filename =  paste0("./Outputs/E4_sens_new/Sens_papercols_supp1.jpg"),
   plot = last_plot(),
   width = 15,
   height = 17,
   dpi = 1200
)  

pall <- ggarrange(plotlist=list(p2,p3,p4)      ,
                  ncol = 1, nrow = 3 )





ggsave(
   filename =  paste0("./Outputs/E4_sens_new/Sens_papercols_supp2.jpg"),
   plot = last_plot(),
   width = 15,
   height = 17,
   dpi = 1200
)  


# 
pleft <- ggarrange(plotlist= list(pblist[[1]]+ggtitle("Protective efficacy half-life: 30-90 days"),
                                  pblist[[2]]+ggtitle("Protective efficacy half-life: 90-120 days")),
                       labels = c("a","b"),common.legend=TRUE,legend="bottom",ncol=1,nrow=2,
                   font.label = list(size = 30, face = "bold"),hjust=-0.1,vjust=0.5)
pright <- ggarrange( plotlist = list(list_hill_ls[[1]]+geom_vline(xintercept=90,colour="black",linetype="dotted",size=1.5),
                                     list_hill_ls[[2]],
                                     list_hill_ls[[3]]) , ncol = 1, nrow = 3, labels = c("c","d","e"),
                     font.label = list(size = 30, face = "bold"),hjust=0.8,vjust=1,heights = c(2,2,2),widths=c(0.9,0.9,0.9) )

pall <- ggarrange(plotlist=list(
ggparagraph(text=" ", face = "italic", size = 0.5, color = "black") ,
ggparagraph(text=" ", face = "italic", size = 0.5, color = "black"),
ggparagraph(text=" ", face = "italic", size = 0.5, color = "black"),
ggparagraph(text=" ", face = "italic", size = 0.5, color = "black"),
   pleft, ggparagraph(text=" ", face = "italic", size = 0.15, color = "black"),
   pright, ggparagraph(text=" ", face = "italic", size = 0.5, color = "black")
),
ncol = 4, nrow = 2, widths = c(1,0.07,0.7,0.07),heights=c(0.07,1))

pall2 <-   annotate_figure(pall,
                           top = text_grob(  "long season; sigmoidal LAIs" ,

                                                   color = "black", face = "bold", size = 35),

)

ggsave(
   filename =  paste0("./Outputs/E4_sens_new/Sens_papertest.jpg"),
   plot = last_plot(),
   width = 10,
   height = 12,
   dpi = 1200
)



# publication plot 
pup <- ggarrange( plotlist = list(list_exp_ls[[1]],
                                  ggparagraph(text=" ", face = "italic", size = 0.5, color = "black") ,
                                  list_exp_ls[[2]],
                                  ggparagraph(text=" ", face = "italic", size = 0.5, color = "black"),
                                  list_exp_ls[[3]],
                                  ggparagraph(text=" ", face = "italic", size = 0.5, color = "black")
) , ncol = 6, nrow = 1, labels = c("a","","b","","c"),
font.label = list(size = 40, face = "bold"),hjust=0,vjust=1,heights = c(2,2,2,2,2,2),widths=c(0.9,0.1,0.9,0.1,0.9,0.1) )

pup2 <- annotate_figure(pup,
                        top = text_grob(  "long season; exponential LAIs" ,
                                          color = "black", face = "bold", size = 35))
pdown <- ggarrange( plotlist = list(list_hill_ls[[1]],
                                    ggparagraph(text=" ", face = "italic", size = 0.5, color = "black") ,
                                    list_hill_ls[[2]],
                                    ggparagraph(text=" ", face = "italic", size = 0.5, color = "black"),
                                    list_hill_ls[[3]],
                                    ggparagraph(text=" ", face = "italic", size = 0.5, color = "black")
) , ncol = 6, nrow = 1, labels = c("d","","e","","f"),
font.label = list(size = 40, face = "bold"),hjust=0,vjust=1,heights = c(2,2,2,2,2,2),widths=c(0.9,0.1,0.9,0.1,0.9,0.1) )

pdown2 <- annotate_figure(pdown,
                          top = text_grob(  "long season; sigmoidal LAIs" ,
                                            color = "black", face = "bold", size = 35))

pall <- ggarrange(plotlist=list(pup2,pdown2)      ,
                  ncol = 1, nrow = 2, widths = c(1),heights=c(1,1) )


ggsave(
   filename =  paste0("./Outputs/E4_sens_new/Sens_papercols.jpg"),
   plot = last_plot(),
   width = 15,
   height = 9,
   dpi = 1200
)  
