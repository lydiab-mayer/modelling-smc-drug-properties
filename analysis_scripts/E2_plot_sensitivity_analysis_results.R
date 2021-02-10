setwd("~/smc_lai/analysis_workflow/analysis_scripts")

library(ggplot2)
library(ggpubr)
library(reshape2)

#S_eff: first order indices
#T_eff: total sensitivity indices (includes interaction)

dir.create(file.path(paste0(mainDir,"/Outputs/E2_Sens/")), showWarnings = FALSE)

Access <- c(0.1,0.5)
LAI_dec <- c("hill","exp")
seasonaility <- c("Mali","Sen")
IntAge <- c(4.9167)

settings <- expand.grid(Access,LAI_dec,seasonaility,IntAge)
names(settings) <- c("Access","LAI_dec","seasonality","IntAge")


EIRs <- c(1,3,5,7,9,10,12,14,16,18,20,22,24,50,100,120,140,160,180,200)

plot_sens_E2 <- function(setting,EIRs) {
  
  dataplotlist_total <- list()
  dataplotlist_first <- list()
 
   save_id <- paste0("E2_sens_",setting[, "seasonality"],"_",setting[, "LAI_dec"],"_",setting[, "Access"],"_",setting[, "IntAge"])
  
  
  for (i in 1:length(EIRs)) {

    EIR <- EIRs[i]
   load(paste('/scicore/home/smith/GROUP/smc_lai/E2_LAI/gp_5/sensitivity/seeds_E2_LAI_',setting[, "seasonality"],"_",
              setting[, "IntAge"]  ,'_',setting[, "LAI_dec"],"_",setting[, "Access"],"_cv_EIR_",EIR,'_sidx.RData', sep = "") ) 
                      
   dataplotlist_total[[i]] <- c(sobol_idx_list$T_eff[-1],EIR)
   dataplotlist_first[[i]] <- c(sobol_idx_list$S_eff[-1], EIR)
   
  }
   
 total_effects <- ldply(dataplotlist_total)
 names(total_effects) <- c(  "Coverage" ,"Halflife" ,"Efficacy","EIR")
 first_effects <- ldply(dataplotlist_first)
 names(first_effects) <- c(  "Coverage" ,"Halflife" ,"Efficacy","EIR")
 
 colors <- c(  "Coverage"="#00AFBB", "Halflife"="#E7B800" ,"Efficacy"="#FC4E07")
 
 
 q <- ggplot() +
   geom_line(data = total_effects,aes(x = EIR, y = Coverage, color = "Coverage", linetype = "sobol_total"),size=1) +
   geom_line(data = total_effects,aes(x = EIR, y = Halflife, color = "Halflife", linetype = "sobol_total"),size=1) + 
   geom_line(data = total_effects,aes(x = EIR, y = Efficacy, color = "Efficacy", linetype = "sobol_total"),size=1) + 
   geom_line(data = first_effects,aes(x = EIR, y = Coverage, color = "Coverage", linetype = "sobol_first_order"),size=1) +
   geom_line(data = first_effects,aes(x = EIR, y = Halflife, color = "Halflife", linetype = "sobol_first_order"),size=1) + 
   geom_line(data = first_effects,aes(x = EIR, y = Efficacy, color = "Efficacy", linetype = "sobol_first_order"),size=1) + 
   
   
     labs( x = "EIR", y = "Sensitivity") +
   scale_linetype_manual(values = c("sobol_first_order" = 1, "sobol_total" = 2),
                         labels = c("first order", "total")) +
   scale_color_manual(values = colors)+
   theme(plot.title = element_text(face="bold"))+
   theme(panel.border = element_blank(), panel.background = element_blank(),
         panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
         legend.text=element_text(size=15),
         axis.line = element_line(colour = "black"),
         axis.text=element_text(size=15,colour="black"),
         axis.text.x=element_text(margin = margin(t = 5)),
         axis.title=element_text(size=15,face="bold"),  
         legend.title=element_text(size=15))+
   theme(legend.background=element_blank()) + theme(legend.key=element_blank(),
                                                    legend.text = element_text(size=15,face="bold"),
                                                    legend.title = element_text(size=15,face="bold") )
 
 
 save_plot(
   filename =  paste0("./Outputs/E2_Sens/Sens_",save_id,".jpg"),
   fig = last_plot(),
   width = 35,
   height = 15,
   dpi = 1200
 )  
 
}
 
for (k in 1:nrow(settings)) {
  setting <- settings[k, ]
  plot_sens_E2(setting,EIRs)
}
