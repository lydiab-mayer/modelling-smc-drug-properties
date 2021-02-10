setwd("~/smc_lai/analysis_workflow/analysis_scripts")

library(ggplot2)
library(ggpubr)
library(reshape2)

#S_eff: first order indices
#T_eff: total sensitivity indices (includes interaction)
<<<<<<< Updated upstream
df <- data.frame(seq(1.0,25,0.5))
=======
df <- data.frame(EIR=c(seq(1.0,10.5,0.5)))
>>>>>>> Stashed changes
colnames(df) <- "EIR"


df$first_order_EIR <- 0
df$first_order_Coverage <- 0
df$first_order_Halflife <- 0
df$first_order_Efficacy <- 0

df$total_EIR <- 0
df$total_Coverage <- 0
df$total_Halflife <- 0
df$total_Efficacy <- 0

LAI_dec <- "exp"

for (k in 1:nrow(df)){ 
<<<<<<< Updated upstream
 # load(paste('/scicore/home/smith/GROUP/smc_lai/E2_LAI/gp_5/sensitivity/seeds_E2_LAI_seasonal_Low_indoor_cv_EIR_',
#             formatC(df[k,]$EIR, digits = 1, format = "f"),'_sidx.RData',
 #            sep = ""))
  load(paste('/scicore/home/smith/GROUP/smc_lai/E3_E4_comparison/sensitivity_analysis/seeds_E3_E4_comparison_seasonal_Low_indoor_cv_EIR_',
            formatC(df[k,]$EIR, digits = 1, format = "f"),'_sidx.RData',
           sep = ""))
=======
  load(paste("/scicore/home/smith/GROUP/smc_lai/E2_LAI/gp_5/sensitivity/seeds_E2_LAI_seasonal_Low_indoor_",LAI_dec,"_cv_EIR_",
             formatC(df[k,]$EIR, digits = 1, format = "f"),"_sidx.RData",
             sep = ""))
 # load(paste('/scicore/home/smith/GROUP/smc_lai/E3_E4_comparison/sensitivity_analysis/seeds_E3_E4_comparison_seasonal_Low_indoor_cv_EIR_',
 #           formatC(df[k,]$EIR, digits = 1, format = "f"),'_sidx.RData',
 #          sep = ""))
>>>>>>> Stashed changes
  df[k,]$first_order_EIR <- sobol_idx_list$S_eff[1]
  df[k,]$first_order_Coverage <- sobol_idx_list$S_eff[2]
  df[k,]$first_order_Halflife <- sobol_idx_list$S_eff[3]
  df[k,]$first_order_Efficacy <- sobol_idx_list$S_eff[4]
  
  df[k,]$total_EIR <- sobol_idx_list$T_eff[1]
  df[k,]$total_Coverage <- sobol_idx_list$T_eff[2]
  df[k,]$total_Halflife <- sobol_idx_list$T_eff[3]
  df[k,]$total_Efficacy <- sobol_idx_list$T_eff[4]
}


df_plot <- data.frame(rep(df$EIR, 4))
colnames(df_plot) <- c("EIR")
df_plot$variable <- rep(c("EIR", "coverage", "halflife", "efficacy"),each = length(df$EIR))
df_plot$sobol_first_order <- c(df$first_order_EIR, df$first_order_Coverage,
                               df$first_order_Halflife,df$first_order_Efficacy)
df_plot$sobol_total <- c(df$total_EIR, df$total_Coverage,
                         df$total_Halflife,df$total_Efficacy)
df_plot <- melt(df_plot, id = c("EIR","variable"))

colnames(df_plot) <- c("EIR","variable","sensitivity_index", "value")

colors <- c("SMC" = "#66FF33", "LAI"="#66CCFF", "LAI2"="#0033FF")


p <- ggplot(data = df_plot[df_plot$variable != "EIR",]) +
 # geom_point(aes(x = EIR, y = value, colour = variable)) + 
  geom_line(aes(x = EIR, y = value, colour = variable, linetype = sensitivity_index),size=1) + 
  
  labs(#title = "Sensitivity analysis",
     #  subtitle = "relative difference between SMC and LAI" , 
       x = "EIR", y = "Sensitivity", linetype = "sensitivity index") +
  scale_linetype_manual(values = c("sobol_first_order" = 1, "sobol_total" = 2),
                        labels = c("first order", "total")) +
  scale_color_manual(values = c("#00AFBB", "#E7B800", "#FC4E07"))+ 
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
<<<<<<< Updated upstream
  filename =  paste0('./Outputs/E3_E4_sensitivity.jpg'),
=======
  filename =  paste0("./Outputs/E2_LAIsensitivity_",LAI_dec,".jpg"),
>>>>>>> Stashed changes
  fig = last_plot(),
  width = 15,
  height = 10,
  dpi = 1200
)  
