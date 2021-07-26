####################################################################################################################
################################################## PLOT OPTIMIZATION 2

# Setup
rm(list = ls())
library(tgp)
library(dplyr)
library(ggrepel)
library(data.table)
library(akima)
library(ggpubr)
library(ggplot2)

set.seed(42)

#########################################
###
### Mock example plots: optimization
### NN 02/06/2021
###
#########################################


# List of values
# exp1 = c("MOCK_EIR15_sig2_april_age10","MOCK_EIR15_sig2_may_age10","MOCK_EIR15_sig2_june_age10")
# exp2 = c("MOCKEIR15sig2aprilage10","MOCKEIR15sig2mayage10","MOCKEIR15sig2juneage10")
# age = c("10","10","10")
# EIR = c("15","15","15")
exp1 = c("MOCK_EIR15_sig2_april_age10","MOCK_EIR15_sig2_may_age10")
exp2 = c("MOCKEIR15sig2aprilage10","MOCKEIR15sig2mayage10")
age = c("10","10")
EIR = c("15","15")

season = c("seas4mo")
predicted = "inc_red_int_5mo"
reductions_list <- c(40,50,60,70)

# Load data
all_data = data.frame()
for(k in 1:length(exp1)){
  for(i in 1:length(season)){
    for(j in 1:length(reductions_list)){
      
      # File path
      file = paste0("/scicore/home/penny/GROUP/M3TPP/",exp1[[k]],"/gp/optimisation/",predicted,"/seeds_",exp2[[k]],"_",season[[i]],"_Equal_",EIR[[k]],"_",age[[k]],"_exp_0.1_",predicted,"_cv_Halflife_cutoff",reductions_list[[j]],"_opt.txt")
      
      # Check if exists
      if(file.exists(file)){
        cat("Exists for",exp1[[k]],season[[i]],reductions_list[[j]],"\n")
        data = as.data.frame(read.table(file, header = T))
        data = data[,1:3]
        data = cbind(data,exp=exp1[[k]],target=reductions_list[[j]],season=season[[i]])
      }else{
        cat("Does not exist for",season[[i]],reductions_list[[j]],"\n")
        data = NULL
      }
      
      # Combine
      all_data = rbind(all_data,data)
    }
  }
}

# Factor
all_data$target = factor(all_data$target, levels = reductions_list)

# Categorize optimal half-life
all_data = all_data %>% 
  mutate(optimal_Halflife_cat = cut(optimal_Halflife, 
                                    breaks = seq(30,180,30),
                                    labels = c("1 - 2",
                                               "2 - 3",
                                               "3 - 4",
                                               "4 - 5",
                                               "5 - 6")))

# Season labels
all_data$season_lab = NA
all_data$season_lab = ifelse(all_data$season == "seas4mo", "4 months", all_data$season_lab)

# Intervention labels
all_data$intervention = NA
all_data$intervention = ifelse(all_data$exp == "MOCK_EIR15_sig2_april_age10", "April", all_data$intervention)
all_data$intervention = ifelse(all_data$exp == "MOCK_EIR15_sig2_may_age10", "May", all_data$intervention)
# all_data$intervention = ifelse(all_data$exp == "MOCK_EIR15_sig2_june_age10", "June", all_data$intervention)
all_data$intervention = factor(all_data$intervention, levels = c("April","May"))

#########################

# Create an Eff-Cov variable
all_data$EffCov = paste0(round(all_data$Efficacy*100, digits=0),"% Efficacy\n",round(all_data$Coverage*100, digits=0),"% Coverage")
unique(all_data$EffCov)

colors0=viridis(5)[2:5]

# For 80% efficacy and 80% coverage
ggplot(data = subset(all_data, season == "seas4mo" & EffCov %in% c("72% Efficacy\n72% Coverage","83% Efficacy\n83% Coverage","94% Efficacy\n94% Coverage")), 
       aes(x=intervention, y=optimal_Halflife, color = EffCov)) +
  facet_wrap(.~target) +
  geom_point(size = 3.5) +
  geom_text(aes(y=optimal_Halflife+15, label=EffCov),hjust="left", vjust=0.5) +
  scale_color_manual(values=colors0, name = "% incidence\nreduction") +
  scale_y_continuous(breaks = seq(30,180,30), limits = c(30,180)) +
  labs(x = NULL,
       y = "Duration (days)") +
  theme(panel.border = element_rect(color = "black", size = 1, fill = NA),
        panel.background = element_rect(fill = NA),
        panel.grid = element_line(color = "gray90"),
        legend.position = "none",
        title = element_text(size=11, face="bold"),
        legend.text=element_text(size=11),
        axis.line = element_line(colour = "black"),
        axis.text=element_text(size=11,colour="black"),
        axis.title=element_text(size=11,face="bold"),
        strip.background = element_rect(fill = "gray95"),
        strip.text = element_text(size=11,face="bold"))



dev.copy(png, file = paste0("/scicore/home/penny/nekkab0000/M3TPP/Optimization_3strategies_april_may_alltargets.png"),
         width = 2*2400, height = 2*1100, res = 300)
dev.off()


#########################

# Create an Eff-Cov variable
all_data$EffCov = paste0(round(all_data$Efficacy*100, digits=0),"% Efficacy & ",round(all_data$Coverage*100, digits=0),"% Coverage")
unique(all_data$EffCov)

colors0=viridis(5)[2:5]

# For 80% efficacy and 80% coverage
subdata=subset(all_data, season == "seas4mo" & target == 60 & Coverage < 1 & Efficacy < 1 & Coverage > 0.8 & Efficacy > 0.8)
unique(subdata$EffCov)
may_EffCov = subdata$EffCov[which(subdata$intervention == "May" & !is.na(subdata$optimal_Halflife))]
subdata = subdata[which(subdata$EffCov %in% may_EffCov),]

ggplot(data = subdata, 
       aes(x=intervention, y=optimal_Halflife, color = EffCov)) +
  geom_point(size = 3.5) +
  geom_text(aes(y=optimal_Halflife+5, label=EffCov), hjust=-0.1, vjust = 1.4, size = 6) +
  # scale_color_manual(values=colors0, name = "% incidence\nreduction") +
  scale_color_viridis(discrete = T, name = "% incidence\nreduction") +
  scale_y_continuous(breaks = seq(30,180,30), limits = c(30,180)) +
  labs(x = NULL,
       y = "Duration (days)") +
  theme(panel.border = element_rect(color = "black", size = 1, fill = NA),
        panel.background = element_rect(fill = NA),
        panel.grid = element_line(color = "gray90"),
        legend.position = "none",
        title = element_text(size=11, face="bold"),
        legend.text=element_text(size=11),
        axis.line = element_line(colour = "black"),
        axis.text=element_text(size=11,colour="black"),
        axis.title=element_text(size=11,face="bold"),
        strip.background = element_rect(fill = "gray95"),
        strip.text = element_text(size=11,face="bold"))



dev.copy(png, file = paste0("/scicore/home/penny/nekkab0000/M3TPP/Optimization_3strategies_april_may_alltargets.png"),
         width = 4500, height = 1800, res = 300)
dev.off()
