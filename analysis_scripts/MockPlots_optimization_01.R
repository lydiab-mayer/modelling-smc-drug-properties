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
### Mock example plots: optimisation
### NN 02/06/2021
###
#########################################


# List of values
# exp1 = c("MOCK_1a","MOCK_1b","MOCK_1c")
# exp2 = c("MOCK1a","MOCK1b","MOCK1c")
# exp1 = c("MOCK_1","MOCK_1b")
# exp2 = c("MOCK1","MOCK1b")
# age = c("5","10","90")
# age = c("10","10")
# EIR = c("20","100")
# exp1="MOCK_EIR15_sig_april_age10"
# exp2="MOCKEIR15sigaprilage10"

# exp1="MOCK_EIR15_sig2_may_age10"
# exp2="MOCKEIR15sig2mayage10"

exp1="MOCK_EIR15_sig2_april_age10"
exp2="MOCKEIR15sig2aprilage10"

age="10"
EIR="15"
exp = exp1

season = c("seas4mo")
predicted = "inc_red_int_5mo"
reductions_list <- seq(10,80,10)

# Load data
all_data = data.frame()
for(k in 1:length(exp1)){
  for(i in 1:length(season)){
    for(j in 1:length(reductions_list)){
      
      # File path
      file = paste0("/scicore/home/penny/GROUP/M3TPP/",exp1,"/gp/optimisation/",predicted,"/seeds_",exp2,"_",season,"_Equal_",EIR,"_",age,"_exp_0.1_",predicted,"_cv_Halflife_cutoff",reductions_list[[j]],"_opt.txt")
      
      # Check if exists
      if(file.exists(file)){
        cat("Exists for",exp1,season,reductions_list[[j]],"\n")
        data = as.data.frame(read.table(file, header = T))
        data=data[,1:3]
        data = cbind(data,exp=exp1,target=reductions_list[[j]],season=season)
      }else{
        cat("Does not exist for",season,reductions_list[[j]],"\n")
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

# Colors
colorsX = viridis(length(seq(10,80,10)))
names(colorsX)=seq(10,80,10)


#################### PLOT: show that optimal half-life increases with seasonal length

index_minEff = which.min(abs(all_data$Efficacy - 0.9))
index_minCov = which.min(abs(all_data$Coverage - 0.85))
value_minEff = all_data$Efficacy[index_minEff]
value_minCov = all_data$Coverage[index_minCov]

# Create an Eff-Cov variable
all_data$EffCov = paste0(round(all_data$Efficacy*100, digits=0),"%-",round(all_data$Coverage*100, digits=0),"%")
unique(all_data$EffCov)

colors0=viridis(5)[2:4]

# For 80% efficacy and 80% coverage
ggplot(data = subset(all_data, target == 40 & EffCov %in% c("72%-72%","83%-83%","94%-94%")), 
       aes(x=season, y=optimal_Halflife, color = EffCov, group = EffCov)) +
  scale_color_manual(values=colors0, name = "% incidence\nreduction") +
  geom_point() +
  geom_text(aes(label=EffCov),hjust=-0.25, vjust=0.5) +
  theme(legend.position = "none") +
  scale_x_discrete(labels = c("4","5","6","9")) +
  labs(x = "Seasonal duration (months)",
       y = "Duration (days)") +
  theme(legend.position = "none",
        panel.border = element_rect(color = "black", size = 1, fill = NA),
        panel.background = element_rect(fill = NA),
        panel.grid = element_line(color = "gray90")) +
  theme(title = element_text(size=11, face="bold"),
        legend.text=element_text(size=11),
        axis.line = element_line(colour = "black"),
        axis.text=element_text(size=11,colour="black"),
        axis.title=element_text(size=11,face="bold"),
        strip.background = element_rect(fill = "gray95"),
        strip.text = element_text(size=11,face="bold"),
        panel.background = element_rect(fill = "gray98"))


###############################################################################

index_minEff80 = which.min(abs(all_data$Efficacy - 0.9))
index_minCov80 = which.min(abs(all_data$Coverage - 0.85))
value_minEff80 = all_data$Efficacy[index_minEff80]
value_minCov80 = all_data$Coverage[index_minCov80]

index_minEff95 = which.min(abs(all_data$Efficacy - 0.95))
index_minCov60 = which.min(abs(all_data$Coverage - 0.8))
value_minEff95 = all_data$Efficacy[index_minEff95]
value_minCov60 = all_data$Coverage[index_minCov60]

colors1=viridis(5)[2:3]
colors2=viridis(5)[4:5]

# Choose incidence reduction
target_lvl = "60"

#################

sub_data = subset(all_data, target == target_lvl & Coverage >= value_minCov80)

ggplot(data= sub_data, 
       aes(x=Efficacy, 
           y=optimal_Halflife, 
           color=factor(Coverage),
           group=factor(Coverage))) +
  geom_point(size = 2) +
  geom_line() 

dev.copy(png, file = paste0("/scicore/home/penny/nekkab0000/M3TPP/Experiments/",exp,"/Outputs/optimization_target",target_lvl,".png"),
         width = 2400, height = 1100, res = 300)
dev.off()

######

GROUP_dr = "/scicore/home/penny/GROUP/M3TPP/"
exp ="MOCK_EIR15_sig2_may_age10"
post4mo_agg = read.table(file = paste0(GROUP_dr,exp,"/postprocessing/agg_MOCKEIR15sig2mayage10_seas4mo_Equal_15_10_exp_0.1.txt"), header = T, sep = "\t")
post4mo_agg$inc_red_int_5mo_round = round(post4mo_agg$inc_red_int_5mo, digits = -1)
post4mo_agg$Coverage_round = round(post4mo_agg$Coverage, digits = 2)
post4mo_agg$Efficacy_round = round(post4mo_agg$Efficacy, digits = 2)
post4mo_agg$Halflife_cat = cut(post4mo_agg$Halflife, seq(30,180,15))
post4mo_agg$Coverage_round_cat = cut(post4mo_agg$Coverage_round, seq(0,1,0.05))
sub_post4mo_agg = subset(post4mo_agg, inc_red_int_5mo_round == 60 & Coverage_round >= 0.80)

ggplot(data=sub_post4mo_agg, 
       aes(x=Efficacy_round, 
           y=Halflife_cat, 
           color=Coverage_round_cat,
           group=Coverage_round_cat)) +
  geom_point(size = 5)

dev.copy(png, file = paste0("/scicore/home/penny/nekkab0000/M3TPP/Experiments/",exp,"/Outputs/postprocessing_target",target_lvl,".png"),
         width = 2400, height = 1100, res = 300)
dev.off()

######################

# Efficacy
p1a=ggplot(data= subset(all_data, target == target_lvl & Coverage %in% c(value_minCov60)), 
       aes(x=Efficacy, y=optimal_Halflife, color=factor(round(Coverage*100, digits =-1)), group=factor(round(Coverage*100, digits =-1)))) +
  geom_line(size = 1.8) +
  # geom_ribbon(aes(ymin=sd_minus,ymax=sd_plus, fill = factor(round(Coverage*100, digits =0))), alpha=0.2, color = NA) +
  scale_color_manual(values = colors1[1], name = "Coverage", labels = function(x) paste0(x, "%")) +
  # scale_fill_manual(values = colors1[1], name = "Coverage [95%CI]", labels = function(x) paste0(x, "%")) +
  scale_y_continuous(breaks = seq(30,180,30), limits = c(30,180)) +
  scale_x_continuous(breaks = seq(0.5,1,0.1), labels = function(x) scales::percent(x, accuracy = 1)) +
  labs(x = "Initial protective efficacy", y = "Duration (days)") +
  coord_cartesian(expand = c(0,0)) +
  theme(legend.position = "right",
        panel.border = element_rect(color = "black", size = 1, fill = NA),
        panel.background = element_rect(fill = NA),
        panel.grid = element_line(color = "gray90")) +
  theme(title = element_text(size=11, face="bold"),
        plot.background = element_blank(),
        legend.text=element_text(size=11),
        axis.line = element_line(colour = "black"),
        axis.text=element_text(size=11,colour="black"),
        axis.title=element_text(size=11,face="bold"),
        strip.background = element_rect(fill = "gray95"),
        strip.text = element_text(size=11,face="bold"),
        panel.grid.minor = element_blank(),
        legend.position = "right")
p1a

# Save
# insert experiment name here
dev.copy(png, file = paste0("/scicore/home/penny/nekkab0000/M3TPP/Experiments/",exp,"/Outputs/new2_Optimization_target",target_lvl,"_Eff_Cov80.png"),
         width = 2400, height = 1100, res = 300)
dev.off()


# Efficacy
p1b=ggplot(data= subset(all_data, target == target_lvl & Coverage %in% c(value_minCov80)), 
           aes(x=Efficacy, y=optimal_Halflife, color=factor(round(Coverage*100, digits =-1)), group=factor(round(Coverage*100, digits =-1)))) +
  geom_line(size = 1.8) +
  # geom_ribbon(aes(ymin=sd_minus,ymax=sd_plus, fill = factor(round(Coverage*100, digits =0))), alpha=0.2, color = NA) +
  scale_color_manual(values = colors1[2], name = "Coverage", labels = "85%") +
  # scale_fill_manual(values = colors1[1], name = "Coverage [95%CI]", labels = function(x) paste0(x, "%")) +
  scale_y_continuous(breaks = seq(30,180,30), limits = c(30,180)) +
  scale_x_continuous(breaks = seq(0.5,1,0.1), labels = function(x) scales::percent(x, accuracy = 1)) +
  labs(x = "Initial protective efficacy", y = "Duration (days)") +
  coord_cartesian(expand = c(0,0)) +
  theme(legend.position = "right",
        panel.border = element_rect(color = "black", size = 1, fill = NA),
        panel.background = element_rect(fill = NA),
        panel.grid = element_line(color = "gray90")) +
  theme(title = element_text(size=11, face="bold"),
        plot.background = element_blank(),
        legend.text=element_text(size=11),
        axis.line = element_line(colour = "black"),
        axis.text=element_text(size=11,colour="black"),
        axis.title=element_text(size=11,face="bold"),
        strip.background = element_rect(fill = "gray95"),
        strip.text = element_text(size=11,face="bold"),
        panel.grid.minor = element_blank(),
        legend.position = "right")
p1b

# Save
# insert experiment name here
dev.copy(png, file = paste0("/scicore/home/penny/nekkab0000/M3TPP/Experiments/",exp,"/Outputs/new_Optimization_target",target_lvl,"_Eff_Cov85.png"),
         width = 2400, height = 1100, res = 300)
dev.off()

#############


# Coverage
p2a=ggplot(data= subset(all_data, target == target_lvl & Efficacy %in% c(value_minEff80)), 
          aes(x=Coverage, y=optimal_Halflife, color=factor(round(Efficacy*100, digits =-1)), group=factor(round(Efficacy*100, digits =-1)))) +
  geom_line(size = 2) +
  # geom_ribbon(aes(ymin=sd_minus,ymax=sd_plus, fill = factor(round(Efficacy*100, digits =0))), alpha=0.2, color = NA) +
  scale_color_manual(values = colors1[1], name = "Efficacy", labels = function(x) paste0(x, "%")) +
  # scale_fill_manual(values = colors1, name = "Efficacy [95%CI]", labels = function(x) paste0(x, "%")) +
  scale_y_continuous(breaks = seq(30,180,30), limits = c(30,180)) +
  scale_x_continuous(breaks = seq(0.5,1,0.1), labels = function(x) scales::percent(x, accuracy = 1)) +
  labs(x = "Coverage", y = "Duration (days)") +
  coord_cartesian(expand = c(0,0)) +
  theme(legend.position = "right",
        panel.border = element_rect(color = "black", size = 1, fill = NA),
        panel.background = element_rect(fill = NA),
        panel.grid = element_line(color = "gray90")) +
  theme(title = element_text(size=11, face="bold"),
        plot.background = element_blank(),
        legend.text=element_text(size=11),
        axis.line = element_line(colour = "black"),
        axis.text=element_text(size=11,colour="black"),
        axis.title=element_text(size=11,face="bold"),
        strip.background = element_rect(fill = "gray95"),
        strip.text = element_text(size=11,face="bold"),
        panel.grid.minor = element_blank(),
        legend.position = "right")
p2a

# Save
dev.copy(png, file = paste0("/scicore/home/penny/nekkab0000/M3TPP/Experiments/",exp,"/Outputs/new_Optimization_target",target_lvl,"_Cov_Eff90.png"),
         width = 2400, height = 1100, res = 300)
dev.off()


# Coverage
p2b=ggplot(data= subset(all_data, target == target_lvl & Efficacy %in% c(value_minEff95)), 
           aes(x=Coverage, y=optimal_Halflife, color=factor(round(Efficacy*100, digits=0)+1), group=factor(round(Efficacy*100, digits =0)+1))) +
  geom_line(size = 2) +
  # geom_ribbon(aes(ymin=sd_minus,ymax=sd_plus, fill = factor(round(Efficacy*100, digits =0))), alpha=0.2, color = NA) +
  scale_color_manual(values = colors1[2], name = "Efficacy", labels = function(x) paste0(x, "%")) +
  # scale_fill_manual(values = colors1, name = "Efficacy [95%CI]", labels = function(x) paste0(x, "%")) +
  scale_y_continuous(breaks = seq(30,180,30), limits = c(30,180)) +
  scale_x_continuous(breaks = seq(0.5,1,0.1), labels = function(x) scales::percent(x, accuracy = 1)) +
  labs(x = "Coverage", y = "Duration (days)") +
  coord_cartesian(expand = c(0,0)) +
  theme(legend.position = "right",
        panel.border = element_rect(color = "black", size = 1, fill = NA),
        panel.background = element_rect(fill = NA),
        panel.grid = element_line(color = "gray90")) +
  theme(title = element_text(size=11, face="bold"),
        plot.background = element_blank(),
        legend.text=element_text(size=11),
        axis.line = element_line(colour = "black"),
        axis.text=element_text(size=11,colour="black"),
        axis.title=element_text(size=11,face="bold"),
        strip.background = element_rect(fill = "gray95"),
        strip.text = element_text(size=11,face="bold"),
        panel.grid.minor = element_blank(),
        legend.position = "right")
p2b

# Save
dev.copy(png, file = paste0("/scicore/home/penny/nekkab0000/M3TPP/Experiments/",exp,"/Outputs/new_Optimization_target",target_lvl,"_Cov_Eff95.png"),
         width = 2400, height = 1100, res = 300)
dev.off()

############################################################ all targets

# Efficacy
p1aa=ggplot(data= subset(all_data, Coverage %in% c(value_minCov60)), 
           aes(x=Efficacy, y=optimal_Halflife, color=target, group=target)) +
  geom_line(size = 2) +
  geom_point(size = 2) +
  # geom_ribbon(aes(ymin=sd_minus,ymax=sd_plus, fill = factor(round(Efficacy*100, digits =0))), alpha=0.2, color = NA) +
  # scale_color_manual(values = colors1[2], name = "Efficacy", labels = function(x) paste0(x, "%")) +
  # scale_fill_manual(values = colors1, name = "Efficacy [95%CI]", labels = function(x) paste0(x, "%")) +
  # scale_color_viridis(discrete = T, na.translate=FALSE) +
  scale_color_manual(values = colorsX) +
  scale_y_continuous(breaks = seq(30,180,30), limits = c(30,180)) +
  scale_x_continuous(breaks = seq(0.5,1,0.1), labels = function(x) scales::percent(x, accuracy = 1)) +
  labs(x = "Efficacy", y = "Duration (days)", title = "80% Coverage", color = "[%] incidence\nreduction") +
  coord_cartesian(expand = c(0,0)) +
  theme(legend.position = "right",
        panel.border = element_rect(color = "black", size = 1, fill = NA),
        panel.background = element_rect(fill = NA),
        panel.grid = element_line(color = "gray90")) +
  theme(title = element_text(size=11, face="bold"),
        plot.background = element_blank(),
        legend.text=element_text(size=11),
        axis.line = element_line(colour = "black"),
        axis.text=element_text(size=11,colour="black"),
        axis.title=element_text(size=11,face="bold"),
        strip.background = element_rect(fill = "gray95"),
        strip.text = element_text(size=11,face="bold"),
        panel.grid.minor = element_blank(),
        legend.position = "right")
p1aa

# Save
dev.copy(png, file = paste0("/scicore/home/penny/nekkab0000/M3TPP/Experiments/",exp1,"/Outputs/new2_Optimization_alltargets_Eff_Cov80.png"),
         width = 2400, height = 1100, res = 300)
dev.off()

# Efficacy
p1bb=ggplot(data= subset(all_data, Coverage %in% c(value_minCov80)), 
            aes(x=Efficacy, y=optimal_Halflife, color=target, group=target)) +
  geom_line(size = 2) +
  geom_point(size = 2) +
  # geom_ribbon(aes(ymin=sd_minus,ymax=sd_plus, fill = factor(round(Efficacy*100, digits =0))), alpha=0.2, color = NA) +
  # scale_color_manual(values = colors1[2], name = "Efficacy", labels = function(x) paste0(x, "%")) +
  # scale_fill_manual(values = colors1, name = "Efficacy [95%CI]", labels = function(x) paste0(x, "%")) +
  # scale_color_viridis(discrete = T, na.translate=FALSE) +
  scale_color_manual(values = colorsX) +
  scale_y_continuous(breaks = seq(30,180,30), limits = c(30,180)) +
  scale_x_continuous(breaks = seq(0.5,1,0.1), labels = function(x) scales::percent(x, accuracy = 1)) +
  labs(x = "Efficacy", y = "Duration (days)", title = "85% Coverage", color = "[%] incidence\nreduction") +
  coord_cartesian(expand = c(0,0)) +
  theme(legend.position = "right",
        panel.border = element_rect(color = "black", size = 1, fill = NA),
        panel.background = element_rect(fill = NA),
        panel.grid = element_line(color = "gray90")) +
  theme(title = element_text(size=11, face="bold"),
        plot.background = element_blank(),
        legend.text=element_text(size=11),
        axis.line = element_line(colour = "black"),
        axis.text=element_text(size=11,colour="black"),
        axis.title=element_text(size=11,face="bold"),
        strip.background = element_rect(fill = "gray95"),
        strip.text = element_text(size=11,face="bold"),
        panel.grid.minor = element_blank(),
        legend.position = "right")
p1bb

# Save
dev.copy(png, file = paste0("/scicore/home/penny/nekkab0000/M3TPP/Experiments/",exp1,"/Outputs/new2_Optimization_alltargets_Eff_Cov85.png"),
         width = 2400, height = 1100, res = 300)
dev.off()

# Coverage
p2aa=ggplot(data= subset(all_data, Efficacy %in% c(value_minEff80)), 
           aes(x=Coverage, y=optimal_Halflife, color=target, group=target)) +
  geom_line(size = 2) +
  geom_point(size = 2) +
  # geom_ribbon(aes(ymin=sd_minus,ymax=sd_plus, fill = factor(round(Efficacy*100, digits =0))), alpha=0.2, color = NA) +
  # scale_color_manual(values = colors1[2], name = "Efficacy", labels = function(x) paste0(x, "%")) +
  # scale_fill_manual(values = colors1, name = "Efficacy [95%CI]", labels = function(x) paste0(x, "%")) +
  # scale_color_viridis(discrete = T, na.translate=FALSE) +
  scale_color_manual(values = colorsX) +
  scale_y_continuous(breaks = seq(30,180,30), limits = c(30,180)) +
  scale_x_continuous(breaks = seq(0.5,1,0.1), labels = function(x) scales::percent(x, accuracy = 1)) +
  labs(x = "Coverage", y = "Duration (days)", title = "90% Efficacy", color = "[%] incidence\nreduction") +
  coord_cartesian(expand = c(0,0)) +
  theme(legend.position = "right",
        panel.border = element_rect(color = "black", size = 1, fill = NA),
        panel.background = element_rect(fill = NA),
        panel.grid = element_line(color = "gray90")) +
  theme(title = element_text(size=11, face="bold"),
        plot.background = element_blank(),
        legend.text=element_text(size=11),
        axis.line = element_line(colour = "black"),
        axis.text=element_text(size=11,colour="black"),
        axis.title=element_text(size=11,face="bold"),
        strip.background = element_rect(fill = "gray95"),
        strip.text = element_text(size=11,face="bold"),
        panel.grid.minor = element_blank(),
        legend.position = "right")
p2aa

# Save
dev.copy(png, file = paste0("/scicore/home/penny/nekkab0000/M3TPP/Experiments/",exp1,"/Outputs/new2_Optimization_alltargets_Cov_Eff90.png"),
         width = 2400, height = 1100, res = 300)
dev.off()




# Coverage
p2bb=ggplot(data= subset(all_data, Efficacy %in% c(value_minEff95)), 
           aes(x=Coverage, y=optimal_Halflife, color=target, group=target)) +
  geom_line(size = 2) +
  geom_point(size = 2) +
  # geom_ribbon(aes(ymin=sd_minus,ymax=sd_plus, fill = factor(round(Efficacy*100, digits =0))), alpha=0.2, color = NA) +
  # scale_color_manual(values = colors1[2], name = "Efficacy", labels = function(x) paste0(x, "%")) +
  # scale_fill_manual(values = colors1, name = "Efficacy [95%CI]", labels = function(x) paste0(x, "%")) +
  # scale_color_viridis(discrete = T, na.translate=FALSE) +
  scale_color_manual(values = colorsX) +
  scale_y_continuous(breaks = seq(30,180,30), limits = c(30,180)) +
  scale_x_continuous(breaks = seq(0.5,1,0.1), labels = function(x) scales::percent(x, accuracy = 1)) +
  labs(x = "Coverage", y = "Duration (days)", title = "95% Efficacy", color = "[%] incidence\nreduction") +
  coord_cartesian(expand = c(0,0)) +
  theme(legend.position = "right",
        panel.border = element_rect(color = "black", size = 1, fill = NA),
        panel.background = element_rect(fill = NA),
        panel.grid = element_line(color = "gray90")) +
  theme(title = element_text(size=11, face="bold"),
        plot.background = element_blank(),
        legend.text=element_text(size=11),
        axis.line = element_line(colour = "black"),
        axis.text=element_text(size=11,colour="black"),
        axis.title=element_text(size=11,face="bold"),
        strip.background = element_rect(fill = "gray95"),
        strip.text = element_text(size=11,face="bold"),
        panel.grid.minor = element_blank(),
        legend.position = "right")
p2bb

# Save
dev.copy(png, file = paste0("/scicore/home/penny/nekkab0000/M3TPP/Experiments/",exp1,"/Outputs/new2_Optimization_alltargets_Cov_Eff95.png"),
         width = 2400, height = 1100, res = 300)
dev.off()


