####################################
####################################
###                              ###
### Analysis of number of seeds  ###
###                              ###
####################################
####################################


### -------------------------------------------------------------------------
###
### M3TPP PROJECT:
### Analysis of range of stochastic variations between seeds
###
### Objective: decide the seed number needed for experiments
###
### Conclusion 1: 5 seeds is sufficient for most experiments since the 
### mean of the outcomes variables was very similar to >5 seeds
###
### Conclusion 2: for final results, 10 or 15 seeds is suitable to minimize 
### or approximate the variation in the SD of the mean values as observed
### in 50 seeds
###
### Results: presentation of plots saved in
### J:\DM\DMID\4310_BMGF_MalariaImmunoTPP\3. Data & Analysis\Seed analysis\
###
### Saved 08.10.2021
### narimane.nekkab@unibas.ch
###
### R version 3.6.0
###
### -------------------------------------------------------------------------


##############
### HEADER ###
##############

# Clear environment
rm(list = ls())

# Set seed for replication
set.seed(42)

# Library
library(dplyr)

# User 
user = strsplit(getwd(), "/", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]][5]

# Working directory
setwd(paste0("/scicore/home/penny/",user,"/M3TPP"))


##################
### EXPERIMENT ###
##################

# Insert experiment name here
exp ="MOCK_seed_test"

# Postprocessing file
pp_file = "/scicore/home/penny/GROUP/M3TPP/MOCK_seed_test/postprocessing/seeds_MOCKseedtest_seas4mo_Equal_15_10_exp_0.1.txt"
pp_results_50seeds = read.table(file = pp_file, header = T, stringsAsFactors = F)
# Select variables
pp_results_50seeds = pp_results_50seeds[,c(1,25:ncol(pp_results_50seeds))]

# Select seeds
pp_results_5seeds = pp_results_50seeds %>% filter(seed %in% seq(1,5))
pp_results_10seeds = pp_results_50seeds %>% filter(seed %in% seq(1,10))
pp_results_15seeds = pp_results_50seeds %>% filter(seed %in% seq(1,15))
pp_results_20seeds = pp_results_50seeds %>% filter(seed %in% seq(1,20))
pp_results_25seeds = pp_results_50seeds %>% filter(seed %in% seq(1,25))
pp_results_30seeds = pp_results_50seeds %>% filter(seed %in% seq(1,30))
pp_results_40seeds = pp_results_50seeds %>% filter(seed %in% seq(1,40))

################
### VARIANCE ###
################

############################ 5 seeds

# Summary statistics
pp_results_5seeds_summary = pp_results_5seeds %>% 
  group_by(Scenario_Name, Coverage, Efficacy, Halflife, seed) %>% 
  gather(key = "outcome", value = "value", -c(Scenario_Name, Coverage, Efficacy, Halflife, seed)) %>% 
  select(-seed) %>% 
  group_by(Scenario_Name, Coverage, Efficacy, Halflife, outcome) %>% 
  summarise(mean = mean(value), sd = sd(value)) %>% 
  ungroup() %>% 
  mutate(n_seed = "5 seeds")

# Visualize variation
ggplot(pp_results_5seeds_summary, 
       aes(x = mean, y = sd, color = Efficacy)) +
  facet_wrap(.~outcome) +
  geom_point() +
  scale_color_viridis()

############################ 10 seeds

# Summary statistics
pp_results_10seeds_summary = pp_results_10seeds %>% 
  group_by(Scenario_Name, Coverage, Efficacy, Halflife, seed) %>% 
  gather(key = "outcome", value = "value", -c(Scenario_Name, Coverage, Efficacy, Halflife, seed)) %>% 
  select(-seed) %>% 
  group_by(Scenario_Name, Coverage, Efficacy, Halflife, outcome) %>% 
  summarise(mean = mean(value), sd = sd(value)) %>% 
  ungroup() %>% 
  mutate(n_seed = "10 seeds")

############################ 15 seeds

# Summary statistics
pp_results_15seeds_summary = pp_results_15seeds %>% 
  group_by(Scenario_Name, Coverage, Efficacy, Halflife, seed) %>% 
  gather(key = "outcome", value = "value", -c(Scenario_Name, Coverage, Efficacy, Halflife, seed)) %>% 
  select(-seed) %>% 
  group_by(Scenario_Name, Coverage, Efficacy, Halflife, outcome) %>% 
  summarise(mean = mean(value), sd = sd(value)) %>% 
  ungroup() %>% 
  mutate(n_seed = "15 seeds")

############################ 20 seeds

# Summary statistics
pp_results_20seeds_summary = pp_results_20seeds %>% 
  group_by(Scenario_Name, Coverage, Efficacy, Halflife, seed) %>% 
  gather(key = "outcome", value = "value", -c(Scenario_Name, Coverage, Efficacy, Halflife, seed)) %>% 
  select(-seed) %>% 
  group_by(Scenario_Name, Coverage, Efficacy, Halflife, outcome) %>% 
  summarise(mean = mean(value), sd = sd(value)) %>% 
  ungroup() %>% 
  mutate(n_seed = "20 seeds")

############################ 25 seeds

# Summary statistics
pp_results_25seeds_summary = pp_results_25seeds %>% 
  group_by(Scenario_Name, Coverage, Efficacy, Halflife, seed) %>% 
  gather(key = "outcome", value = "value", -c(Scenario_Name, Coverage, Efficacy, Halflife, seed)) %>% 
  select(-seed) %>% 
  group_by(Scenario_Name, Coverage, Efficacy, Halflife, outcome) %>% 
  summarise(mean = mean(value), sd = sd(value)) %>% 
  ungroup() %>% 
  mutate(n_seed = "25 seeds")

############################ 30 seeds

# Summary statistics
pp_results_30seeds_summary = pp_results_30seeds %>% 
  group_by(Scenario_Name, Coverage, Efficacy, Halflife, seed) %>% 
  gather(key = "outcome", value = "value", -c(Scenario_Name, Coverage, Efficacy, Halflife, seed)) %>% 
  select(-seed) %>% 
  group_by(Scenario_Name, Coverage, Efficacy, Halflife, outcome) %>% 
  summarise(mean = mean(value), sd = sd(value)) %>% 
  ungroup() %>% 
  mutate(n_seed = "30 seeds")


############################ 40 seeds

# Summary statistics
pp_results_40seeds_summary = pp_results_40seeds %>% 
  group_by(Scenario_Name, Coverage, Efficacy, Halflife, seed) %>% 
  gather(key = "outcome", value = "value", -c(Scenario_Name, Coverage, Efficacy, Halflife, seed)) %>% 
  select(-seed) %>% 
  group_by(Scenario_Name, Coverage, Efficacy, Halflife, outcome) %>% 
  summarise(mean = mean(value), sd = sd(value)) %>% 
  ungroup() %>% 
  mutate(n_seed = "40 seeds")

############################ 50 seeds

# Summary statistics
pp_results_50seeds_summary = pp_results_50seeds %>% 
  group_by(Scenario_Name, Coverage, Efficacy, Halflife, seed) %>% 
  gather(key = "outcome", value = "value", -c(Scenario_Name, Coverage, Efficacy, Halflife, seed)) %>% 
  select(-seed) %>% 
  group_by(Scenario_Name, Coverage, Efficacy, Halflife, outcome) %>% 
  summarise(mean = mean(value), sd = sd(value)) %>% 
  ungroup() %>% 
  mutate(n_seed = "50 seeds")

# Visualize variation
ggplot(pp_results_50seeds_summary, 
       aes(x = mean, y = sd, color = Efficacy)) +
  facet_wrap(.~outcome) +
  geom_point() +
  scale_color_viridis()

############################ Combine data

# Combine dataframes
pp_results_all_seed_summaries = rbind(pp_results_5seeds_summary,
                                      pp_results_10seeds_summary,
                                      pp_results_15seeds_summary,
                                      pp_results_20seeds_summary,
                                      pp_results_25seeds_summary,
                                      pp_results_30seeds_summary,
                                      pp_results_40seeds_summary,
                                      pp_results_50seeds_summary)

# Prep for plotting
pp_results_all_seed_summaries = pp_results_all_seed_summaries %>% 
  mutate(n_seed = factor(n_seed, levels = c("5 seeds","10 seeds","15 seeds",
                                            "20 seeds","25 seeds","30 seeds",
                                            "40 seeds","50 seeds")))

# Visualize variation: all points
ggplot(pp_results_all_seed_summaries, 
       aes(x = mean, y = sd, group = n_seed, color = n_seed)) +
  facet_wrap(.~outcome) +
  geom_point() +
  scale_color_viridis(discrete = T, alpha = 0.5)

# Visualize variation
ggplot(pp_results_all_seed_summaries, 
       aes(x = n_seed, y = sd, color = n_seed)) +
  facet_wrap(.~outcome) +
  geom_point() +
  scale_color_viridis(discrete = T) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

# Ridgeline
library(ggridges)
ggplot(pp_results_all_seed_summaries, 
       aes(x = mean, y = n_seed, fill = n_seed)) +
  facet_wrap(.~outcome) +
  geom_density_ridges_gradient() +
  theme_ridges() +
  scale_fill_viridis(discrete = T)

# Save
dev.copy(png, file = paste0("/scicore/home/penny/nekkab0000/M3TPP/Experiments/",exp,"/Mean_density_ridgeline.png"),
         width = 3000, height = 2000, res = 300)
dev.off()

# Ridgeline
ggplot(pp_results_all_seed_summaries, 
       aes(x = sd, y = n_seed, fill = n_seed)) +
  facet_wrap(.~outcome) +
  geom_density_ridges_gradient() +
  theme_ridges() +
  scale_fill_viridis(discrete = T)

# Save
dev.copy(png, file = paste0("/scicore/home/penny/nekkab0000/M3TPP/Experiments/",exp,"/SD_density_ridgeline.png"),
         width = 3000, height = 2000, res = 300)
dev.off()

# # Save
# dev.copy(png, file = paste0("/scicore/home/penny/nekkab0000/M3TPP/Experiments/",exp,"/SD_per_nseed_per_scenario.png"),
#          width = 3000, height = 2000, res = 300)
# dev.off()


# Visualize variation: loess smoothing
ggplot(pp_results_all_seed_summaries, 
       aes(x = mean, y = sd, group = n_seed, color = n_seed)) +
  facet_wrap(.~outcome) +
  # geom_point() +
  scale_color_viridis(discrete = T) +
  geom_smooth(se = T)

# Save
dev.copy(png, file = paste0("/scicore/home/penny/nekkab0000/M3TPP/Experiments/",exp,"/SD_per_nseed_per_scenario_loess.png"),
         width = 3000, height = 2000, res = 300)
dev.off()

# Visualize variation: histogram
ggplot(pp_results_all_seed_summaries, 
       aes(x = mean, fill = outcome)) +
  facet_grid(n_seed~outcome) +
  geom_histogram(position="identity") +
  scale_fill_viridis(discrete = T) +
  theme(legend.position = "none")

# Save
dev.copy(png, file = paste0("/scicore/home/penny/nekkab0000/M3TPP/Experiments/",exp,"/meanDistr_per_nseed_per_scenario_histogram.png"),
         width = 3000, height = 2000, res = 300)
dev.off()

# Visualize variation: histogram
ggplot(pp_results_all_seed_summaries, 
       aes(x = sd, fill = outcome)) +
  facet_grid(n_seed~outcome) +
  geom_histogram(position="identity") +
  scale_fill_viridis(discrete = T) +
  theme(legend.position = "none")

# Save
dev.copy(png, file = paste0("/scicore/home/penny/nekkab0000/M3TPP/Experiments/",exp,"/sdDistr_per_nseed_per_scenario_histogram.png"),
         width = 3000, height = 2000, res = 300)
dev.off()

############################ Summary

pp_results_all_seed_summaries_mean = pp_results_all_seed_summaries %>% 
  group_by(outcome, n_seed) %>% 
  summarise(avg_mean = mean(mean), avg_sd = mean(sd)) %>% 
  ungroup()

# Visualize variation: mean
ggplot(pp_results_all_seed_summaries_mean, 
       aes(x = avg_mean, y = n_seed)) +
  facet_wrap(.~outcome) +
  geom_point() +
  scale_color_viridis(discrete = T, alpha = 0.25) +
  labs(x = "average means") 

# Save
dev.copy(png, file = paste0("/scicore/home/penny/nekkab0000/M3TPP/Experiments/",exp,"/avgMeans_per_nseed.png"),
         width = 3000, height = 2000, res = 300)
dev.off()


# Visualize variation: sd
ggplot(pp_results_all_seed_summaries_mean, 
       aes(x = avg_sd, y = n_seed)) +
  facet_wrap(.~outcome) +
  geom_point() +
  scale_color_viridis(discrete = T, alpha = 0.25) +
  labs(x = "average sds") 

# Save
dev.copy(png, file = paste0("/scicore/home/penny/nekkab0000/M3TPP/Experiments/",exp,"/avgSDs_per_nseed.png"),
         width = 3000, height = 2000, res = 300)
dev.off()

####################### variation by grouped outcome values

pp_results_all_seed_summaries_outcomeCat = pp_results_all_seed_summaries %>% 
  mutate(mean_cat = cut(mean, breaks = seq(0,100,5))) %>% 
  select(-c(mean, Scenario_Name, Coverage, Efficacy, Halflife)) %>% 
  group_by(outcome, n_seed, mean_cat) %>%
  summarise(mean = mean(sd), sd=sd(sd)) %>% 
  ungroup()

# Visualize variation: mean
ggplot(pp_results_all_seed_summaries_outcomeCat, 
       aes(x = mean_cat, y = mean, group = n_seed, color = n_seed)) +
  facet_wrap(.~outcome) +
  geom_point() +
  geom_line() +
  scale_color_viridis(discrete = T) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  labs(y = "sd")

# Save
dev.copy(png, file = paste0("/scicore/home/penny/nekkab0000/M3TPP/Experiments/",exp,"/avgSDs_outcomeCat_per_nseed.png"),
         width = 3000, height = 2000, res = 300)
dev.off()
  
# Visualize variation: sd
ggplot(pp_results_all_seed_summaries_outcomeCat, 
       aes(x = mean_cat, y = mean, group = n_seed, color = n_seed)) +
  facet_grid(n_seed~outcome) +
  geom_point() +
  geom_line() +
  scale_color_viridis(discrete = T) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  labs(y = "sd") +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd))

# Save
dev.copy(png, file = paste0("/scicore/home/penny/nekkab0000/M3TPP/Experiments/",exp,"/avgSDs_outcomeCat_per_nseed_errorbars.png"),
         width = 3000, height = 2000, res = 300)
dev.off()

###############################
### VARIANCE FOR 1 SCENARIO ###
###############################

# Postprocessing file
pp_file = "/scicore/home/penny/GROUP/M3TPP/MOCK_seed_test/postprocessing/seeds_MOCKseedtest_seas4mo_Equal_15_10_exp_0.1.txt"
pp_results_50seeds = read.table(file = pp_file, header = T, stringsAsFactors = F)
# Select variables
pp_results_50seeds = pp_results_50seeds[,c(1,25:ncol(pp_results_50seeds))]

# List of scenarios
scenarios = unique(pp_results_50seeds$Scenario_Name)

# Pick 1 scenario
pp_results_50seeds_s1 = pp_results_50seeds %>% filter(Scenario_Name %in% scenarios[1:5])

# Select seeds
pp_results_5seeds_s1 = pp_results_50seeds_s1 %>% filter(seed %in% seq(1,5)) %>% mutate(n_seed = "5 seeds")
pp_results_10seeds_s1 = pp_results_50seeds_s1 %>% filter(seed %in% seq(1,10)) %>% mutate(n_seed = "10 seeds")
pp_results_15seeds_s1 = pp_results_50seeds_s1 %>% filter(seed %in% seq(1,15)) %>% mutate(n_seed = "15 seeds")
pp_results_20seeds_s1 = pp_results_50seeds_s1 %>% filter(seed %in% seq(1,20)) %>% mutate(n_seed = "20 seeds")
pp_results_25seeds_s1 = pp_results_50seeds_s1 %>% filter(seed %in% seq(1,25)) %>% mutate(n_seed = "25 seeds")
pp_results_30seeds_s1 = pp_results_50seeds_s1 %>% filter(seed %in% seq(1,30)) %>% mutate(n_seed = "30 seeds")
pp_results_40seeds_s1 = pp_results_50seeds_s1 %>% filter(seed %in% seq(1,40)) %>% mutate(n_seed = "40 seeds")
pp_results_50seeds_s1 = pp_results_50seeds_s1 %>% filter(seed %in% seq(1,50)) %>% mutate(n_seed = "50 seeds")

# Combine dataframes
pp_results_all_seeds_s1 = rbind(pp_results_5seeds_s1,
                                pp_results_10seeds_s1,
                                pp_results_15seeds_s1,
                                pp_results_20seeds_s1,
                                pp_results_25seeds_s1,
                                pp_results_30seeds_s1,
                                pp_results_40seeds_s1,
                                pp_results_50seeds_s1)

# Summary statistics
pp_results_all_seeds_s1_summary = pp_results_all_seeds_s1 %>% 
  select(-c(seed)) %>% 
  group_by(Scenario_Name, n_seed) %>% 
  gather(key = "outcome", value = "value", -c(Scenario_Name, Coverage, Efficacy, Halflife, n_seed)) %>% 
  group_by(Scenario_Name, n_seed, outcome) %>% 
  summarise(mean = mean(value), sd = sd(value)) %>% 
  ungroup() %>% 
  mutate(n_seed = factor(n_seed, levels = c("5 seeds","10 seeds","15 seeds",
                                            "20 seeds","25 seeds","30 seeds",
                                            "40 seeds","50 seeds")))

# Visualize variation
ggplot(pp_results_all_seeds_s1_summary, 
       aes(x = mean, y = sd, color = n_seed)) +
  facet_grid(Scenario_Name~outcome) +
  geom_point() +
  scale_color_viridis(discrete = T)

# Save
dev.copy(png, file = paste0("/scicore/home/penny/nekkab0000/M3TPP/Experiments/",exp,"/SD_mean_per_nseed_5scenarios.png"),
         width = 3000, height = 2000, res = 300)
dev.off()

# Visualize variation
ggplot(pp_results_all_seeds_s1_summary, 
       aes(x = n_seed, y = sd, color = n_seed)) +
  facet_grid(Scenario_Name~outcome) +
  geom_point() +
  scale_color_viridis(discrete = T) +
  theme(axis.text.x = element_text(angle = 90, hjust = 0.5))

# Save
dev.copy(png, file = paste0("/scicore/home/penny/nekkab0000/M3TPP/Experiments/",exp,"/SD_by_nseed_5scenarios.png"),
         width = 3000, height = 2000, res = 300)
dev.off()
