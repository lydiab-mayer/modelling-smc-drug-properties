#################################################################
#################################################################
######
###### Plotting time series of incidence and prevalence
###### Created 26.07.2021
###### Narimane Nekkab
###### narimane.nekkab@swisstph.ch
######
###### Adapted from 
###### lydia.burgert@unibas.ch 
######
#################################################################
#################################################################

rm(list = ls())

#################
### LIBRARIES ###

library(tgp)
library(dplyr)
library(ggrepel)
detach(package:plyr)
library(zoo)

######################
### SETUP & PARAMS ###

set.seed(42)

# User 
user = strsplit(getwd(), "/", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]][5]

# Working directory and group directory
setwd(paste0("/scicore/home/penny/",user,"/M3TPP"))
GROUP_dr = "/scicore/home/penny/GROUP/M3TPP/"

# Source function scripts
source(paste0("./analysisworkflow/analysis_scripts/supp/import_functions.R"))

# insert experiment name here
# exp ="MOCK_EIR8_sig_april_age10"
# exp ="MOCK_EIR15_sig_june_age10"
# exp ="MOCK_EIR15_sig_may_age10"
# exp ="MOCK_EIR15_sig2_may_age10"
exp ="MOCK_EIR15_sig2_april_age10"

# Import parameter table for experiment 
param_table_file <- paste0(GROUP_dr,exp,"/param_tab.txt")
param_table <- read.table(param_table_file, sep= "\t", header = TRUE, as.is = TRUE, stringsAsFactors = FALSE)
# Coverage & Eff combos
unique(param_table$Coverage)
unique(param_table$Efficacy)

# choose which categories to plot 
load(paste0(GROUP_dr,exp,"/param_ranges_cat.RData"))

# define follow up and year before intervention
follow_up=5
years_before_interv=5

# Set coverage and efficacy min. values for plotting
minCov = 0.78
maxCov = 0.82
minEff = 0.88
maxEff = 0.92

# minCov = 0.68
# maxCov = 0.72
# minEff = 0.83
# maxEff = 0.87

# minCov = 0.68
# maxCov = 0.72
# minEff = 0.88
# maxEff = 0.92

##############
### SUBSET ###

# Subset
param_table = param_table[which(param_table$Coverage > minCov &
                                  param_table$Coverage < maxCov &
                                  param_table$Efficacy > minEff &
                                  param_table$Efficacy < maxEff &
                                  param_table$Seasonality == "seas4mo"),]

###################
### DATE UPLOAD ###

# Loop through all scenarios
all_data_lists = foreach(i = 1:length(unique(param_table$Scenario_Name))) %do% {
  
  # set simulation to plot 
  scenario_name <- unique(param_table$Scenario_Name)[i] #input scenario name from param_table$Scenario_Name column
  
  # Get data
  results <- import_cont_4var(exp,
                         scenario_name,
                         param_table_file,
                         seeds=max(param_table$SEED),
                         follow_up,
                         years_before_interv)
  
  # Combine
  results = results
  
}

# Cut & relabel
prevalence_210_per=lapply(all_data_lists, `[[`, 1)
prevalence_210_per=dplyr::bind_rows(prevalence_210_per, .id = "sim")
prevalence_int_per=lapply(all_data_lists, `[[`, 2)
prevalence_int_per=dplyr::bind_rows(prevalence_int_per, .id = "sim")
incidence_int_per=lapply(all_data_lists, `[[`, 3)
incidence_int_per=dplyr::bind_rows(incidence_int_per, .id = "sim")
incidence_int_5mo_per=lapply(all_data_lists, `[[`, 4)
incidence_int_5mo_per=dplyr::bind_rows(incidence_int_5mo_per, .id = "sim")

# Add cov/eff/half info
all_values = NULL
for(i in 1:nrow(prevalence_210_per)){
  sim = as.numeric(prevalence_210_per$sim[i])
  scenario = unique(param_table$Scenario_Name)[sim]
  values = unique(param_table[which(param_table$Scenario_Name == scenario),c("Coverage","Halflife","Efficacy")])
  all_values[[i]] = values
}
all_values = do.call(rbind, all_values)

# Merge
prevalence_210_per = cbind(prevalence_210_per, all_values)
prevalence_int_per = cbind(prevalence_int_per, all_values)
incidence_int_per = cbind(incidence_int_per, all_values)
incidence_int_5mo_per = cbind(incidence_int_5mo_per, all_values)

###########################
### DATA TRANSFORMATION ###

##############################
# Moving average of incidence

# Smoothing distance 
distance = 3

# Moving average
prevalence_int_per_mavg_minCov_minEff = prevalence_int_per %>% 
  # Filter
  filter(Coverage > minCov & Coverage < maxCov &
           Efficacy > minEff & Efficacy < maxEff) %>% 
  # Mean of 10 sims
  dplyr::select(-sim) %>% 
  group_by(time,year,Coverage,Halflife,Efficacy) %>% 
  summarise_all(list(mean)) %>% 
  # New incidence value
  group_by(year,Coverage,Halflife,Efficacy) %>% 
  arrange(time) %>% 
  mutate(year_2 = 2029 + year) %>% 
  group_by(Coverage,Halflife,Efficacy) %>% 
  mutate(time_2 = (time-min(time))*5/30,
         prev_ma = rollmean(prev, distance, na.pad=TRUE, align="right"),
         Halflife_round = round(as.numeric(Halflife),digits=-1)) %>% 
  ungroup()

# Moving average
prevalence_210_per_mavg_minCov_minEff = prevalence_210_per %>% 
  # Filter
  filter(Coverage > minCov & Coverage < maxCov &
           Efficacy > minEff & Efficacy < maxEff) %>% 
  # Mean of 10 sims
  dplyr::select(-sim) %>% 
  group_by(time,year,Coverage,Halflife,Efficacy) %>% 
  summarise_all(list(mean)) %>% 
  # New incidence value
  group_by(year,Coverage,Halflife,Efficacy) %>% 
  arrange(time) %>% 
  mutate(year_2 = 2029 + year) %>% 
  group_by(Coverage,Halflife,Efficacy) %>% 
  mutate(time_2 = (time-min(time))*5/30,
         prev_ma = rollmean(prev, distance, na.pad=TRUE, align="right"),
         Halflife_round = round(as.numeric(Halflife),digits=-1)) %>% 
  ungroup()

# Moving average
incidence_int_per_mavg_minCov_minEff = incidence_int_per %>% 
  # Filter
  filter(Coverage > minCov & Coverage < maxCov &
           Efficacy > minEff & Efficacy < maxEff) %>% 
  # Mean of 10 sims
  dplyr::select(-sim) %>% 
  group_by(time,year,Coverage,Halflife,Efficacy) %>% 
  summarise_all(list(mean)) %>% 
  # New incidence value
  group_by(Coverage,Halflife,Efficacy) %>% 
  arrange(time) %>% 
  mutate(time_2 = (time-min(time))*5/30,
         cpp_ma = rollmean(cpp, distance, na.pad=TRUE, align="right"),
         Halflife_round = round(as.numeric(Halflife),digits=-1)) %>% 
  ungroup()


################################################################### 
############## PREVALENCE PLOT

# Function
getPlot0 = function(data, subset_list, month1, month1name, month2name){
  
  # Subset
  subset_data = data[which(data$Halflife_round %in% subset_list),]
  
  # Order 
  subset_data$Halflife_round_fct = factor(subset_data$Halflife_round, levels = subset_list)
  
  # Add colors
  colors1 = viridis(4)[1:length(subset_list)]
  names(colors1) = subset_list
  
  # Incidence moving average + intervention + axis limits + labs
  p0=ggplot() +
    scale_x_continuous(breaks = sort(c(seq(2,2+16*73,73), # 1 mark per year
                                       seq(2+5*73+month1*6,2+16*73+month1*6,73))-7), # 1 mark for 5 mo
                       limits = c(2+4*73,2+10*73)-7,
                       labels = c("2030","2031","2032","2033","2034",
                                  "2035",month1name,
                                  "2036",month1name,
                                  "2037",month1name,
                                  "2038",month1name,
                                  "2039",month1name,
                                  "2040",month1name,
                                  "2041",month1name,
                                  "2042",month1name,
                                  "2043",month1name,
                                  "2044",month1name,
                                  "2045",month1name,
                                  "2046",month1name)) +
    scale_y_continuous(limits = c(0,100)) +
    geom_vline(xintercept = c(seq(2+5*73+month1*6,2+16*73+month1*6,73)-7),
               color = "red",size=0.7,alpha=0.8,linetype="dashed") +
    geom_vline(xintercept = c(seq(2,2+16*73,73)-7),
               color = "gray70",size=0.7,alpha=0.8,linetype="solid") +
    geom_line(data=subset_data, aes(x=time, y=prev_ma*100, group=Halflife_round_fct, 
                                    color=Halflife_round_fct),size=1.5,alpha=1) +
    geom_hline(yintercept = 0, color = "black", size = 1) +
    geom_vline(xintercept = 2+4*73-7, color = "black", size = 1) +
    coord_cartesian(expand = c(0,0)) +
    labs(x = NULL, y = expression(paste(bolditalic("Pf"),bold("PR"["2-10"]),bold(" [%]" ))), color = "Duration (days)") +
    scale_color_manual(values = colors1, drop = T) +
    theme_minimal() +
    theme(panel.grid = element_blank(),
          text = element_text(size = 12),
          axis.text.x = element_text(angle = 90, vjust = 0.5),
          legend.title = element_text(face = "bold"))
  return(p0)
}

# Choose halflives to plot
unique(prevalence_210_per_mavg_minCov_minEff$Halflife)
unique(prevalence_210_per_mavg_minCov_minEff$Halflife_round)
# halflife_list = c(80,100,160)
halflife_list = c(100,130,180)

# Plot
dataX = prevalence_210_per_mavg_minCov_minEff[which(prevalence_210_per_mavg_minCov_minEff$Halflife > 100),]
n = 1
p_0 = getPlot0(dataX, halflife_list[1:n], 4, "April", "September")
p_0

# # Save
# dev.copy(png, file = paste0("/scicore/home/penny/nekkab0000/M3TPP/Experiments/",exp,"/Outputs/new_Timeseries_4mo_may_minCov_minEff_prev210yo_",n,".png"),
#          width = 2500, height = 1100, res = 300)
# dev.off()

################################################################### 
############## INCIDENCE PLOTS

# Function
getPlot1_old = function(data, subset_list, month1, month1name, month2name){
  
  # Subset
  subset_data = data[which(data$Halflife_round %in% subset_list),]
  
  # Order 
  subset_data$Halflife_round_fct = factor(subset_data$Halflife_round, levels = subset_list)
  
  # Add colors
  colors1 = viridis(4)[1:length(subset_list)]
  names(colors1) = subset_list
  
  # Incidence moving average + intervention + axis limits + labs
  p1=ggplot() +
    scale_x_continuous(breaks = sort(c(seq(2,2+16*73,73), # 1 mark per year
                                       seq(2+5*73+month1*6,2+16*73+month1*6,73), # 1 mark per april post-2035
                                       2+9*73+(month1+5)*6)-7), # 1 mark for 5 mo
                       limits = c(2+4*73,2+10*73)-7,
                       labels = c("2030","2031","2032","2033","2034",
                                  "2035",month1name,
                                  "2036",month1name,
                                  "2037",month1name,
                                  "2038",month1name,
                                  "2039",month1name,month2name,
                                  "2040",month1name,
                                  "2041",month1name,
                                  "2042",month1name,
                                  "2043",month1name,
                                  "2044",month1name,
                                  "2045",month1name,
                                  "2046",month1name)) +
    # scale_y_continuous(limits = c(0,100)) +
    geom_vline(xintercept = c(seq(2+5*73+month1*6,2+16*73+month1*6,73)-7),
               color = "red",size=0.7,alpha=0.8,linetype="dashed") +
    geom_vline(xintercept = c(2+9*73+(month1+5)*6)-7,
               color = "blue",size=0.7,alpha=0.8,linetype="dashed") +
    geom_vline(xintercept = c(seq(2,2+16*73,73)-7),
               color = "gray70",size=0.7,alpha=0.8,linetype="solid") +
    # geom_hline(yintercept = c(seq(0,100,25)),
    #            color = "gray70",size=0.7,alpha=0.8,linetype="solid") +
    geom_line(data=subset_data, aes(x=time, y=cpp_ma, group=Halflife_round_fct, 
                                    color=Halflife_round_fct),size=1.5,alpha=1) +
    geom_hline(yintercept = 0, color = "black", size = 1) +
    geom_vline(xintercept = 2+4*73-7, color = "black", size = 1) +
    coord_cartesian(expand = c(0,0)) +
    labs(x = NULL, y = "Incidence (infections person/year)", color = "Half-life (days)") +
    scale_color_manual(values = colors1, drop = T) +
    theme_minimal() +
    theme(panel.grid = element_blank(),
          text = element_text(size = 12),
          axis.text.x = element_text(angle = 90, vjust = 0.5),
          legend.title = element_text(face = "bold"),
          axis.title.y = element_text(face = "bold"))
  return(p1)
}
getPlot1 = function(data, subset_list, month1, month1name, month2name){
  
  # Subset
  subset_data = data[which(data$Halflife_round %in% subset_list),]
  
  # Order 
  subset_data$Halflife_round_fct = factor(subset_data$Halflife_round, levels = subset_list)
  
  # Add colors
  colors1 = viridis(4)[1:length(subset_list)]
  names(colors1) = subset_list
  
  # Incidence moving average + intervention + axis limits + labs
  p1=ggplot() +
    scale_x_continuous(breaks = sort(c(seq(2,2+16*73,73), # 1 mark per year
                                       seq(2+5*73+month1*6,2+16*73+month1*6,73), # 1 mark per april post-2035
                                       2+9*73+(month1+5)*6)-7), # 1 mark for 5 mo
                       limits = c(2+4*73,2+10*73)-7,
                       labels = c("2030","2031","2032","2033","2034",
                                  "2035",month1name,
                                  "2036",month1name,
                                  "2037",month1name,
                                  "2038",month1name,
                                  "2039",month1name,month2name,
                                  "2040",month1name,
                                  "2041",month1name,
                                  "2042",month1name,
                                  "2043",month1name,
                                  "2044",month1name,
                                  "2045",month1name,
                                  "2046",month1name)) +
    # scale_y_continuous(limits = c(0,100)) +
    geom_vline(xintercept = c(seq(2+5*73+month1*6,2+16*73+month1*6,73)-7),
               color = "red",size=0.7,alpha=0.8,linetype="dashed") +
    geom_vline(xintercept = c(2+9*73+(month1+5)*6)-7,
               color = "blue",size=0.7,alpha=0.8,linetype="dashed") +
    geom_vline(xintercept = c(seq(2,2+16*73,73)-7),
               color = "gray70",size=0.7,alpha=0.8,linetype="solid") +
    # geom_hline(yintercept = c(seq(0,100,25)),
    #            color = "gray70",size=0.7,alpha=0.8,linetype="solid") +
    geom_line(data=subset_data, aes(x=time, y=cpp_ma*1000, group=Halflife_round_fct, 
                                    color=Halflife_round_fct),size=1.5,alpha=1) +
    geom_hline(yintercept = 0, color = "black", size = 1) +
    geom_vline(xintercept = 2+4*73-7, color = "black", size = 1) +
    coord_cartesian(expand = c(0,0)) +
    labs(x = NULL, y = expression(paste(bold("Incidence "),"(infections per 1000 children/5-days)")), color = "Duration (days)") +
    scale_color_manual(values = colors1, drop = T) +
    theme_minimal() +
    theme(panel.grid = element_blank(),
          text = element_text(size = 12),
          axis.text.x = element_text(angle = 90, vjust = 0.5),
          legend.title = element_text(face = "bold"),
          axis.title.y = element_text(size = 10))
  return(p1)
}

# Plots
dataY = incidence_int_per_mavg_minCov_minEff[which(incidence_int_per_mavg_minCov_minEff$Halflife > 100),]
p_1 = getPlot1(dataY, halflife_list[1:n], 4, "April", "September")
p_1

# # Save
# dev.copy(png, file = paste0("/scicore/home/penny/nekkab0000/M3TPP/Experiments/",exp,"/Outputs/new_Timeseries_4mo_may_minCov_minEff_inc_",n,".png"),
#          width = 2500, height = 1100, res = 300)
# dev.off()


################################################################### 
############## DROP in INCIDENCE (1 year eval period)

# Calculate drop in incidence
inc_int_red_per_year = incidence_int_per_mavg_minCov_minEff %>%
  group_by(sim, year) %>%
  arrange(time) %>% 
  mutate(intinc_cum = cumsum(intinc)) %>% 
  ungroup() %>%
  dplyr::select(sim, year, intinc_cum, intn, Halflife_round) %>%
  group_by(sim, year) %>%
  filter(intinc_cum == max(intinc_cum)) %>% 
  ungroup() %>% 
  dplyr::select(-sim) %>% 
  group_by(year, Halflife_round) %>%
  summarise_all(list(mean)) %>%
  ungroup() %>% 
  mutate(inc_perperson_peryear = (intinc_cum / intn)) %>% 
  # Calculate pre-intervention average inc red
  group_by(Halflife_round) %>% 
  mutate(preint = ifelse(year <= 5, "pre","post")) %>% 
  group_by(Halflife_round, preint) %>% 
  mutate(preint_avg_inc = mean(inc_perperson_peryear)) %>% 
  group_by(year, Halflife_round,intinc_cum,intn,inc_perperson_peryear) %>% 
  spread(key = "preint", value = preint_avg_inc) %>% 
  group_by(Halflife_round) %>% 
  mutate(post = max(post), pre = max(pre),
         inc_red_all_years = (pre - post) / pre,
         inc_red_per_year = (pre - inc_perperson_peryear) / pre) %>% 
  ungroup()


################################################################### 
############## DROP in INCIDENCE (5 month eval period)


# Calculate drop in incidence
inc_int_red_per_5mo = incidence_int_per_mavg_minCov_minEff %>%
  filter(Halflife > 100) %>% 
  # Average all values across same sims
  dplyr::select(time, year, Coverage, Halflife, Efficacy, Coverage, intinc, intn, cpp, cppcum) %>% 
  group_by(time, year, Halflife, Coverage, Efficacy) %>% 
  summarise_all(list(mean)) %>% 
  ungroup() %>% 
  # Calculate time year select pre and post intervetion years
  filter(year >= 5) %>% 
  arrange(time) %>% 
  group_by(year,Halflife, Coverage, Efficacy) %>% 
  mutate(pre_int_y1_status = ifelse(year == 5, "pre","post"),
         timeyear = time-min(time),
         post_int_5mo_status = ifelse(timeyear %in% seq(18,48),
                                      "5mo_interval",pre_int_y1_status)) %>% 
  filter(post_int_5mo_status %in% c("5mo_interval")) %>% 
  # Get new cumcpp values
  arrange(timeyear) %>% 
  group_by(Halflife,year) %>% 
  mutate(cumcpp_new = cumsum(cpp)) %>% 
  # Get only last obs per year
  filter(timeyear == 48) %>% 
  mutate(status = ifelse(pre_int_y1_status == "post",paste0("5mo_",as.character(year)),pre_int_y1_status)) %>% 
  dplyr::select(-pre_int_y1_status,-timeyear,-post_int_5mo_status,-intinc,-intn,-cppcum,-cpp,-time) %>% 
  group_by(status, Halflife, Coverage, Efficacy) %>% 
  # summarise_all(list(mean)) %>% 
  spread(key = "status", value = "cumcpp_new") %>% 
  # Calculate reduction in cpp 
  group_by(Halflife, Coverage, Efficacy,year,pre) %>% 
  gather(key = "5mo", value = "cumcpp_new", -c(Halflife, Coverage, Efficacy,year,pre)) %>% 
  group_by(Halflife, Coverage, Efficacy) %>% 
  mutate(pre = max(pre)) %>% 
  filter(!is.na(cumcpp_new)) %>% 
  mutate(reduction = ((pre - cumcpp_new) / pre) * 100) %>% 
  ungroup()


################################################################### 
############## DROP in PREVALENCE (1 year eval period)


# Calculate drop in prev
prev_210_red_per_year = prevalence_210_per_mavg_minCov_minEff %>%
  filter(Halflife > 100) %>% 
  group_by(year) %>%
  arrange(time) %>% 
  mutate(avg_prev1 = mean(prev)) %>% 
  ungroup() %>%
  dplyr::select(year, avg_prev1, Halflife_round) %>%
  distinct() %>% 
  # Calculate pre-intervention average inc red
  group_by(Halflife_round) %>% 
  mutate(preint = ifelse(year <= 5, "pre","post")) %>% 
  group_by(Halflife_round,preint) %>% 
  mutate(avg_prev2 = mean(avg_prev1)) %>%
  spread(key = "preint", value = avg_prev2) %>% 
  group_by(Halflife_round) %>% 
  mutate(post = max(post), pre = max(pre),
         prev_red_all_years = (pre - post) / pre,
         prev_red_per_year = (pre - avg_prev1) / pre) %>% 
  ungroup()


################################################################### 
############## DROP in PREVALENCE (5 month eval period)

# Calculate drop in incidence
prev_210_red_per_5mo = prevalence_210_per_mavg_minCov_minEff %>%
  # Average all values across same sims
  dplyr::select(time, year, Coverage, Halflife, Efficacy, Coverage, prev) %>% 
  group_by(time, year, Halflife, Coverage, Efficacy) %>% 
  summarise_all(list(mean)) %>% 
  ungroup() %>% 
  # Classify time period to calculate pre and post cpp
  filter(year >= 5) %>% 
  arrange(time) %>% 
  group_by(year,Halflife, Coverage, Efficacy) %>% 
  mutate(pre_int_y1_status = ifelse(year == 5, "pre","post"),
         timeyear = time-min(time),
         post_int_5mo_status = ifelse(pre_int_y1_status == "post" & timeyear %in% seq(18,48),
                                      "5mo_interval",pre_int_y1_status)) %>% 
  filter(post_int_5mo_status %in% c("5mo_interval","pre")) %>% 
  mutate(status = ifelse(post_int_5mo_status == "5mo_interval",paste0("5mo_",as.character(year)),post_int_5mo_status)) %>% 
  dplyr::select(-pre_int_y1_status,-timeyear,-post_int_5mo_status,-time) %>% 
  group_by(status, Halflife, Coverage, Efficacy) %>% 
  summarise_all(list(mean)) %>% 
  spread(key = "status", value = "prev") %>% 
  # Calculate reduction in cpp 
  group_by(Halflife, Coverage, Efficacy,year,pre) %>% 
  gather(key = "5mo", value = "yearly_prev", -c(Halflife, Coverage, Efficacy,year,pre)) %>% 
  group_by(Halflife, Coverage, Efficacy) %>% 
  mutate(pre = max(pre)) %>% 
  filter(!is.na(yearly_prev)) %>% 
  mutate(reduction = ((pre - yearly_prev) / pre) * 100) %>% 
  ungroup()
