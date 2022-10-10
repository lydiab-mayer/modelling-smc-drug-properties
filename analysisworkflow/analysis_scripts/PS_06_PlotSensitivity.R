############################################################
# PS_06_PlotSensitivity
#
# Visualises relationships between emulator input and predictor variables
# Note: This script depends on outputs of the script 4_sensitivityanalysis_workflow.R
#
# Written by Lydia Braunack-Mayer
############################################################

# ----------------------------------------------------------
# Setup
# ----------------------------------------------------------

rm(list = ls())

# !!! Insert your experiment name here as a string, e.g. "MyExperiment" !!!
exp <- "..."

# !!! Insert your predicted parameters here. Note that this must match with one column name in post-processing files !!!
pred_list <- "..."

library(ggplot2)
library(dplyr)
library(tidyr)
library(ggalluvial)
# library(wesanderson) # Load for access to nice colour palettes

user <- strsplit(getwd(), "/", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]][5]
GROUP_dr <- "/scicore/home/penny/GROUP/M3TPP/"
setwd(paste0("/scicore/home/penny/", user, "/M3TPP"))

if (!dir.exists(paste0("./Experiments/", exp, "/Outputs"))) dir.create(paste0("./Experiments/", exp, "/Outputs"))

load(paste0(GROUP_dr, exp, "/param_ranges.RData"))
param_ranges_cont


# ----------------------------------------------------------
# Define data to plot
# ----------------------------------------------------------

# Import settings
setting <- Sys.glob(paste0(GROUP_dr, exp, "/gp/trained/sensitivity/seeds*"))
setting <- setting[grepl(pred_list, setting, fixed = TRUE)]
setting_id <- sub(paste0("_", pred_list, ".*"), "", sub(".*seeds_", "", setting))

# Import total effect sizes for each setting
df <- data.frame("S_eff" = c(), "T_eff" = c(), scenario = c())

for (i in 1:length(setting)) {

  load(setting[i]) #loads list called sobol_idx_list
  
  sobol_idx_list <- as.data.frame(sobol_idx_list)
  
  sobol_idx_list$S_eff <- sobol_idx_list$S_eff/sum(sobol_idx_list$S_eff) # rescale so total = 1
  sobol_idx_list$T_eff <- sobol_idx_list$T_eff/sum(sobol_idx_list$T_eff) # rescale so total = 1
  
  sobol_idx_list$scenario <- setting_id[i]
  sobol_idx_list$parameter <- rownames(param_ranges_cont)
  
  df <- rbind(df, sobol_idx_list)
  
}

# Import median impact for each setting
df_impact <- data.frame()

for (i in 1:length(setting_id)) {
  
  temp <- read.table(paste0(GROUP_dr, exp, "/postprocessing/agg_", setting_id[i], ".txt"), header = TRUE, sep = "")
  temp <- temp %>%
    group_by(Seasonality, Biting_pattern, EIR, MaxAge, Decay, Access, Timing) %>%
    summarise(inc_red_int_Med = median(inc_red_int_Avg),
              sev_red_int_Med = median(sev_red_int_Avg),
              mor_red_int_Med = median(mor_red_int_Avg))
  temp$scenario <- setting_id[i]
  
  df_impact <- rbind(df_impact, as.data.frame(temp))
  
}

# Scale total effects by median impact for each setting
df <- merge(df, df_impact, by = "scenario")
df$T_eff_scaled <- df$T_eff * df$inc_red_int_Med

# ----------------------------------------------------------
# Format data for plotting - THIS MUST BE ADJUSTED FOR EACH EXPERIMENT
# ----------------------------------------------------------
# 
# df$parameter <- as.factor(df$parameter); df$scenario <- as.factor(df$scenario)
# df <- df %>%
#   separate(col = scenario, 
#            into = c("Experiment", "Seasonality", "System", "EIR", "Agegroup", "Decay", "Access", "Timing"),
#            sep = "_",
#            remove = FALSE)
# 
# df$Seasonality <- ifelse(df$Seasonality == "sharpseasonal", "SHORT SEASON", "LONG SEASON")
# df$Access <- ifelse(df$Access == 0.04, "LOW", "HIGH")
# df$Agegroup <- paste0("CHILDREN 3M TO ", df$Agegroup, "Y")
# df$EIR <- ifelse(df$EIR == 1, "LOW",
#                  ifelse(df$EIR == 4, "MODERATE", "HIGH"))
# df$parameter <- ifelse(df$parameter == "Coverage1", "Coverage of intervention cohort",
#                        ifelse(df$parameter == "Coverage2", "Coverage of each round of treatment",
#                               ifelse(df$parameter == "Efficacy", "Intervention initial efficacy",
#                                      "Intervention duration of protection halflife")))
# 
# df$Seasonality <- factor(df$Seasonality, levels = c("SHORT SEASON", "LONG SEASON"))
# df$Access <- factor(df$Access, levels = c("LOW", "HIGH"))
# df$EIR <- factor(df$EIR, levels = c("LOW", "MODERATE", "HIGH"))


# ----------------------------------------------------------
# Define plot settings
# ----------------------------------------------------------

# !!! Define name of your plot !!!
plot_name <- "..."

cols <- c(...)

# ----------------------------------------------------------
# Generate plot
# ----------------------------------------------------------

p <- ggplot(df, aes(x = ..., y = T_eff_scaled, alluvium = parameter))

p <- p + geom_alluvium(aes(fill = parameter), colour = "white", alpha = 1, decreasing = FALSE)

p <- p + geom_flow(fill = "white", alpha = 0.5)

p <- p + facet_grid(... ~ ... + ...)

p <- p + theme(panel.border = element_blank(), 
               panel.background = element_blank(),
               panel.grid = element_blank(),
               text = element_text(family = "Courier"),
               strip.background = element_blank(),
               axis.line = element_blank(),
               axis.ticks = element_blank(),
               axis.text.y = element_blank(),
               axis.title = element_text(face="bold"),
               legend.key = element_blank(),
               legend.title = element_text(face = "bold"))

p <- p + scale_fill_manual(values = cols)

p <- p + labs(x = "...",
              y = "RELATIVE IMPORTANCE",
              title = plot_name,
              fill = "PARAMETER")

p

ggsave(filename= paste0("./Experiments/", exp, "/Outputs/PS_06_", plot_name, ".jpg"),
       plot = last_plot(),
       width = 15,
       height = 9,
       dpi = 800)

