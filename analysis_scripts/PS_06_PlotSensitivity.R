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
exp <- "iTPP3_tradeoffs"

# !!! Insert your predicted parameters here. Note that this must match with one column name in post-processing files !!!
pred_list <- "inc_red_int_Avg"

library(ggplot2)
library(dplyr)
library(tidyr)
library(ggalluvial)

user <- strsplit(getwd(), "/", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]][5]
GROUP_dr <- "/scicore/home/penny/GROUP/M3TPP/"
setwd(paste0("/scicore/home/penny/", user, "/M3TPP"))

if (!dir.exists(paste0("./Experiments/", exp, "/Outputs"))) dir.create(paste0("./Experiments/", exp, "/Outputs"))

load(paste0(GROUP_dr, exp, "/param_ranges.RData"))
param_ranges_cont


# ----------------------------------------------------------
# Define data to plot
# ----------------------------------------------------------

setting <- Sys.glob(paste0(GROUP_dr, exp, "/gp/trained/sensitivity/seeds*"))
setting <- setting[grepl(pred_list, setting, fixed = TRUE)]
setting_id <- sub(paste0("_", pred_list, ".*"), "", sub(".*seeds_", "", setting))

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

df$parameter <- as.factor(df$parameter); df$scenario <- as.factor(df$scenario)
df <- df %>%
  separate(col = scenario, 
           into = c("Experiment", "Seasonality", "System", "EIR", "Agegroup", "Decay", "Access", "Timing"),
           sep = "_",
           remove = FALSE)

df$Seasonality <- ifelse(df$Seasonality == "sharpseasonal", "Short season", "Long season")
df$Access <- ifelse(df$Access == 0.04, "Low", "High")
df$Agegroup <- paste0("Children 3m to ", df$Agegroup, "y")
df$EIR <- ifelse(df$EIR == 1, "Low",
                 ifelse(df$EIR == 4, "Moderate", "High"))

df$Seasonality <- factor(df$Seasonality, levels = c("Short season", "Long season"))
df$Access <- factor(df$Access, levels = c("Low", "High"))
df$EIR <- factor(df$EIR, levels = c("Low", "Moderate", "High"))

# ----------------------------------------------------------
# Define plot settings
# ----------------------------------------------------------

# !!! Define name of your plot !!!
plot_name <- "Test plot 2 for sensitivity analysis"


# ----------------------------------------------------------
# Generate plot
# ----------------------------------------------------------

p <- ggplot(df, aes(x = EIR, y = T_eff, alluvium = parameter))

p <- p + geom_alluvium(aes(fill = parameter), colour = "white", alpha = 1, decreasing = FALSE)

p <- p + facet_grid(Agegroup ~ Access + Seasonality)

p <- p + theme(panel.border = element_blank(), 
               panel.background = element_blank(),
               panel.grid = element_blank(),
               text = element_text(colour = "#323f4f", family = "Simplon Norm (Body)"),
               strip.background = element_blank(),
               axis.line = element_blank(),
               axis.ticks = element_blank(),
               axis.title = element_text(face="bold"),
               legend.key = element_blank())

p <- p + scale_y_continuous(labels = scales::percent_format(accurary = 1L))

p <- p + labs(x = "Scenario",
              y = "Importance",
              title = plot_name,
              fill = "Parameter")

p

ggsave(filename= paste0("./Experiments/", exp, "/Outputs/PS_06_", plot_name, ".jpg"),
       plot = last_plot(),
       width = 15,
       height = 9,
       dpi = 800)

