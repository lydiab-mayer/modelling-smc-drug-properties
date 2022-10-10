############################################################
# DS_01_PrepareOptimisation
#
# Imports, combines and exports optimisation results
# Note: This script depends on outputs of the script 6_grid_optimization_workflow_original.R
#
# Written by Lydia Braunack-Mayer
############################################################

# ----------------------------------------------------------
# Setup
# ----------------------------------------------------------

rm(list = ls())

# !!! Insert your experiment name here as a string, e.g. "MyExperiment" !!!
exp <- "iTPP3_tradeoffs_4rounds"

# !!! Insert your predicted parameters here. Note that this must match with one column name in post-processing files !!!
pred <- c("inc_red_int_Tot")

# Load packages
library(dplyr)
library(tidyr)

# Set working directory
user <- strsplit(getwd(), "/", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]][5]
GROUP_dr <- "/scicore/home/penny/GROUP/M3TPP/"
setwd(paste0("/scicore/home/penny/", user, "/M3TPP"))


# ----------------------------------------------------------
# Import data
# ----------------------------------------------------------

# Import settings
setting <- Sys.glob(paste0(GROUP_dr, exp, "/gp/GP_grid_optimization/", pred, "/opt*"))[1:3]
setting_id <- sub(paste0("_", pred, ".rds"), "", basename(setting))

# Import optimisation results for each setting
out <- data.frame()

for (i in 1:length(setting)) {
  
  # Read data
  df <- readRDS(setting[i])
  df <- df$scenarios
  
  # Format target_range
  df$target_range <- floor(df$mean)
  
  # Add scenario id
  df$scenario <- setting_id[i]
  
  # Expand scenario id into separate variables
  df <- df %>%
    separate(col = scenario, 
             into = c("Opt", "Experiment", "Seasonality", "System", "EIR", "Agegroup", "Decay", "Access", "Timing"),
             sep = "_",
             remove = FALSE)
  
  out <- rbind(out, df)
}


# ----------------------------------------------------------
# Write outputs to file
# ----------------------------------------------------------

write.csv(out, paste0(GROUP_dr, exp, "/gp/GP_grid_optimization/", pred, "/opt_", pred, ".csv"))
write.csv(df, paste0(GROUP_dr, exp, "/gp/GP_grid_optimization/", pred, "/opt_", pred, unique(df$scenario), ".csv"))
