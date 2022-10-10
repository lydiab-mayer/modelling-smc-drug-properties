############################################################
# PS_07_PlotOptimisation
#
# Visualises relationships between emulator input and optimised variables
# Note: This script depends on outputs of the script 5_optimisation_workflow.R
#
# Written by Lydia Braunack-Mayer
############################################################

# ----------------------------------------------------------
# Setup
# ----------------------------------------------------------

rm(list = ls())

# !!! Insert your experiment name here as a string, e.g. "MyExperiment" !!!
exp <- "iTPP3_tradeoffs_4rounds"

# !!! Insert your predicted parameter here. Note that this must match with one column name in post-processing files !!!
pred <- "sev_red_int_Tot"

library(ggplot2)
library(tidyr)
library(dplyr)
library(scales)

user <- strsplit(getwd(), "/", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]][5]
GROUP_dr <- "/scicore/home/penny/GROUP/M3TPP/"
setwd(paste0("/scicore/home/penny/", user, "/M3TPP/"))

if (!dir.exists(paste0("./Experiments/", exp, "/Outputs"))) dir.create(paste0("./Experiments/", exp, "/Outputs"))


# ----------------------------------------------------------
# Define data to plot
# ----------------------------------------------------------

setting <- Sys.glob(paste0(GROUP_dr, exp, "/gp/GP_grid_optimization/", pred, "/*"))
(setting_id <- sub(".rds", "", sub("opt_", "", basename(setting))))

out <- data.frame()

for (i in 1:length(setting)) {
  print(paste0("Reading data for setting ", setting_id[i]))
  df <- readRDS(setting[i])$scenarios
  df$scenario <- sub(paste0("_", pred), "", setting_id[i])
  out <- rbind(out, df)
}

# out <- out %>%
#   separate(col = scenario,
#            into = c("Experiment", "Seasonality", "System", "EIR", "Agegroup", "Decay", "Access", "Timing"),
#            sep = "_",
#            remove = FALSE)

saveRDS(out, paste0("./Experiments/", exp, "/Outputs/", exp, "_", pred, "_opt_data.rds"))

