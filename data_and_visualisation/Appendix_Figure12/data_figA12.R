############################################################
#
# Visualises SMC deployment dynamics
#
# Written by Lydia Braunack-Mayer
############################################################

# ----------------------------------------------------------
# Setup
# ----------------------------------------------------------

# Clear workspace
rm(list = ls())

# !!! Insert your experiment name here as a string, e.g. "MyExperiment" !!!
exp <- "iTPP3_ChemoBlood_TreatLiver_4rounds"

# Load required libraries
library(dplyr)

# Define user and group
user <- strsplit(getwd(), "/", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]][5]
GROUP_dr <- "/scicore/home/penny/GROUP/M3TPP/"

# Set working directory
setwd(paste0("/scicore/home/penny/",user,"/M3TPP"))

# Source helper functions
source(paste0("./analysisworkflow/analysis_scripts/import_functions.R"))

# Create folder to store outputs
if (!dir.exists(paste0("./Experiments/",exp,"/Outputs"))) dir.create(paste0("./Experiments/",exp,"/Outputs"))

# Import parameter table
param_table <- read.table(paste0(GROUP_dr,exp,"/param_tab.txt"), sep = "\t", header = TRUE, as.is = TRUE, stringsAsFactors = FALSE)


# ----------------------------------------------------------
# Define data to plot
# ----------------------------------------------------------

# !!! Define which setting you want to plot !!! Should be specified as integer (e.g. '1') or integer vector (e.g. 'c(1, 3, 5)')
(setting <- unique(param_table[, !(names(param_table) == "SEED")]))
setting_id <- 7:8
setting[setting_id, ]

# !!! Define which random seeds you want to visualise !!! Should be defined as an integer or integer vector
seeds <- 1:5

# !!! Define which time steps you want to visualise !!! Should be defined as an integer vector
timesteps <- 1:73

# Calculate average input vs. simulated EIR across all seeds for this setting
df <- import_EIRs_cat(exp = exp, 
                      scenario_id = setting[setting_id, "Scenario_Name"], 
                      seeds = seeds, 
                      timesteps = timesteps)

c(paste0("Total input EIR: ", sum(df$Average$input.EIR)), paste0("Total simulated EIR: ", sum(df$Average$simulated.EIR)))

df$Average$scenario_id <- as.factor(df$Average$scenario_id)
df$Average$scenario_id <- recode(df$Average$scenario_id, 
                                 "iTPP3_ChemoBlood_TreatLiver_4rounds_7" = "THREE MONTHS",
                                 "iTPP3_ChemoBlood_TreatLiver_4rounds_8" = "FIVE MONTHS")


# ----------------------------------------------------------
# Write data to file
# ----------------------------------------------------------

saveRDS(df, "./SMC_TPP/data_and_visualisation/Appendix_Figure12/data_figA12.rds")
