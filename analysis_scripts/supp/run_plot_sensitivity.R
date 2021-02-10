##########################
# call the plotting functions for sensitivity analysis results
#
# created 06.12.2019
# monica.golumbeanu@swisstph.ch
##########################

# Load the necessary plotting functions
source("~/MMC/TPP/scripts_v38/plotting/plot_GP_sensitivity.R")
library(stringr)

# Retrieve the command argments and define function inputs
# args = commandArgs(TRUE)
# sim_dir = args[1]
# follow_up = args[2]
# plot_title = args[3]

# For testing
sim_dir = "~/MMC/TPP/simulations/PEV_twice_3years/"
follow_up = 6
plot_title = "PEV twice 3 years"

sens_dir = paste0(sim_dir, "gp_", follow_up, "/sensitivity/")
param_ranges_file = paste0(sim_dir, "param_ranges.RData")
param_table = read.table(paste0(sim_dir, "param_tab.txt"), header = TRUE)
model_pattern = str_remove(param_table$Scenario_Name[1], "_1")
plot_dir = paste0("~/MMC/TPP/figures/sensitivity_analysis/", model_pattern, "/")

dir.create(plot_dir, showWarnings = FALSE)
plot_sens_GP(sens_dir, plot_dir, param_ranges_file, plot_title, paste0(model_pattern, "_", follow_up))

# plot_sens_GP("~/MMC/TPP/simulations/MAB_once_3years_avg_prev/gp_4/sensitivity/", 
#              plot_dir, "~/MMC/TPP/simulations/MAB_once_3years_avg_prev/param_ranges.RData",
#              "", "MAB_once_3years_avg_prev_4")
# plot_sens_GP("~/MMC/TPP/simulations/MAB_once_3years_avg_prev/gp_6/sensitivity/",
#              plot_dir, "~/MMC/TPP/simulations/MAB_once_3years_avg_prev/param_ranges.RData",
#              "", "MAB_once_3years_avg_prev_6")
