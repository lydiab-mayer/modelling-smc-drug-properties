##################################
# Run the plotting for various models/settings
#
# created 02.12.2019
# monica.golumbeanu@swisstph.ch
##################################

# Load the necessary plotting functions
source("~/MMC/TPP/scripts_v38/plotting/plot_GP_performance.R")
library(stringr)
library(ggplot2)

# Retrieve the command argments and define function inputs
# args = commandArgs(TRUE)
# sim_dir = args[1]
# test_dir = args[2]
# follow_up = args[3]
# plot_title = args[4]

# For testing:
sim_dir = "~/MMC/TPP/simulations/PEV_twice_3years/"
test_dir = "~/MMC/TPP/simulations/test_sets/test_PEV_twice_3years/"
follow_up = 6
plot_title = "PEV twice 3 years"

training_dir = paste0(sim_dir, "gp_", follow_up, "/trained/")
as_dir = paste0(sim_dir, "gp_", follow_up, "/as/")
train_processing_dir = paste0(sim_dir, "postprocessing_", follow_up, "/")
test_processing_dir = paste0(test_dir, "postprocessing_", follow_up, "/")
param_ranges_file = paste0(sim_dir, "param_ranges.RData")
param_table = read.table(paste0(sim_dir, "param_tab.txt"), header = TRUE)
model_pattern = str_remove(param_table$Scenario_Name[1], "_1")
model_name = paste0("seeds_", model_pattern)
model_test_name = paste0("seeds_test_", model_pattern)
plot_dir = paste0("~/MMC/TPP/figures/gp_performance/", model_pattern, "/")

## PLOTS GENERATION
dir.create(plot_dir, showWarnings = FALSE)
# Plot performance on training set (no averaging of seeds)
# plot_performance_GP(training_dir, plot_dir, plot_title, follow_up, model_name)

# Plot performance on training set (seed averaging)
plot_performance_GP_test(training_dir, train_processing_dir, param_ranges_file, plot_dir, plot_title, 
                         paste0("train_", follow_up), model_name)

# Plot performance on the test set (seed averaging)
plot_performance_GP_test(training_dir, test_processing_dir, param_ranges_file, plot_dir, plot_title, 
                         paste0("test_", follow_up), model_test_name)

# Plot performance after adaptive sampling on training set (seed averaging)
plot_performance_GP_test(as_dir, train_processing_dir, param_ranges_file, plot_dir, plot_title, 
                         paste0("train_as_", follow_up), model_name)

# Plot performance after adaptive sampling on test set (seed averaging)
plot_performance_GP_test(as_dir, test_processing_dir, param_ranges_file, plot_dir, plot_title, 
                         paste0("test_as_", follow_up), model_test_name)

# plot_dir = "~/MMC/TPP/figures/gp_performance/MAB_once_3_years/"
# plot_performance_GP("~/MMC/TPP/simulations/MAB_once_3years_avg_prev/gp_4/trained/",
#                     plot_dir, "MAB once 3 years", "4", "seeds_MAB_once_3_years")
# plot_performance_GP("~/MMC/TPP/simulations/MAB_once_3years_avg_prev/gp_6/trained/",
#                     plot_dir, "MAB once 3 years", "6", "seeds_MAB_once_3_years")
# 
# # # Plot performance on test set (averaging of seeds)
# plot_performance_GP_test("~/MMC/TPP/simulations/MAB_once_3years_avg_prev/gp_4/trained/",
#                          "~/MMC/TPP/simulations/test_MAB_once_3years_avg_prev/postprocessing_4/",
#                          "~/MMC/TPP/simulations/MAB_once_3years_avg_prev/param_ranges.RData",
#                          plot_dir, "MAB once 3 years", "test_4", "seeds_MAB_once_3_years")
# plot_performance_GP_test("~/MMC/TPP/simulations/MAB_once_3years_avg_prev/gp_6/trained/",
#                          "~/MMC/TPP/simulations/test_MAB_once_3years_avg_prev/postprocessing_6/",
#                          "~/MMC/TPP/simulations/MAB_once_3years_avg_prev/param_ranges.RData",
#                          plot_dir, "MAB once 3 years", "test_6", "seeds_MAB_once_3_years")
# 
# # Plot the cv results on training set
# plot_performance_GP_test("~/MMC/TPP/simulations/MAB_once_3years_avg_prev/gp_4/trained/",
#                          "~/MMC/TPP/simulations/MAB_once_3years_avg_prev/postprocessing_4/",
#                          "~/MMC/TPP/simulations/MAB_once_3years_avg_prev/param_ranges.RData",
#                          plot_dir, "MAB once 3 years", "train_4", "seeds_MAB_once_3_years")
# plot_performance_GP_test("~/MMC/TPP/simulations/MAB_once_3years_avg_prev/gp_6/trained/",
#                          "~/MMC/TPP/simulations/MAB_once_3years_avg_prev/postprocessing_6/",
#                          "~/MMC/TPP/simulations/MAB_once_3years_avg_prev/param_ranges.RData",
#                          plot_dir, "MAB once 3 years", "train_6", "seeds_MAB_once_3_years")
# 
# # Plot the adaptive sampling results on test set
# plot_performance_GP_test("~/MMC/TPP/simulations/MAB_once_3years_avg_prev/gp_4/as/",
#                          "~/MMC/TPP/simulations/test_MAB_once_3years_avg_prev/postprocessing_4/",
#                          "~/MMC/TPP/simulations/MAB_once_3years_avg_prev/param_ranges.RData",
#                          plot_dir, "MAB once 3 years", "test_as_4", "seeds_MAB_once_3_years")
# plot_performance_GP_test("~/MMC/TPP/simulations/MAB_once_3years_avg_prev/gp_6/as/",
#                          "~/MMC/TPP/simulations/test_MAB_once_3years_avg_prev/postprocessing_6/",
#                          "~/MMC/TPP/simulations/MAB_once_3years_avg_prev/param_ranges.RData",
#                          plot_dir, "MAB once 3 years", "test_as_6", "seeds_MAB_once_3_years")
# 
# # Plot the adaptive sampling results on training set
# plot_performance_GP_test("~/MMC/TPP/simulations/MAB_once_3years_avg_prev/gp_4/as/",
#                          "~/MMC/TPP/simulations/MAB_once_3years_avg_prev/postprocessing_4/",
#                          "~/MMC/TPP/simulations/MAB_once_3years_avg_prev/param_ranges.RData",
#                          plot_dir, "MAB once 3 years", "train_as_4", "seeds_MAB_once_3_years")
# plot_performance_GP_test("~/MMC/TPP/simulations/MAB_once_3years_avg_prev/gp_6/as/",
#                          "~/MMC/TPP/simulations/MAB_once_3years_avg_prev/postprocessing_6/",
#                          "~/MMC/TPP/simulations/MAB_once_3years_avg_prev/param_ranges.RData",
#                          plot_dir, "MAB once 3 years", "train_as_6", "seeds_MAB_once_3_years")
# 
