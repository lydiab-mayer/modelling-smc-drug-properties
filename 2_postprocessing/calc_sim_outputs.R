##############################
# calc_sim_outputs.R - script which post-processes the OpenMalaria 
#                      simulation result and calculates necessary outputs
#
# INPUTS:
#       om_results_folder: path of the folder containing the OpenMalaria simulation output files
#       split_file: path of the file containing the scenario parameter specifications for each scenario to be processed
#       dest_dir: path of the folder where output files will be saved
#
# OUTPUTS:
#   In the folder specified at dest_dir, the script will create two files:
#       - one file with results aggregated across seeds (agg_*.*)
#       - one file with results for all seeds (seed_*.*)
#
# created on 04.01.2019
# monica.golumbeanu@unibas.ch
##############################

library(dplyr)    
library(rapportools)

source("~/smc_lai/analysis_workflow/2_postprocessing/postprocessing_resources.R")

##### Main part of script: #####

# Only for testing and developing:
#om_results_folder = "/scicore/home/smith/GROUP/smc_lai/E2_LAI/om/"
#split_file = "/scicore/home/smith/GROUP/smc_lai/E2_LAI/postprocessing_5/split/E2_LAI_Sen_4.9167_hill_0.5.txt"
#dest_dir = "/scicore/home/smith/GROUP/smc_lai/E2_LAI/postprocessing_5/"
#follow_up = 5

# Read in command arguments
args = commandArgs(TRUE)
om_results_folder = args[1]
split_file = args[2]
dest_dir = args[3]
follow_up = as.integer(args[4])

# Create output file names
split_name = basename(split_file)
dest_table_agg = paste(dest_dir, "agg_", split_name, sep="")
dest_table_seeds = paste(dest_dir, "seeds_", split_name, sep="")

# Postprocess the OpenMalaria simulations
postprocess_OM(results_folder = om_results_folder, param_table_file = split_file,
               final_table_dest = dest_table_agg, final_seed_table_dest = dest_table_seeds,
               follow_up)
