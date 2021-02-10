#############################
# Sensitivity analysis for GP prevalence reduction 
#
# created 28.11.2018
# monica.golumbeanu@unibas.ch
#############################

source("/scicore/home/penny/burlyd00/smc_lai/analysis_workflow/3_GP_train/GP_toolbox.R")

args = commandArgs(TRUE)
gp_file = args[1]
ranges_file = args[2]
results_folder = args[3]
# print(cv_file)
# print(ranges_file)
# print(EIR_fixed_lvl)
# print(results_folder)

# For testing:
 # gp_file = "/scicore/home/penny/GROUP/smc_lai/E3E4disc_comparison/gp_trained/gp_trained_pppy_LAI_Mali_4.9167_hill_0.1_150.RData"
 #ranges_file = "/scicore/home/penny/GROUP/smc_lai/E4_SMCLAIdisc/param_ranges.RData"
#  results_folder = "/scicore/home/penny/GROUP/smc_lai/E3E4_lb/sensitivity/"


exp_name = tools::file_path_sans_ext(basename(gp_file))
sidx_file = paste(results_folder, exp_name, "_sidx_lowls.RData", sep="")

# Load the GP and parameter ranges
gp_result_name = load(gp_file)
gp_result = get(gp_result_name)
load(ranges_file)

param_col <- c("Coverage_SMC", "Efficacy","Halflife","Coverage")

param_ranges <- param_ranges[param_col, ]
param_ranges[3,] <- c(30,90)
# Calculate the Sobol indices
sobol_idx_list = calc_sobol_idx(gp_result, param_ranges, num_points = 500000)

sobol_idx_list$seasonality = strsplit(exp_name,"_")[[1]][5]
sobol_idx_list$LAI_dec = strsplit(exp_name,"_")[[1]][7]
sobol_idx_list$Access = strsplit(exp_name,"_")[[1]][8]
sobol_idx_list$EIR = strsplit(exp_name,"_")[[1]][9]

save(sobol_idx_list, file = sidx_file)

