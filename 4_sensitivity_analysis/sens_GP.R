#############################
# Sensitivity analysis for GP prevalence reduction 
#
# created 28.11.2018
# monica.golumbeanu@unibas.ch
#############################

source("/scicore/home/smith/burlyd00/smc_lai/analysis_workflow/3_GP_train/GP_toolbox.R")

args = commandArgs(TRUE)
gp_file = args[1]
ranges_file = args[2]
EIR_fixed_lvl = args[3]
results_folder = args[4]
# print(cv_file)
# print(ranges_file)
# print(EIR_fixed_lvl)
# print(results_folder)

# For testing:
 # gp_file = "/scicore/home/smith/GROUP/smc_lai/E2_LAI/gp_5/trained/pppy_y10_all/seeds_E2_LAI_Sen_4.9167_hill_0.5_cv.RData"
 # ranges_file = "/scicore/home/smith/GROUP/smc_lai/E2_LAI/param_ranges.RData"
 # results_folder = "/scicore/home/smith/GROUP/smc_lai/E2_LAI/gp_5/sensitivity/"
 # EIR_fixed_lvl = 200

exp_name = tools::file_path_sans_ext(basename(gp_file))
sidx_file = paste(results_folder, exp_name, "_EIR_", EIR_fixed_lvl, "_sidx.RData", sep="")

# Load the GP and parameter ranges
gp_result_name = load(gp_file)
gp_result = get(gp_result_name)
load(ranges_file)
#param_ranges <- param_ranges[-which(param_ranges[,1]==param_ranges[,2]),]

param_ranges["EIR",] = c(EIR_fixed_lvl, EIR_fixed_lvl)

# Calculate the Sobol indices
sobol_idx_list = calc_sobol_idx(gp_result$GP_model, param_ranges, num_points = 50000)
sobol_idx_list$EIR = EIR_fixed_lvl
sobol_idx_list$seasonality = gp_result$seasonality
sobol_idx_list$LAI_dec = gp_result$LAI_dec
sobol_idx_list$Access = gp_result$Access
sobol_idx_list$IntAge = gp_result$IntAge

save(sobol_idx_list, file = sidx_file)

