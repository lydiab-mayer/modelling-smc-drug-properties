#############################
# Sensitivity analysis for GP prevalence reduction 
#
# created 28.11.2018
# monica.golumbeanu@unibas.ch
#############################


user <- strsplit(getwd(), "/", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]][5]
source(paste0("/scicore/home/penny/",user,"/M3TPP/analysisworkflow/3_GP_train/GP_toolbox.R"))


args = commandArgs(TRUE)
gp_file = args[1]
ranges_file = args[2]
results_folder = args[3]


exp_name = tools::file_path_sans_ext(basename(gp_file))
sidx_file = paste(results_folder, exp_name, "_sidx.RData", sep="")

# Load the GP and parameter ranges
gp_result_name = load(gp_file)
gp_result =  cv_result$GP_model
load(ranges_file)
param_ranges <- param_ranges_cont

# Calculate the Sobol indices
sobol_idx_list = calc_sobol_idx(gp_result, param_ranges, num_points = 50000)


save(sobol_idx_list, file = sidx_file)

