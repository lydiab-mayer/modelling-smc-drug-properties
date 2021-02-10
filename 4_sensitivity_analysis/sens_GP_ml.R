#############################
# Sensitivity analysis for GP prevalence reduction 
#
# created 28.11.2018
# monica.golumbeanu@unibas.ch
#############################

source("/scicore/home/smith/laagmi01/smc_lai/analysis_workflow/3_GP_train/GP_toolbox.R")

args = commandArgs(TRUE)
gp_file = args[1]
ranges_file = args[2]
EIR_fixed_lvl = args[3]
results_folder = args[4]

print(gp_file)
print(ranges_file)
print(EIR_fixed_lvl)
print(results_folder)

exp_name = tools::file_path_sans_ext(basename(gp_file))
sidx_file = paste(results_folder, exp_name, "_EIR_", EIR_fixed_lvl, "_sidx.RData", sep="")

# Load the GP and parameter ranges
gp_result_name = load(gp_file)
gp_result = get(gp_result_name)
load(ranges_file)
param_ranges["EIR",] = c(EIR_fixed_lvl, EIR_fixed_lvl)

#not all parameters are varied. choose here for now.
#come back to this later
param_ranges <- param_ranges[c("EIR","Coverage","Halflife","Efficacy"),]

# Calculate the Sobol indices
sobol_idx_list = calc_sobol_idx(gp_result$GP_model, param_ranges, num_points = 200000)
sobol_idx_list$EIR = EIR_fixed_lvl
sobol_idx_list$seasonality = gp_result$seasonality
sobol_idx_list$biting_pattern = gp_result$biting_pattern
save(sobol_idx_list, file = sidx_file)

