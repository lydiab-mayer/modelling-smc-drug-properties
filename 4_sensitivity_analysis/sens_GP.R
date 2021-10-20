#############################
# Sensitivity analysis for GP prevalence reduction 
#
# created 28.11.2018
# monica.golumbeanu@unibas.ch
#
# Updated September 2021
# lydia.braunack-mayer@swisstph.ch
#############################


user <- strsplit(getwd(), "/", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]][5]
source(paste0("/scicore/home/penny/",user,"/M3TPP/analysisworkflow/3_GP_train/GP_toolbox.R"))


args = commandArgs(TRUE)
gp_file = args[1]
ranges_file = args[2]
results_folder = args[3]
scale = args[4]

# Sample arguments, retained here for testing
gp_file <- "/scicore/home/penny/GROUP/M3TPP/iTPP3_tradeoffs/gp/trained/inc_red_int_Avg/seeds_iTPP3tradeoffs_sharpseasonal_Mali_15_10_exp_0.241193660515256_May_inc_red_int_Avg_cv.RData"
ranges_file <- "/scicore/home/penny/GROUP/M3TPP/iTPP3_tradeoffs/param_ranges.RData"
results_folder <- "/scicore/home/penny/GROUP/M3TPP/iTPP3_tradeoffs/gp/trained/"
scale <- TRUE

exp_name = tools::file_path_sans_ext(basename(gp_file))
sidx_file = paste(results_folder, exp_name, "_sidx.RData", sep="")

# Load the GP and parameter ranges
gp_result_name = load(gp_file)
gp_result =  cv_result$GP_model
load(ranges_file); param_ranges <- param_ranges_cont

# If inputs have been scaled to c(0, 1) when training the emulator, update parameter ranges to match
if (scale) {
  param_ranges[, 1] <- 0 
  param_ranges[, 2] <- 1 
}

# Calculate the Sobol indices
sobol_idx_list = calc_sobol_idx(GP_model = gp_result, param_spec = param_ranges, num_points = 50000)


save(sobol_idx_list, file = sidx_file)