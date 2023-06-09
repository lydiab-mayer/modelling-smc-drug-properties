#############################
# script train_GP.R
#
# Trains a gaussian process emulator on a database of simulations prepared by postprocessing
# 
#
# INPUT: 
#       split_file: file pathway with the data points for a setting
#       results_folder: folder file pathway where the trained object will be saved
#       predicted: name of the quantity to be predicted (e.g. prev_red_int for prevalence reduction)
#       ranges_file: file with parameter names and ranges of their values
#       lower: lower input parameter bounds passed to mleHetGP() for mle optimisation
#       upper: upper input parameter bounds passed to mleHetGP() for mle optimisation
#       scale: TRUE/FALSE indicator of whether input data should be scaled to c(0, 1)
#
# OUTPUT: 
#       Rdata file containing a list with the following attributes:
#         'train_data', the dataset used to perform gaussian process regression
#         'test_data', the dataset set aside as test data
#         'GP_model', the gaussian process regression model fitted to train_data
#       jpeg plot visualising performance attributes of the gaussian process regression
#
# Created 28.11.2018
# monica.golumbeanu@unibas.ch
#
# Adapted October 2021
# lydia.braunack-mayer@swisstph.ch
#
#############################

##############
### HEADER ###
##############

# Source helper functions
source("../../../analysisworkflow/3_GP_train/GP_toolbox.R")

# Define script inputs
args = commandArgs(TRUE)
split_file = args[1]
results_folder = args[2]
predicted = args[3]
ranges_file = args[4]
lower = args[5]
upper = args[6]
scale = args[7]
hetGP = args[8]

# # Sample script inputs, retained here to facilitate testing
# split_file <- "/scicore/home/penny/GROUP/M3TPP/iTPP3_ChemoBlood_TreatLiver_4rounds/postprocessing/seeds_iTPP3ChemoBloodTreatLiver4rounds_seas5mo_Mali_1_5_0.24_May_0.02083134.txt"
# results_folder <- "/scicore/home/penny/GROUP/M3TPP/iTPP3_ChemoBlood_TreatLiver_4rounds/gp/trained/inc_red_int_Tot/"
# predicted <- "inc_red_int_Tot"
# ranges_file <- "/scicore/home/penny/GROUP/M3TPP/iTPP3_ChemoBlood_TreatLiver_4rounds/param_ranges.RData"
# lower <- "0.001/0.001/0.001/0.001/0.001"
# upper <- "10/10/10/10/10"
# scale <- TRUE
# hetGP <- TRUE


#######################
### SET UP EMULATOR ###
#######################

# Format script inputs
input_data <- read.table(split_file, sep="\t", header = TRUE, as.is = TRUE)
exp_name <- tools::file_path_sans_ext(basename(split_file))
cv_file <- paste(results_folder, exp_name,"_",predicted, "_cv.RData", sep="")
ranges <- load(ranges_file)
param_ranges <- param_ranges_cont
input_parameters <- rownames(param_ranges_cont)
lower <- as.numeric(strsplit(lower, "/")[[1]])
upper <- as.numeric(strsplit(upper, "/")[[1]])

# Check for errors in inputs
if((predicted %in%colnames(input_data)) == FALSE) {
  stop(paste("Column ", predicted, "not found.", sep=""))
}

if(sum(is.na(input_data)) != 0) {
  cat("Warning: One or more inputs are missing, proceeding with complete cases only")
  input_data <- input_data[complete.cases(input_data), ]
}

# Define data used to train the gaussian process regression
input_data <- input_data[, c(input_parameters, predicted)]

# Set up scaling
if (scale == TRUE) {
  scale <- param_ranges 
} else {
  scale <- NULL
}


######################
### TRAIN EMULATOR ###
######################

# Train gaussian process regression
cv_result <- cv_train_matern(input_data = input_data, lower = lower, upper = upper, scale = scale, hetGP = hetGP)

# Generate plot of gaussian process regression performance
save <- paste0(results_folder, exp_name, "_", predicted, "_cv.jpg")
test_GP_plot(cv_result, save)

# Store output
cv_result$input_parameters <- input_parameters
cv_result$predicted <- predicted
save(cv_result, file = cv_file)
