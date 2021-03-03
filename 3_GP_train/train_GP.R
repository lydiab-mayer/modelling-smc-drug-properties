#############################
# Training GP on a database of simulations
#
# INPUT: 
#       split_file: file with the data points for a setting
#       results_folder: folder where the trained object will be saved
#       predicted: name of the quantity to be predicted (e.g. prev_red for prevalence reduction)
#       ranges_file: file with parameter names and ranges of their values
#
# OUTPUT: an Rdata file comtaining a list with the following attributes:
#           model: the trained GP model
#           EIR: the setting EIR
#           seasonality: the setting seasonality
#           biting_pattern: the mosquito biting pattern
#
# created 28.11.2018
# monica.golumbeanu@unibas.ch
#############################

source("../../../analysisworkflow/3_GP_train/GP_toolbox.R")
library(sjPlot)

args = commandArgs(TRUE)
split_file = args[1]
results_folder = args[2]
predicted = args[3]
ranges_file = args[4]


OM_result = read.table(split_file, sep="\t", header = TRUE, as.is = TRUE)
exp_name = tools::file_path_sans_ext(basename(split_file))
cv_file = paste(results_folder, exp_name,"_",predicted, "_cv.RData", sep="")
ranges = load(ranges_file)
param_ranges <- param_ranges_cont
  
if((predicted %in%colnames(OM_result)) == FALSE) {
    stop(paste("Column ", predicted, "not found.", sep=""))
}


# Select only the entries where the prevalence pre-deployment was not 0

#to_del <-OM_result[which( OM_result$iprev_y5_1 == 0),"Scenario_Name"]
#if (length(to_del) >0) {input_data = OM_result[-which(OM_result$Scenario_Name  %in%  to_del), ]
#} else {input_data = OM_result
#}

input_data <- OM_result

# Select relevant columns from dataframe such that the last column is the predicted one:


# 5-fold cross-validation

input_parameters <- rownames(param_ranges_cont)

#split input_data into train and test
n_seeds <- length(unique(input_data$seed))

print(c("n_seeds: ", n_seeds))


input_data <- input_data[,c(input_parameters,predicted)]
#train GP
cv_result = cv_train_matern(input_data, 5)

save= paste(results_folder, exp_name,"_",predicted, "_cv.jpg", sep="")

test_GP_plot(cv_result,save)

#store output
cv_result$input_parameters <- input_parameters
cv_result$predicted <- predicted

save(cv_result, file = cv_file)
