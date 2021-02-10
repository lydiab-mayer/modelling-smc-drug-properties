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

source("/scicore/home/smith/laagmi01/smc_lai/analysis_workflow/3_GP_train/GP_toolbox.R")

args = commandArgs(TRUE)
seeds_file = args[1]
results_folder = args[2]
predicted = args[3]
ranges_file = args[4]


input_data = read.table(seeds_file, sep="\t",header = T, as.is	= T)

exp_name = tools::file_path_sans_ext(basename(seeds_file))
cv_file = paste(results_folder, exp_name, "_cv.RData", sep="")
ranges = load(ranges_file)

if((predicted %in%colnames(input_data)) == FALSE) {
    stop(paste("Column ", predicted, "not found.", sep=""))
}

print(head(input_data))
print(rownames(param_ranges))
print(predicted)

# 5-fold cross-validation
# ml: param ranges file contains parameters that are currently not varied
# hardcode parameter choice for now. come back to this later.

input_parameters <- rownames(param_ranges)
input_parameters <- input_parameters[c(1,2,4,5)]

#split input_data into train and test
n_seeds <- length(unique(input_data$seed))
print(c("n_seeds: ", n_seeds))

input_data$id <- rep(seq(1,nrow(input_data)/n_seeds),each = n_seeds)
n_points = round(nrow(input_data)/n_seeds*0.8)
print(c("n_points: ", n_points))


index_train = sample(nrow(input_data)/n_seeds, n_points)
index_test = setdiff(1:nrow(input_data)/n_seeds, index_train)

train_data <-input_data[input_data$id %in% index_train,]
test_data <- input_data[input_data$id %in% index_test,]

print(c(nrow(train_data),nrow(test_data)))
print(head(train_data[, c(input_parameters, predicted)]))


#train GP
GP_model = train_GP_matern(train_data[, c(input_parameters, predicted)])
cv_result <- list(GP_model = GP_model)

#store output
cv_result$input_parameters <- input_parameters
cv_result$predicted <- predicted

cv_result$train_data <- train_data
cv_result$test_data <- test_data
cv_result$seasonality = unique(input_data$Seasonality)
cv_result$biting_pattern = unique(input_data$Biting_pattern)


save(cv_result, file = cv_file)


