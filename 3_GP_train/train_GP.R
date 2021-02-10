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

source("/scicore/home/smith/burlyd00/smc_lai/analysis_workflow/3_GP_train/GP_toolbox.R")
library(sjPlot)

args = commandArgs(TRUE)
split_file = args[1]
results_folder = args[2]
predicted = args[3]
ranges_file = args[4]

# Only for testing:
  # split_file = "/scicore/home/smith/GROUP/smc_lai/E2_LAI/postprocessing_5/seeds_E2_LAI_Sen_4.9167_hill_0.5.txt"
  # results_folder = "/scicore/home/smith/GROUP/smc_lai/E2_LAI/gp_5/trained/"
  # predicted = "incred"
  # ranges_file = "/scicore/home/smith/GROUP/smc_lai/E2_LAI/param_ranges.RData"

# Only for testing:
 # split_file = "/scicore/home/smith/burlyd00/smc_lai/E3E4_comp/postprocessing_5/seeds_E3E4_comp_Sen_4.9167_hill_0.5.txt"
 # results_folder = "/scicore/home/smith/burlyd00/smc_lai/E3E4_comp/gp_5/trained/"
 # predicted = "UL"
 # ranges_file = "/scicore/home/smith/burlyd00/smc_lai/E3E4_comp/param_ranges.RData"



OM_result = read.table(split_file, sep="\t", header = TRUE, as.is = TRUE)
exp_name = tools::file_path_sans_ext(basename(split_file))
cv_file = paste(results_folder, exp_name,"_",predicted, "_cv.RData", sep="")
ranges = load(ranges_file)

if(predicted=="incred") {
  OM_result$incred <- ((OM_result$clin_y1/OM_result$mean_popint) -  OM_result$pppy_y10_all)/( OM_result$clin_y1/OM_result$mean_popint) 
  OM_result$incred <-  ifelse(OM_result$incred<0,0,OM_result$incred)
}
  


if((predicted %in%colnames(OM_result)) == FALSE) {
    stop(paste("Column ", predicted, "not found.", sep=""))
}


# Select only the entries where the prevalence pre-deployment was not 0

to_del <-OM_result[which( OM_result$iprev_y5_1 == 0),"Scenario_Name"]
if (length(to_del) >0) {input_data = OM_result[-which(OM_result$Scenario_Name  %in%  to_del), ]
} else {input_data = OM_result
}

# Select relevant columns from dataframe such that the last column is the predicted one:
#param_ranges <- param_ranges[-which(param_ranges[,1]==param_ranges[,2]),]

# 5-fold cross-validation

input_parameters <- rownames(param_ranges)

#split input_data into train and test
n_seeds <- length(unique(input_data$seed))

if(n_seeds==5.5) {n_seeds=1}
print(c("n_seeds: ", n_seeds))

input_data$id <- rep(seq(1,nrow(input_data)/n_seeds),each = n_seeds)
n_points = round(nrow(input_data)/n_seeds*0.8)
print(c("n_points: ", n_points))


index_train = sample(nrow(input_data)/n_seeds, n_points)
index_test = setdiff(1:nrow(input_data)/n_seeds, index_train)

train_data <-input_data[input_data$id %in% index_train,]
test_data <- input_data[input_data$id %in% index_test,]

print(c(nrow(train_data),nrow(test_data)))


#train GP
GP_model = train_GP_matern(train_data[, c(input_parameters, predicted)])

test_data2 =   test_data %>% group_by(Scenario_Name) %>% summarise_all(funs(if(is.numeric(.)) median(., na.rm = TRUE) else unique(.)))
test_GP_plot(GP_model=GP_model, test=test_data2[,c(rownames(param_ranges), predicted)],save=paste(results_folder, exp_name,"_",predicted, "_R2.jpg", sep=""))



cv_result <- list(GP_model = GP_model)

#store output
cv_result$input_parameters <- input_parameters
cv_result$predicted <- predicted

cv_result$train_data <- train_data
cv_result$test_data <- test_data
cv_result$seasonality = unique(input_data$Seasonality)
cv_result$LAI_dec = unique(input_data$Decay_Scen)
cv_result$Access = unique(input_data$Access)
cv_result$IntAge = unique(input_data$IntAge)


save(cv_result, file = cv_file)


# testing of GP 
# input_data$id <- rep(seq(1,nrow(input_data)/seeds),each=seeds)
# 
# seeds = 5
# num_points = round(nrow(input_data)/seeds*0.75)
# 
# index_train = sample(nrow(input_data)/seeds, num_points)
# index_test = setdiff(1:nrow(input_data)/seeds, index_train)
# 
# 
# train_data <-input_data[input_data$id %in% index_train,]
# test_data <- input_data[input_data$id %in% index_test,]
# test_data2 =   test_data %>% group_by(Scenario_Name) %>% summarise_at(c(names(test_data)[which(names(input_data)=="kdecay"):length(names(test_data) ) ])
#                                                                       ,mean,na.rm=TRUE)
# 
# trained_model = train_GP(train_data[ , c(rownames(param_ranges), predicted)] ) 
# test_model<- test_GP(GP_model=trained_model, train_data=train_data[,c(rownames(param_ranges), predicted)], test_data=test_data2[,c(rownames(param_ranges), predicted)])
# test_GP_plot(GP_model=trained_model, test=test_data2[,c(rownames(param_ranges), predicted)]) 
# 
# 
# 
# train_model_matern <- train_GP_matern(train_data[ ,c(rownames(param_ranges), predicted)]) 
# test_model_matern  <- test_GP(GP_model=train_model_matern, train_data=train_data[,c(rownames(param_ranges), predicted)],
#                               test_data=test_data2[,c(rownames(param_ranges), predicted)])
# test_GP_plot(GP_model=train_model_matern, test=test_data2[,c(rownames(param_ranges), predicted)]) 
# 
# 
