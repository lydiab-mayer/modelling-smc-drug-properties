####################################################
# Adaptive sampling scheme for a multi-dimensional space
#
# 
# monica.golumbeanu@unibas.ch
###################################################

library(hetGP)
library(tgp)
library(dplyr)
library(magrittr)
library(rapportools)
source("~/smc_lai/analysis_workflow/2_postprocessing/postprocessing_resources.R")

# Function that returns a set of points to be explored with OpenMalaria and added to the training set
# Method: uniform sampling across the parameter space
select_points_uniform = function(GP_model, param_ranges, num_points) {
    # select a random set of points uniformly distributed across the parameter space
    X_random_samples = lhs(100000, param_ranges) 
    colnames(X_random_samples) = rownames(param_ranges)
    
    # Predict the output with the GP model on the new points and evaluate the variance   
    pred_obj = predict(x = X_random_samples, object = GP_model)
    new_data = cbind.data.frame(X_random_samples, pred_obj$sd2, row.names = NULL)
    colnames(new_data) = c(colnames(X_random_samples), "variance")
    
    # Select a set of points with the largest posterior predictive variance
    XX = new_data[order(new_data$variance, decreasing = TRUE),]
    return(list(param_tab = XX[1:num_points,], variances_vec = c(min(pred_obj$sd2), mean(pred_obj$sd2), max(pred_obj$sd2))))
}

# Create parameter table from newly sampled points to be used for building corresponding scenarios 
create_param_table_from_samples = function(sampled_points, seasonality, biting,Decay_Scen,maxGroup,Access, n_seeds) {
    # remove variance column
    sampled_points$variance = NULL
    # load seasonality monthly values and mosquito biting patterns
    seasons =read.table("~/smc_lai/analysis_workflow/resource_files/seasonality.txt", sep="\t", header = TRUE)
    #biting_patterns = read.table("~/smc_lai/analysis_workflow/resource_files/biting.txt", sep="\t", header = TRUE)
    biting_pattern <- data.frame(Biting_pattern=c("Mali"),indoor=c(0.6),outdoor=c(0.4))
    # for now: only low indoor biting and seasonal transmission 
    LAIdecay <- data.frame(fundecay=c("weibull","hill"),kdecay=c(1,8 ),Decay_Scen=c("exp","hill" ) )#"weibull", 2
    
    IntAge = data.frame(IntAge=c(4.9167),maxGroup=c(3))#,9.9167 and ,4
    
    # add setting parameters
    as_param_tab = cbind.data.frame(seasons[which(seasons$Seasonality == seasonality),], 
                                    biting_pattern[which(biting_pattern$Biting_pattern == biting),],
                                    LAIdecay[which(LAIdecay$Decay_Scen == Decay_Scen),],
                                    IntAge[which(IntAge$maxGroup == maxGroup),],
                                    Access,
                                    sampled_points, row.names = NULL)
    # add seeds and define scenarios names
    SEED = c(1:n_seeds)
    scenarios_names = paste("Scenario", 1:nrow(as_param_tab), sep="_")
    # construct parameter table
    as_param_tab = cbind.data.frame(scenarios_names, as_param_tab, row.names = NULL)
    colnames(as_param_tab)[1] = "Scenario_Name"
    as_param_tab = merge(as_param_tab, as.data.frame(SEED))
    
    return(as_param_tab)
}

# Function that waits until all jobs have finished (marked in a log file)
wait_for_jobs = function(log_file, n_lines) {
    while((!file.exists(log_file))) {
        # waiting for log file to be created
        Sys.sleep(1)
    }
    while(nrow(read.table(log_file)) < n_lines ) {
        # waiting for scenarios to be created
        Sys.sleep(1)
    }
    return(0)
}

args = commandArgs(TRUE)
as_folder = args[1]
gp_file = args[2]
parameter_ranges = args[3]
follow_up = args[4]
predicted = args[5]

# For testing:
as_folder = "/scicore/home/smith/GROUP/smc_lai/Test_SMC/gp_5/as/"

 gp_file = "/scicore/home/smith/GROUP/smc_lai/Test_SMC/gp_5/trained/agg_E3_E4_comparison_seasonal_Low_indoor_exp_cv.RData"
 parameter_ranges = "/scicore/home/smith/GROUP/smc_lai/Test_SMC/param_ranges.RData"
 scaffold_file1 =  "/scicore/home/smith/burlyd00/smc_lai/E4_SMCLAI/scaffold_old.xml"
 scaffold_file2 =  "/scicore/home/smith/burlyd00/smc_lai/E3_SMCSMC/scaffold_old.xml"
 
 follow_up = 5
 predicted = "prop"

 
 
 as_folder2 = paste0(as_folder, "as_smcsmc/")
as_folder1 = paste0(as_folder, "as_smclai/")
                                         
## Load the GP model and parameter ranges
load(gp_file)
GP_updated = cv_result$GP_model
load(parameter_ranges)

## Initialization
# Set number of seeds for OM runs
n_seeds = 20
# Set number of new sampled points for updating the GP model
n_samples = 600
# Initalize variances vector
variances = NULL
om_results_dir1 = paste0(as_folder1, "om/")
om_scenarios_dir1 = paste0(as_folder1, "scenarios/")

om_results_dir2 = paste0(as_folder2, "om/")
om_scenarios_dir2 = paste0(as_folder2, "scenarios/")

## Adaptive sampling
for(as_runs in 1:10) {
    cat("Adaptive sampling run ", as_runs)
    # Remove log file before resampling
    #system(paste("rm", log_file))
    
    # Obtain new samples and create the parameter table with new data points
    new_points_obj = select_points_uniform(GP_updated, param_ranges, n_samples)
    param_tab = create_param_table_from_samples(new_points_obj$param_tab, cv_result$seasonality, 
                                                cv_result$biting_pattern,cv_result$Decay_Scen,cv_result$maxGroup,cv_result$Access, n_seeds) 
    
    param_tab$Drug_halflife = 5
    param_tab$Drug_efficacy= 0.95
    table_file1 = paste(as_folder1, "param_tab.txt", sep="")
    table_file2 = paste(as_folder2, "param_tab.txt", sep="")
    
    write.table(param_tab, table_file1, sep = "\t", quote = FALSE, col.names = TRUE,
                row.names = FALSE)
    write.table(param_tab, table_file2, sep = "\t", quote = FALSE, col.names = TRUE,
                row.names = FALSE)
    variances = rbind(variances, new_points_obj$variances_vec)
    
    # Create the new scenarios and run OpenMalaria simulations
    file.copy(scaffold_file1, as_folder1,overwrite=TRUE)
    file.copy(scaffold_file2, as_folder2,overwrite=TRUE)
    
    # To make the workflow more generic, the OM workflow folder could be passed as an argument to this script
    system(paste("cd /scicore/home/smith/burlyd00/smc_lai/analysis_workflow/1_OM_basic_workflow/; bash OM_base_workflow.sh", as_folder1,"5"))
    
    system(paste("cd /scicore/home/smith/burlyd00/smc_lai/analysis_workflow/1_OM_basic_workflow/; bash OM_base_workflow.sh", as_folder2,"5"))
    
    
    # Run postprocessing, get the new training data and update GP
    processing_results1 = postprocess_OM_as(paste0(as_folder1,"om/"), param_tab, follow_up)
    processing_results2 = postprocess_OM_as(paste0(as_folder2,"om/"), param_tab, follow_up)
    
    save(processing_results1, file= paste0(as_folder1,"postProc.RData"))
    
    save(processing_results2, file= paste0(as_folder2,"postProc.RData"))
    
    processing_results <- aggregated_OM
    new_train_data = processing_results[, c(rownames(param_ranges), predicted)]
    n = ncol(new_train_data)
    prdata = find_reps(X = as.matrix(new_train_data[, 1:(n-1)]), Z = as.matrix(new_train_data[, n]),
                       rescale = FALSE, normalize = FALSE)
    new_GP_updated = update(GP_updated, prdata$X0, prdata$Z0, maxiter = 0, maxit=0)
    
    train_data2 <- data.frame(rbind(train_data[, c(input_parameters, predicted)],new_train_data))
    GP_model2 = train_GP_matern(train_data2)
    
    test_data2 =   test_data %>% group_by(Scenario_Name) %>% summarise_all(funs(if(is.numeric(.)) median(., na.rm = TRUE) else unique(.)))
    test_GP_plot(GP_model=GP_model2, test=test_data2[,c(rownames(param_ranges), predicted)],save=paste(results_folder, exp_name, "_R2.jpg", sep=""))
    
    test_GP_plot(GP_model=new_GP_updated, test=test_data2[,c(rownames(param_ranges), predicted)],save=paste(results_folder, exp_name, "_R2.jpg", sep=""))
    
    GP_updated = new_GP_updated 
    rm(new_GP_updated)
    # Cleanup of simulations necessary before the next run of OpenMalaria
    system(paste("rm -r", om_results_dir))
    system(paste("rm -r", om_scenarios_dir))
}

# Save the new GP model object to a file
as_result = list(GP_model = GP_updated, seasonality = cv_result$seasonality, 
             biting_pattern = cv_result$biting_pattern)
exp_name = tools::file_path_sans_ext(basename(gp_file))
as_file = paste(dirname(as_folder), "/", exp_name, "_as.RData", sep="")
save(as_result, file = as_file)

# Save the vector of variances
var_file = paste0(as_folder, "variances.RData")
save(variances, file = var_file)

