rm(list = ls())

setwd("~/smc_lai/analysis_workflow/analysis_scripts")

library(hetGP)
library(ggplot2)
library(ggpubr)
library(reshape2)
library(sensitivity)
library(tgp)
library(gtools)
library(iml)


EIR <- 5.0

#use existing parameter table and small training set for testing
load('/scicore/home/smith/GROUP/smc_lai/E3_E4_comparison/gp_5/trained/rel_diff_pppy_y10_all/seeds_E3_E4_comparison_seasonal_Low_indoor_cv.RData')
training_data <- cv_result$train_data[1:100,]

param_col <- c("Coverage","Halflife","Efficacy")

#train GP 1: upper limit
response_col_upper_limit <- c("UL")
prdata = find_reps(X = as.matrix(training_data[, param_col]),
                   Z = as.matrix(training_data[, response_col_upper_limit]),
                   rescale = FALSE, normalize = FALSE)

GP_model_upper_limit = mleHetGP(X = list(X0 = as.matrix(prdata$X0),
                                         Z0 = as.matrix(prdata$Z0), 
                                         mult = prdata$mult),
                                Z = prdata$Z, 
                                lower = rep(0.0001, length(param_col)), 
                                upper = rep(10, length(param_col) ),
                                covtype = "Gaussian")


#train GP 2: cutoff
response_col_cutoff <- c("UL")
prdata = find_reps(X = as.matrix(training_data[, param_col]),
                   Z = as.matrix(training_data[, response_col_cutoff]),
                   rescale = FALSE, normalize = FALSE)

GP_model_cutoff = mleHetGP(X = list(X0 = as.matrix(prdata$X0),
                                    Z0 = as.matrix(prdata$Z0), 
                                    mult = prdata$mult),
                           Z = prdata$Z, 
                           lower = rep(0.0001, length(param_col)), 
                           upper = rep(10, length(param_col) ),
                           covtype = "Gaussian")



# shapley effects

calc_shapley_effects = function(GP_model_upper_limit, GP_model_cutoff, param_spec, num_points){
  S_mat = T_mat = NULL
  print(param_spec)
  # define wrapper for the GP_model prediction to be called by shapley function
  GP_f = function(X){
    upper_limit <- predict(x = as.matrix(X), object = GP_model_upper_limit)
    cutoff <- predict(x = as.matrix(X), object = GP_model_cutoff)
    out <- as.numeric(upper_limit$mean < cutoff$mean)
    
    return(out)
  }
  
  d <- nrow(param_spec)
  
  # Xall: a function to generate a n-sample of a d-dimensional input vector 
  # (following the required joint distribution). Here: joint distribution = uniform?
  Xall <- function(n) {
    a <- runif(n, param_spec[1,1], param_spec[1,2])
    b <- runif(n, param_spec[2,1], param_spec[2,2])
    c <- runif(n, param_spec[3,1], param_spec[3,2])
    return(matrix(cbind(a,b,c), nrow = n))
  }
  
  # Xset: Xset(n, Sj, Sjc, xjc) is a function to generate a n-sample of a d-dimensional 
  # input vector corresponding to the indices in Sj conditional on the input values xjc with 
  # the index set Sjc (following the required joint distribution). Here: joint distribution = uniform?
  
  Xset <- function(n, Sj, Sjc, xjc){
    a <- runif(n, param_spec[1,1], param_spec[1,2])
    b <- runif(n, param_spec[2,1], param_spec[2,2])
    c <- runif(n, param_spec[3,1], param_spec[3,2])
    m <- matrix(cbind(a,b,c), nrow = n)
    return(m[,Sj])
  }
  
  shapley <- shapleyPermEx(model = GP_f, Xall, Xset,d=d, Nv=1e4, No = 1e3, Ni = 3)
  
  return(shapley)
}



# shapley indices for non-inferiority
param_ranges <- matrix(nrow = 3, ncol = 2)
rownames(param_ranges) <- c("Coverage","Halflife","Efficacy")
param_ranges["Coverage", ] <- c(0,1)
param_ranges["Halflife", ] <- c(50,130)
param_ranges["Efficacy", ] <- c(0,1)


a <- calc_shapley_effects(GP_model_upper_limit, GP_model_cutoff, param_ranges)


