######################################
######################################
###                                ###
### STEP 5: OPTIMIZATION ANALYSIS  ###
###                                ###
######################################
######################################


### -------------------------------------------------------------------------
###
### M3TPP PROJECT:
### Main script for running optimization analyses
### 
### Original script:
### Created 12.02.2021
### lydia.burgert@unibas.ch 
###
### Adapted script:
### Saved 01.09.2021
### narimane.nekkab@unibas.ch
###
### R version 3.6.0
###
### -------------------------------------------------------------------------

##############
### HEADER ###
##############

# Clear environment
rm(list = ls())

# Set seed for replication
set.seed(42)

# Library (double-check which are no longer required and remove)
library(tgp)
library(hetGP)
library(ggplot2)
library(viridis)
library(sensitivity)
library(multisensi)
library(lhs)
library(dplyr)
library(reshape2)
library(gridExtra)

# User 
user = strsplit(getwd(), "/", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]][5]

# Working directory
setwd(paste0("/scicore/home/penny/",user,"/M3TPP"))

# Source function scripts
source(paste0("./analysisworkflow/5_optimization/genoptimizationscripts.R"))


##################
### EXPERIMENT ###
##################

# Insert experiment name here
exp ="..."


##################
### PARAMETERS ###
##################

# Specify the following 3 variables:
# --> 1 predicted outcome
# --> multiple desired health targets (% reduction) 
# --> 1 parameter to be optimized (Coverage, Halflife, Efficacy) --> verify spelling

###############
### OUTCOME ###

# Choose 1 outcome variable for optimization from full list (to be updated)
pred_list = c("prev_red_all","prev_red_210","prev_red_int","inc_red_05",
              "inc_red_int","inc_red_all","inc_red_int_5mo")
predicted = pred_list[1]

###############
### TARGETS ###

# Choose health targets (% reduction) for loop
targets_list = seq(10,80,10)

##########################
### OPTIMIZED VARIABLE ###

# Choose 1 variable to optimize
# Needs to be a predictor in the gp -> continuous parameter sampled when simulating
optimized="Halflife"
# optimized="Coverage"
# optimized="Efficacy"

###################
### GRID POINTS ###

# Grid points (number of breaks of continuous variables)
# Increasing number increases number of finite values & run time
n_gridpoints = 10


########################
### RUN OPTIMIZATION ###
########################

# loop run for each target
for(i in targets_list){
  targets = i
  print(i)
  # Run
  genoptimizationscripts(exp, predicted, targets, optimized, n_gridpoints)
}


