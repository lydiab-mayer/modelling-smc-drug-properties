##############################
# Main script for running GP training functions on post-processed OM simulations 
# 
#
# created 12.02.2021
#lydia.burgert@unibas.ch 
#############################

# Setup
rm(list = ls())
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
set.seed(42)

# User 
user = strsplit(getwd(), "/", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]][5]

# Working directory
setwd(paste0("/scicore/home/penny/",user,"/M3TPP"))

# Source function scripts
source(paste0("./analysisworkflow/5_optimization/genoptimizationscripts.R"))

# insert experiment name here
exp ="..."

# specify the predicted outcome, the desired reduction and parameter to be optimized  

# same name as in the gp trained folder 
predicted = "prevred_int_y10"
reductions <- c(50)

# needs to be a predictor in the gp -> continuous parameter sampled when simulating
optimized="Halflife"

# Run
genoptimizationscripts(exp, predicted, reductions,optimized)



