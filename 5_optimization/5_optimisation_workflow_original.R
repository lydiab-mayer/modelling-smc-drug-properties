##############################
# Main script for running GP training functions on post-processed OM simulations 
# 
#
# created 12.02.2021
#lydia.burgert@unibas.ch 
#############################
rm(list = ls())
setwd("~/M3TPP")
set.seed(42)

source("./analysisworkflow/5_optimization/genoptimizationscripts.R")

library(tgp)
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

# insert experiment name here
exp ="..."

# specify the predicted outcome, the desired reduction and parameter to be optimized  

# same name as in the gp trained folder 
predicted = "prevred_int_y10"
reductions <- c(50)

# needs to be a predictor in the gp -> continuous parameter sampled when simulating
optimized="Efficacy"



genoptimizationscripts(exp, predicted, reductions,optimized)



