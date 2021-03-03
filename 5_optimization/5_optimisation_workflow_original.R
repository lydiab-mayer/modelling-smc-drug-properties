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
exp ="E0_MAB"

# specify the predicted outcome and the desired reduction 

predicted = "prevred_int_y10"
reductions <- c(10,30)




genoptimizationscripts(exp, predicted, reductions)



