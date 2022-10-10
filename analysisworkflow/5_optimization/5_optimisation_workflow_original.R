##############################
# Main script for running GP training functions on post-processed OM simulations 
# 
#
# created 12.02.2021
#lydia.burgert@unibas.ch 
#
# updated September 2021
# lydia.braunack-mayer@swisstph.ch
#############################

# Setup
rm(list = ls())
# library(tgp)
# library(hetGP)
# library(ggplot2)
# library(viridis)
# library(sensitivity)
# library(multisensi)
# library(lhs)
# library(dplyr)
# library(reshape2)
# library(gridExtra)
set.seed(42)

# User 
user = strsplit(getwd(), "/", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]][5]

# Working directory
setwd(paste0("/scicore/home/penny/", user, "/M3TPP"))

# Source function scripts
source(paste0("./analysisworkflow/5_optimization/genoptimizationscripts.R"))

# insert experiment name here
exp ="..."

# specify the predicted outcome, the desired reduction and parameter to be optimized  

# same name as in the gp trained folder 
predicted = "..."
reductions_list = seq(10, 80, 10)

# needs to be a predictor in the gp -> continuous parameter sampled when simulating
optimized = "Halflife"

# Grid points (number of breaks of continuous variables)
n_gridpoints = 10

# specify whether, when training the GP, input parameters were scaled to c(0, 1) (TRUE) or used on their original scale (FALSE)
scale = TRUE

# loop run for each
for(i in reductions_list){
  reductions = i
  print(i)
  # Run
  genoptimizationscripts(exp, predicted, reductions, optimized, n_gridpoints, scale)
}


