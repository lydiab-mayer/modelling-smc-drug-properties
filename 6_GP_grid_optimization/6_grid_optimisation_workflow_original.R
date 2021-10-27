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
set.seed(42)

# User 
user = strsplit(getwd(), "/", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]][5]

# Working directory
setwd(paste0("/scicore/home/penny/", user, "/M3TPP"))

# Source function scripts
source(paste0("./analysisworkflow/6_GP_grid_optimization/gengridoptimizationscripts.R"))

# insert experiment name here
exp = "..."

# specify the predicted outcome, the desired reduction and parameter to be optimized  

# same name as in the gp trained folder 
pred_list = c("inc_red_int_Avg", "sev_red_int_Avg", "mor_red_int_Avg")


# Set grid of parameter values to search
ngrid <- c(NA, NA, NA)
# ngrid = c("Coverage" = ..., "Halflife" = ..., "Efficacy" = ...)

# Set target range size (1 is per 1% jumps, 10 by 10% etc.)
target_range_size = 10

# specify whether, when training the GP, input parameters were scaled to c(0, 1) (TRUE) or used on their original scale (FALSE)
scale = TRUE

# loop run for each
for(i in pred_list){
  pred = i
  print(i)
  # Run
  genoptimizationscripts(exp, pred, ngrid, target_range_size, scale)
}


