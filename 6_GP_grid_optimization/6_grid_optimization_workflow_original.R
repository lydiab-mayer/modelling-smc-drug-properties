##################################################
##################################################
###                                            ###
### STEP 6: GP-BASED GRID SEARCH OPTIMIZATION  ###
###                                            ###
##################################################
##################################################


### -------------------------------------------------------------------------
###
### M3TPP PROJECT:
### Additional script for running optimization procedure
### of key performance characteristics from a grid search method
### using a pre-train GP emulator
### 
### Original script:
### Created 29.10.2021
### lydia.braunack-mayer@swisstph.ch
### 
### R version 3.6.0
###
### -------------------------------------------------------------------------


# Setup
rm(list = ls())
set.seed(42)

# User 
user <- strsplit(getwd(), "/", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]][5]

# Working directory
setwd(paste0("/scicore/home/penny/", user, "/M3TPP"))

# Source function scripts
source(paste0("./analysisworkflow/6_GP_grid_optimization/gengridoptimizationscripts.R"))

# insert experiment name here
exp <- "iTPP3_tradeoffs"

# same name as in the gp trained folder 
pred_list <- c("inc_red_int_Avg")

# specify whether, when training the GP, input parameters were scaled to c(0, 1) (TRUE) or used on their original scale (FALSE)
scale <- TRUE

# Set grid of parameter values to search
ngrid <- "10/10/10/10"
# ngrid = c("Coverage" = ..., "Halflife" = ..., "Efficacy" = ...)

# Set target range size (1 is per 1% jumps, 10 by 10% etc.)
target_range_size <- 10

# loop run for each
for(i in pred_list){
  pred <- i
  print(i)
  # Run
  gengridoptimizationscripts(exp, pred, scale, ngrid, target_range_size)
}


