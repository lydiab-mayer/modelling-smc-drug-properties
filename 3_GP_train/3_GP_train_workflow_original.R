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
source(paste0("./analysisworkflow/3_GP_train/genGPtrainscripts.R"))

# insert experiment name here
exp ="..."

# specify predicted parameter
pred_list = c("prev_red_all","prev_red_210","prev_red_int","inc_red_05","inc_red_int","inc_red_all","inc_red_int_5mo")

# specify lower and upper bounds for the vector of length scale parameters used to perform maximum likelihood estimation for the GP in the format "#/#/#"
# For example, lower = "0/0/0" and upper = "1/20/2 sets c(0, 1) as bounds for length scale parameter 1, c(0, 20) for length scale parameter 2 and c(0, 2) for length scale parameter 3
lower = "0.001/0.001/0.001" # default values chosen as close to 0 (exact value of 0 leads to non-convergence)
upper = "10/10/10" # default values chosen suitable for most problems, can be adjusted in case of non-convergence

# specify whether input parameters should be scaled to c(0, 1) (TRUE) or used on their original scale (FALSE)
scale = TRUE

# loop run for each
for(i in pred_list){
  predicted = i
  print(i)
  genGPtrainscripts(exp, predicted, lower, upper, scale)
}



