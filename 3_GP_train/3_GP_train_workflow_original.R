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

source("./analysisworkflow/3_GP_train/genGPtrainscripts.R")

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

# specify predicted parameter
predicted = "prevred_int_y10"

# specify followup

genGPtrainscripts(exp, predicted)



