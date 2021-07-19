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
source(paste0("./analysisworkflow/4_sensitivity_analysis/gensensanalysisscripts.R"))

# insert experiment name here
exp ="..."

# specify predicted parameter
pred_list = c("prev_red_all","prev_red_210","prev_red_int","inc_red_05","inc_red_int","inc_red_all","inc_red_int_5mo")

# loop run for each
for(i in pred_list){
  predicted = i
  print(i)
  gensensanalysisscripts(exp, predicted)
}






