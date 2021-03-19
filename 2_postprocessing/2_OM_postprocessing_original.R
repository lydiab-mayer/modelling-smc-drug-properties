##############################
# Main script for specifying parameter values and running OpenMalaria simulations on the cluster. 
# 
#
# created 12.02.2021
#lydia.burgert@unibas.ch 
#############################

# Setup
rm(list = ls())
library(tgp)
set.seed(42)

# User 
user = strsplit(getwd(), "/", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]][5]

# Working directory
setwd(paste0("/scicore/home/penny/",user,"/M3TPP"))

# Source function scripts
source(paste0("./analysisworkflow/2_postprocessing/genOMpostprocscripts.R"))

# insert experiment name here
exp ="..."

# are you post-processing an implementation or trial setting?

setting= "implementation"
followup = 5
trialweeks=NA


#setting="trial"
#followup= 5
#trialweeks= c(42, 66)

genOMpostprocscripts(exp, setting, followup,trialweeks)


