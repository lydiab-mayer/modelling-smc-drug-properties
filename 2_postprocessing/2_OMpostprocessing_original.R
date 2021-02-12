##############################
# Main script for specifying parameter values and running OpenMalaria simulations on the cluster. 
# 
#
# created 12.02.2021
#lydia.burgert@unibas.ch 
#############################
rm(list = ls())
setwd("~/M3TPP")
set.seed(42)
source("./analysisworkflow/2_postprocessing/genOMpostprocscripts.R")

library(tgp)



# insert experiment name here
exp ="..."

# are you post-processing an implementation or trial setting?

setting= "implementation"
followup = 5
trialweeks=NA


#setting="trial"
#followup= 5
#trialweeks= c(42, 66)

genOMpostprocscripts(setting, followup,trialweeks)


