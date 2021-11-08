#####################################
#####################################
###                               ###
### STEP 4: SENSITIVITY ANALYSIS  ###
###                               ###
#####################################
#####################################


### -------------------------------------------------------------------------
###
### M3TPP PROJECT:
### Main script for running sensitivity analyses
### 
### Original script:
### Created 12.02.2021
### lydia.burgert@unibas.ch 
###
### Adapted script:
### Saved 01.09.2021
### narimane.nekkab@unibas.ch
###
### R version 3.6.0
###
### -------------------------------------------------------------------------

##############
### HEADER ###
##############

# Clear environment
rm(list = ls())

# Set seed for replication
set.seed(42)

# User 
user = strsplit(getwd(), "/", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]][5]

# Working directory
setwd(paste0("/scicore/home/penny/",user,"/M3TPP"))

# Source function scripts
source(paste0("./analysisworkflow/4_sensitivity_analysis/gensensanalysisscripts.R"))


##################
### EXPERIMENT ###
##################

# Insert experiment name here
exp ="..."


##################
### PARAMETERS ###
##################

# Specify parameter to be predicted
# Make sure GP exists or will get error

# All current variables for loop --> may be updated in future
pred_list = c("prev_red_all","prev_red_210","prev_red_int","inc_red_05",
              "inc_red_int","inc_red_all","inc_red_int_5mo")

# Or choose 1
# pred_list = "prev_red_all"

###############
### SCALING ###
###############

# Specify whether, when training the GP, input parameters were scaled to c(0, 1) (TRUE) or used on their original scale (FALSE)
scale = TRUE

####################
### RUN ANALYSIS ###
####################

# loop run for each
for(i in pred_list){
  predicted = i
  print(i)
  gensensanalysisscripts(exp, predicted)
}






