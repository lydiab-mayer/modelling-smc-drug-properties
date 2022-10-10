################################
################################
###                          ###
### STEP 2: POST-PROCESSING  ###
###                          ###
################################
################################


### -------------------------------------------------------------------------
###
### M3TPP PROJECT:
### Main script for running post-processing of OM simulations 
### to aggregate data which will be used to train GP
### 
### Original script:
### Created 12.02.2021
### lydia.burgert@unibas.ch 
###
### Adapted script:
### Saved 22.10.2021
### lydia.braunack-mayer@swisstph.ch
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

# Library
library(dplyr)

# User 
user = strsplit(getwd(), "/", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]][5]

# Working directory
setwd(paste0("/scicore/home/penny/", user, "/M3TPP"))

# Source function scripts
source(paste0("./analysisworkflow/2_postprocessing/genOMpostprocscripts.R"))


##################
### EXPERIMENT ###
##################

# Insert experiment name here
exp = "..."

##################
### PARAMETERS ###
##################

# Identify the minimum intervention age as a numeric, e.g. 0.25
min_int <- ...

# Identify the start date of your monitoring period as "YYYY-MM-DD", e.g. "2030-01-01"
date <- "..."

# Identify the first month of intervention deployment as the month's 3 letter abbreviation, e.g. "Sep"
fmonth <- "..."

# Identify the number of months the intervention is deployed for as a numeric, e.g. 5
months <- ...

# Identify the year to treat as the counterfactual as a numeric, i.e. the baseline year before intervention, e.g. 2034
year_counterfactual <- ...

# Identify the year to treat as the intervention as a numeric, i.e. the year for which incidence reduction will be calculated, e.g. 2039 
year_intervention <- ...

# Write parameters to file
if (!dir.exists(paste0("/scicore/home/penny/GROUP/M3TPP/", exp, "/postprocessing/"))) dir.create(paste0("/scicore/home/penny/GROUP/M3TPP/", exp, "/postprocessing/"))
pp_param <- data.frame(min_int, date, fmonth, months, year_counterfactual, year_intervention)
saveRDS(pp_param, paste0("/scicore/home/penny/GROUP/M3TPP/", exp, "/postprocessing/pp_param_values.RDS"))

##########################
### RUN POSTPROCESSING ###
##########################

# Run
genOMpostprocscripts(exp, date, fmonth, months, year_counterfactual, year_intervention, min_int)


