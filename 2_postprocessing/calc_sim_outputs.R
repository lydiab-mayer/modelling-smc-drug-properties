################################
### STEP 2: POST-PROCESSING  ###
################################

# -------------------------------------------------------------------------------------------------------------
#
# Support script for running post-processing of OM simulations to aggregate data which will be used to train GP
# 
# Original script:
# Created 15.10.2019
# lydia.braunack-mayer@swisstph.ch 
#
# Adapted from monica.golumbeanu@unibas.ch
#
# R version 3.6.0
#
# -------------------------------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------------------------------
# SET UP
# -------------------------------------------------------------------------------------------------------------

# Define user
user <- strsplit(getwd(), "/", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]][5]

# Source helper functions
source(paste0("/scicore/home/penny/", user, "/M3TPP/analysisworkflow/2_postprocessing/postprocessing_resources.R"))

# Read in command arguments
args <- commandArgs(TRUE)
dir <- args[1]
split_file <- args[2]
date <- args[3]
fmonth <- args[4]
months <- args[5]
year_counterfactual <- args[6]
year_intervention <- args[7]
min_int <- args[8]

#Sample command arguments, retained here for testing
dir <- "/scicore/home/penny/GROUP/M3TPP/iTPP3_tradeoffs/om/"
split_file <- "/scicore/home/penny/GROUP/M3TPP/iTPP3_tradeoffs/postprocessing/split/iTPP3tradeoffs_sharpseasonal_Mali_15_10_exp_0.241193660515256_May.txt"
date <- "2030-01-01"
fmonth <- "May"
months <- 5
year_counterfactual <- 2034
year_intervention <- 2039
min_int <- 0.25

cat("Command arguments:")
print(paste("dir:", dir))
print(paste("Split file:", split_file))
print(paste("Date of first monitoring:", date))
print(paste("First month of intervention:", fmonth))
print(paste("Number of months intervention is deployed:", months))
print(paste("Year for counterfactual outcomes:", year_counterfactual))
print(paste("Year for intervention outcomes:", year_intervention))
print(paste("Minimum intervention age:", min_int))


# -------------------------------------------------------------------------------------------------------------
# PERFORM POST-PROCESSING
# -------------------------------------------------------------------------------------------------------------

# Postprocess the OpenMalaria simulations
cat("Run OM postprocessing function:")
postprocess.om(dir = dir, 
               param.file = split_file,
               date = date,
               fmonth = fmonth,
               months = months,
               year.counterfactual = year_counterfactual,
               year.intervention = year_intervention,
               min.int = min_int)
