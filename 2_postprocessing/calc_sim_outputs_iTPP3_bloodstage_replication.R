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
source(paste0("/scicore/home/penny/", user, "/M3TPP/analysisworkflow/2_postprocessing/postprocessing_resources_iTPP3_bloodstage_replication.R"))

# Read in command arguments
args <- commandArgs(TRUE)
dir <- args[1]
split_file <- args[2]
date <- args[3]
timesteps <- args[4]

# # Sample command arguments, retained here for testing
# dir <- "/scicore/home/penny/GROUP/M3TPP/iTPP3_bloodstage_replication/om/"
# split_file <- "/scicore/home/penny/GROUP/M3TPP/iTPP3_bloodstage_replication/postprocessing/split/iTPP3bloodstagereplication_350_1_0.020831339.txt"
# date <- "1970-01-01"
# timesteps <- "2906-2918"

cat("Command arguments:")
print(paste("dir:", dir))
print(paste("Split file:", split_file))
print(paste("Date of first monitoring:", date))
print(paste("Time steps monitored:", timesteps))

timesteps <- do.call(as.numeric, strsplit(timesteps, "-"))
timesteps <- timesteps[1]:timesteps[2]

# -------------------------------------------------------------------------------------------------------------
# PERFORM POST-PROCESSING
# -------------------------------------------------------------------------------------------------------------

# Postprocess the OpenMalaria simulations
cat("Run OM postprocessing function:")
postprocess.om(dir = dir, 
               param.file = split_file,
               date = date,
               timesteps = timesteps)
