################################
### STEP 2: POST-PROCESSING  ###
################################

# -------------------------------------------------------------------------------------------------------------
#
# Helper functions for running post-processing of OM simulations to aggregate data which will be used to train GP
# 
# Original script:
# Created 25.3.2021
# lydia.braunack-mayer@swisstph.ch 
#
# R version 3.6.0
#
# -------------------------------------------------------------------------------------------------------------


# -------------------------------------------------------------------------------------------------------------
# SET UP
# -------------------------------------------------------------------------------------------------------------

# Import packages
library(dplyr)
library(tidyr)
library(stringr)

options(dplyr.summarise.inform = FALSE)


# -------------------------------------------------------------------------------------------------------------
# DEFINE HANDY HELPER FUNCTIONS REFERRED TO WITHIN KEY CALCULATIONS
# -------------------------------------------------------------------------------------------------------------

extract.agegroups <- function(path) {
  # Function to extract age groups from an OpenMalaria xml
  #
  # Inputs: path, a file pathway to an OpenMalaria xml
  # Outputs: a vector containing the age groups contained within the xml
  
  #require(stringr)
  
  file <- readLines(path)
  
  oline <- grep("<ageGroup lowerbound=\"0\">", file)[2]
  cline <- grep("</ageGroup>", file)[2]
  
  groups <- str_extract_all(file[(oline + 1):(cline - 1)], "[0-9.]+")
  
  out <- c(0)
  
  for(i in 1:length(groups)) {
    out <- c(out, as.numeric(groups[[i]]))
  }
  
  return(out)
}


time.to.date <- function(t, time.step, date) {
  # Function to convert the time step in an OpenMalaria outputs file to a date
  #
  # Inputs: 
  # t: the time steps to convert to date format, a numeric or numeric vector
  # time.step: the time step used in OpenMalaria, a numeric in days. Usually 1 or 5
  # date: the starting date of your OpenMalaria survey period in the format "%Y%m%d"
  #
  # Outputs: a vector containing the dates corresponding to each time step
  
  # Format inputs
  date <- as.Date(date, format = "%Y-%m-%d")
  
  # Create reference sequence of dates assuming no leap years
  ref <- seq(from = date,
             to = date + max(t)*time.step + ceiling(max(t)*time.step/365/4),
             by = "day")
  ref <- ref[strftime(ref, "%m-%d") != "02-29"]
  
  # Return sequence of dates incremented by time.step
  out <- ref[t*time.step]

  # Return outputs
  return(out)
  
}
#time.to.date(t = 1:1095, time.step = 5, date = "2030-01-01")


# -------------------------------------------------------------------------------------------------------------
# DEFINE HELPER FUNCTIONS
# -------------------------------------------------------------------------------------------------------------

# # Sample inputs, retained here for testing
# om.result = read.table("/scicore/home/penny/GROUP/M3TPP/iTPP3_bloodstage_replication/om/iTPP3_bloodstage_replication_1_1_out.txt", header = FALSE)
# measure = 14
# age.group = 2:3
# time.step = 5
# date = "1970-01-01"
# 
# calculate.cpp.outcome(om.result, measure, age.group, time.step, date)

calculate.cpp.outcome <- function(om.result, measure, age.group, time.step, date){
  
  # Function to calculate cases per person for a given age group
  #
  # Inputs: 
  # om.result: the outputs of a single OpenMalaria simulation
  # measure: the measure of infections used to count cases per person
  # age.group: an integer or integer vector containing the age groups of interest
  # time.step: the time step used in OpenMalaria, a numeric in days. Usually 1 or 5
  # date: the starting date of your OpenMalaria survey period in the format "%Y%m%d"
  #
  # Outputs: data frame containing cases per person at each timestep
  
  # Load required packages
  #require(dplyr)
  #require(tidyr)

  # Format OpenMalaria output file
  colnames(om.result) <- c("time", "age_group", "measure", "value")
  
  # Error monitoring
  if (!(measure %in% om.result$measure)) {
    print(paste0("Data for measure ", measure, " could not be found. Check that this measure has been included in the SurveyOptions section of your xml."))
  }
  
  # Translate from OpenMalaria 5-day time steps to years and months
  om.result$date <- time.to.date(om.result$time, time.step = time.step, date = date)
  om.result$month <- format(om.result$date, "%b")
  om.result$year <- as.numeric(format(om.result$date, "%Y"))
  
  # Remove first time step from OpenMalaria outputs
  om.result <- om.result[om.result$time != 1, ]
  
  # Remove values for age groups other than those specified
  om.result <- om.result[om.result$age_group %in% age.group, ]
  
  # Remove measures other than that population size and the specified outcome measure
  om.result <- om.result[om.result$measure %in% c(0, measure), ]
  
  # Summarise all measures by summing up over age groups
  om.result <- om.result[, -which(names(om.result) %in% c("age_group"))] %>%
    group_by(measure, time, date, month, year) %>%
    summarise(value = sum(value))

  # Transform to wide format
  om.result <- pivot_wider(om.result, 
                           id_cols = c(time, date, month, year),
                           names_from = measure,
                           values_from = value,
                           names_prefix = "measure")
  om.result <- as.data.frame(om.result)
  
  # Calculate cases per person in target population
  om.result$cpp <- om.result[, paste0("measure", measure)] / om.result$measure0
  
  return(om.result)
  
}


# # Sample arguments, retained here for testing
# om.outcome <- calculate.cpp.outcome(om.result = read.table("/scicore/home/penny/GROUP/M3TPP/iTPP3_bloodstage_replication/om/iTPP3_bloodstage_replication_18_1_out.txt", header = FALSE),
#                                     measure = 3,
#                                     age.group = 2:3,
#                                     time.step = 5,
#                                     date = "1970-01-01")
# id <- "cpp_"
# timesteps <- 2906:2918 # 15 Oct 2009 to 15 Dec 2009
# 
# format.cpp(om.outcome, id, timesteps)

format.cpp <- function(om.outcome, id, timesteps) {

  # Function to calculate cumulative cases per person over given time steps
  #
  # Inputs: 
  # om.outcome: outputs of the function calculate.cpp.outcome
  # id: a string that will be used to name columns of the function outputs, e.g. "cpp_"
  # timesteps: vector of time steps to evaluate the outcome at
  #
  # Outputs: data frame containing a single row with the protective efficacy per time step
  
  #require(tidyr)
  
  # Set up data
  om.outcome <- om.outcome[om.outcome$time %in% timesteps, ]

  # Calculate cumulative outcomes
  om.outcome$cpp <- cumsum(om.outcome$cpp)
  
  # Format function outputs
  om.outcome <- om.outcome[, c("time", "cpp")]
  om.outcome <- pivot_wider(om.outcome,
                            id_cols = time,
                            names_from = time,
                            values_from = cpp,
                            names_prefix = id)
  
  # Return outputs
  return(om.outcome)
  
}

# # Sample arguments, retained here for testing
# dir <- "/scicore/home/penny/GROUP/M3TPP/iTPP3_bloodstage_replication/"
# param.file <- "/scicore/home/penny/GROUP/M3TPP/iTPP3_bloodstage_replication/postprocessing/split/iTPP3bloodstagereplication_350_1_0.020831339.txt"
# param.table <- read.table(param.file, sep = "\t", as.is = TRUE, header = TRUE, stringsAsFactors = FALSE)
# scenario.params <- param.table[1, ]
# om.file <- paste(dir, "om/", param.table[1, ]$Scenario_Name, "_", param.table[1, ]$SEED, "_out.txt", sep = "")
# om.result <- read.table(om.file, sep = "\t")
# date <- "1970-01-01"
# timesteps <- 2906:2918
# 
# report.results(dir, om.result, date, timesteps, scenario.params)

report.results <- function(dir, om.result, date, timesteps, scenario.params) {
  
  # Define age groups 
  #age.groups <- extract.agegroups(paste0(dir, "scaffold.xml")) # All age groups
  age.int <- 2:3
  #age.210 <- seq(which(age.groups == 2), which(age.groups == 10) - 1) # Children 2 to 10 years old
  #age.05 <- seq(which(age.groups == 0), which(age.groups == 5) - 1) # Children 0 to 5 years old

  # Calculate cases per person in intervention group
  om.outcome <- calculate.cpp.outcome(om.result = om.result, measure = 14, age.group = age.int, time.step = 5, date = date)
  cpp <- format.cpp(om.outcome = om.outcome, id = "cpp_", timesteps = timesteps)
  rm(om.outcome)
  
  # Return outputs
  out <- cbind.data.frame(scenario.params$Scenario_Name, scenario.params$SEED, cpp)
  colnames(out) <- c("Scenario_Name", "seed", names(cpp))
  return(out)
  
}


# -------------------------------------------------------------------------------------------------------------
# DEFINE WRAPPER FUNCTIONS TO CALCULATE OUTCOMES FOR MULTIPLE SIMULATIONS
# -------------------------------------------------------------------------------------------------------------

# # Sample arguments, retained here for testing
# dir <- "/scicore/home/penny/GROUP/M3TPP/iTPP3_bloodstage_replication/om/"
# param.file <- "/scicore/home/penny/GROUP/M3TPP/iTPP3_bloodstage_replication/postprocessing/split/iTPP3bloodstagereplication_350_1_0.020831339.txt"
# date <- "1970-01-01"
# timesteps <- 2906:2918

postprocess.om <- function(dir, param.file, date, timesteps) {
  
  results.folder <- dir
  dir <- paste0(dirname(dir), "/")
  dest.agg <- paste0(dir, "postprocessing/agg_", basename(param.file))
  dest.seed <- paste0(dir, "postprocessing/seeds_", basename(param.file))
  
  # Set up function
  param.table <- read.table(param.file, sep = "\t", as.is = TRUE, header = TRUE, stringsAsFactors = FALSE)
  om.outcome <- NULL
  
  # For each row in param.table, process the corresponding OpenMalaria simulation
  for(i in 1:nrow(param.table)) {
    
    skip <- FALSE
    
    print(i)
    om.file <- paste(dir, "om/", param.table[i, ]$Scenario_Name, "_", param.table[i, ]$SEED, "_out.txt", sep = "")
    
    tryCatch(if(file.exists(om.file) & file.info(om.file)$size > 0) {
      
      # Read in file
      om.result <- read.table(om.file, sep = "\t")
      
      out <- report.results(dir = dir, 
                            om.result = om.result,
                            date = date,
                            timesteps = timesteps, 
                            scenario.params = param.table[i, ])

      om.outcome <- data.frame(rbind(om.outcome, out), stringsAsFactors = FALSE)}, 
      error = function(e) {skip <<- TRUE})
    
    if (skip) {next}
  }
  
  # Summarize results over each seed
  om.agg <- om.outcome %>% 
    group_by(Scenario_Name) %>% 
    summarise_at(c(names(om.outcome)[(which(names(om.outcome) == "seed") + 1):length(names(om.outcome))]), median, na.rm = TRUE)
  
  # Prepare results
  no.seed.tab <- unique(param.table[, -c(which(colnames(param.table) == "SEED"))])
  seed.tab <- merge(no.seed.tab, om.outcome, by = c("Scenario_Name"))
  tab <- merge(no.seed.tab, om.agg, by = c("Scenario_Name"))
  
  # Write result tables (summarized and with seeds) to files
  write.table(tab, dest.agg, sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
  write.table(seed.tab, dest.seed, sep = "\t", col.names = TRUE, row.names = FALSE, quote = FALSE)
  
}
