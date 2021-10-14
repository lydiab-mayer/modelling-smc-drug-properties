# -------------------------------------------------------------------------------------------------------------
#
# Handy functions and mini scripts to accompany the M3TPP project team
# 
# Created 14.10.2021
# lydia.braunack-mayer@swisstph.ch 
#
# R version 3.6.0
#
# -------------------------------------------------------------------------------------------------------------



# TRANSLATE THE PEAK OF A FOURIER SERIES ----------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------

fourier_translation <- function(coef, t) {
  
  # Given a set of courier coefficients c(a1, b1, a2, b2) and a desired timeshift t, this helper function calculates the new
  # fourier coefficients translated by t
  # The form of the fourier equation is described here: https://github.com/SwissTPH/openmalaria/wiki/ScenarioTransmission
  
  names(coef) <- c("a1", "b1", "a2", "b2")
  
  a1_hat <- coef["a1"]*cos(2*pi*t/365) + coef["b1"]*sin(2*pi*t/365)
  b1_hat <- coef["b1"]*cos(2*pi*t/365) - coef["a1"]*sin(2*pi*t/365)
  a2_hat <- coef["a2"]*cos(4*pi*t/365) + coef["b2"]*sin(4*pi*t/365)
  b2_hat <- coef["b2"]*cos(4*pi*t/365) - coef["a2"]*sin(4*pi*t/365)
  
  return(c(a1_hat, b1_hat,  a2_hat, b2_hat))
  
}

seasonality <- c(-1.3715341885884602, 1.196649392445882, 0.15281548102696738, 0.16501147548357645)
translation <- -92
fourier_translation(seasonality, translation)

seasonality <- c(-2.3171836535135903, 2.114041487375895,	0.384947935740153,	0.3511592149734497)
fourier_translation(seasonality, translation)

fourier <- function(coef, t) {
  
  # Given a set of courier coefficients c(a1, b1, a2, b2), this function returns the value of the fourier series at time t
  # The form of the fourier equation is described here: https://github.com/SwissTPH/openmalaria/wiki/ScenarioTransmission
  
  names(coef) <- c("a1", "b1", "a2", "b2")
  
  f_t <- exp(coef["a1"]*cos(2*pi*t/365) + coef["b1"]*sin(2*pi*t/365) + coef["a2"]*cos(4*pi*t/365) + coef["b2"]*sin(4*pi*t/365))
  f_t <- unname(f_t)
  
  return(f_t)
  
}

plot(1:365, fourier(seasonality, 1:365), type = "l")
lines(1:365, fourier(fourier_translation(seasonality, translation), 1:365))

# -------------------------------------------------------------------------------------------------------------



# TRANSLATE OPENMALARIA'S TIME STEPS TO THEIR CORRESPONDING DATE ----------------------------------------------
# -------------------------------------------------------------------------------------------------------------

time.to.date <- function(t, time.step, date) {
  # Function to convert the time step in an OpenMalaria outputs file to a date
  #
  # Inputs: 
  # t: the time steps to convert to date format, a numeric or numeric vector
  # time.step: the time step used in OpenMalaria, a numeric in days. Usually 1 or 5
  # date: the starting date of your OpenMalaria survey period in the format "%Y%m%d"
  #
  # Outputs: a vector containing the dates corresponding to each time step
  
  date <- as.Date(date, format = "%Y-%m-%d")
  out <- date + (t - 1)*time.step
  return(out)
}

time.to.date(1:73, time.step = 5, date = "2030-01-01")

# -------------------------------------------------------------------------------------------------------------



# EXTRACT AGE GROUPS FROM AN OPENMALARIA XML ------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------

extract.agegroups <- function(path) {
  # Function to extract age groups from an OpenMalaria xml
  #
  # Inputs: path, a file pathway to an OpenMalaria xml
  # Outputs: a vector containing the age groups contained within the xml
  
  require(stringr)
  
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

# -------------------------------------------------------------------------------------------------------------



# CALCULATE PREVALENCE PRIOR TO INTERVENTION ------------------------------------------------------------------
# -------------------------------------------------------------------------------------------------------------

# Note that for this script to work, your postprocessing must be complete. Postprocessing must also export
# three non-standard outcome measures: prev_all_before, prev_210_beg, prev_int_beg

### SETUP ###

# Define experiment name
exp <- "..."

# Define user
user = strsplit(getwd(), "/", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]][5]

# Define postprocessing directory
dir <- paste0("/scicore/home/penny/GROUP/M3TPP/", exp, "/postprocessing/")

# Load pacakges
require(dplyr)


### IMPORT DATA ###

# Collate results of postprocessing

files <- list.files(dir, pattern ="^[agg]")

out <- data.frame()

for (i in 1:length(files)) {
  agg_file <- read.table(paste0(dir, files[i]), header = TRUE, sep = "")
  out <- rbind(out, agg_file)
}

head(out)


### CALCUALTE OUTPUTS ###

# For each scenario, calculate mean prevalence

out <- out %>%
  group_by(Seasonality, Biting_pattern, EIR, MaxAge, Decay, Access, Timing) %>%
  summarise(prev_all_before = mean(prev_all_before),
            prev_210_beg = mean(prev_210_beg),
            prev_int_beg = mean(prev_int_beg))


### STORE OUTPUTS ###

write.csv(out, file = paste0("/scicore/home/penny/", user, "/M3TPP/Experiments/", exp, "/Outputs/Prevalence_prior_to_intervention.csv"))

# -------------------------------------------------------------------------------------------------------------