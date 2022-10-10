############################################################
# iTPP3_tradeoffs_SPAQ
#
# Given trained GP emulators, this script identifies the likely impact of SPAQ
#
# Written by Lydia Braunack-Mayer
############################################################

# ----------------------------------------------------------
# Setup
# ----------------------------------------------------------

rm(list = ls())

# !!! Insert your experiment name here as a string, e.g. "MyExperiment" !!!
exp_list <- c("iTPP3_tradeoffs_3rounds", "iTPP3_tradeoffs_4rounds", "iTPP3_tradeoffs_5rounds")

# !!! Insert your predicted parameter here. Note that this must match with one column name in post-processing files !!!
pred_list <- c("inc_red_int_Tot", "sev_red_int_Tot", "mor_red_int_Tot")

library(tidyr)
library(dplyr)
library(hetGP)

user <- strsplit(getwd(), "/", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]][5]
GROUP_dr <- "/scicore/home/penny/GROUP/M3TPP/"
setwd(paste0("/scicore/home/penny/", user, "/M3TPP/"))


# ----------------------------------------------------------
# Generate predictions for SPAQ
# ----------------------------------------------------------

for(exp in exp_list) {
  
  # Set up
  if (!dir.exists(paste0("./Experiments/", exp, "/Outputs"))) dir.create(paste0("./Experiments/", exp, "/Outputs"))
  
  load(paste0(GROUP_dr, exp, "/param_ranges.RData"))
  param_ranges_cont
  
  
  # Import settings
  
  setting <- c()
  
  for(pred in pred_list) {
    temp <- Sys.glob(paste0(GROUP_dr, exp, "/gp/trained/", pred, "/*"))
    temp <- temp[grepl(".RData", temp, fixed = TRUE)]
    setting <- c(setting, temp)
  }
  
  (setting_id <- sub("_red_int_Tot_cv.RData", "", sub("seeds_", "", basename(setting))))
  setting_id <- sub("agg_", "", setting_id)
  
  
  # Generate database of predictions
  
  SPAQ <- as.matrix(c("Coverage1" = 1, "Coverage2" = 0.95, "Halflife" = 31, "Efficacy" = 1))
  SPAQ <- (SPAQ - param_ranges_cont[, 1]) / (param_ranges_cont[, 2] - param_ranges_cont[, 1])
  
  df <- c()
  
  for (i in 1:length(setting)) {
   
    load(setting[i]) # loads a list called cv_result
    
    SPAQ_pred <- predict(x = t(SPAQ), object = cv_result$GP_model)
    
    df <- c(df, paste0(round(SPAQ_pred$mean, 0), "%"))
  
  }
  
  df <- data.frame("Scenario" = setting_id, "Prediction" = df)
  
  
  # Format database of predictions
  
  df <- df %>%
    separate(col = Scenario, 
             into = c("Experiment", "Seasonality", "System", "EIR", "Agegroup", "Decay", "Access", "Timing", "Outcome"),
             sep = "_",
             remove = FALSE)
  
  
  # Merge in baseline prevalence
  
  prev <- read.csv(paste0("./Experiments/iTPP3_tradeoffs_4rounds/Outputs/Prevalence_prior_to_intervention.csv"))
  prev <- prev[, c("Seasonality", "EIR", "Access", "annual_prev")]
  prev$annual_prev <- paste0(round(prev$annual_prev*100, 0), "%")
  df <- merge(df, prev, by = c("Seasonality", "EIR", "Access"))
  
  
  # Reformat resulting data frame
  
  df$Outcome <- factor(df$Outcome, levels = c("inc", "sev", "mor"))
  df$Outcome <- recode(df$Outcome,
                         "inc" = "Clinical incidence reduction",
                         "sev" = "Severe disease reduction",
                         "mor" = "Mortality reduction")
  
  df$Seasonality <- factor(df$Seasonality, levels = c("sharpseasonal", "wideseasonal"))
  df$Seasonality <- recode(df$Seasonality,
                            "sharpseasonal" = "Short seasonal profile",
                            "wideseasonal" = "Long seasonal profile")
  
  df$Agegroup <- factor(df$Agegroup, levels = c("5", "10"))
  df$Agegroup <- recode(df$Agegroup,
                         "5" = "Children aged between 3 and 59 months",
                         "10" = "Children aged between 3 and 119 months")
  
  df$Access <- factor(df$Access, levels = c("0.04", "0.24"))
  df$Access <- recode(df$Access,
                       "0.04" = "Low access to care",
                       "0.24" = "High access to care")
  
  df <- df[, !(names(df) == "Scenario")] %>%
    pivot_wider(names_from = Outcome, values_from = "Prediction")
  
  df <- df[order(df$Access, df$Experiment, df$Seasonality, df$Agegroup, as.numeric(df$EIR)), ]
  
  
  # Write to file
  write.csv(df, paste0("./Experiments/", exp, "/Outputs/iTPP3_tradeoffs_SPAQtable.csv"))
}

