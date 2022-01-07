############################################################
# iTPP3_tradeoffs_cutoff_table
#
# Generates table of cut-off criteria for all experiments
# Note: This script depends on outputs of the script 5_optimisation_workflow.R
#
# Written by Lydia Braunack-Mayer
############################################################

# ----------------------------------------------------------
# Setup
# ----------------------------------------------------------

rm(list = ls())

# !!! Insert your experiment name here as a string, e.g. "MyExperiment" !!!
exp_list <-c("iTPP3_tradeoffs_3rounds", "iTPP3_tradeoffs_4rounds", "iTPP3_tradeoffs_5rounds") 

# !!! Insert your predicted parameter here. Note that this must match with one column name in post-processing files !!!
pred_list <- c("inc_red_int_Tot", "sev_red_int_Tot", "mor_red_int_Tot")

library(tidyr)
library(dplyr)
library(hetGP)

user <- strsplit(getwd(), "/", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]][5]
GROUP_dr <- "/scicore/home/penny/GROUP/M3TPP/"
setwd(paste0("/scicore/home/penny/", user, "/M3TPP/"))

out <- matrix(nrow = 1, ncol = 6)
out[1, ] <- c("Scenario", "Target Impact", "Program Coverage", "Round Coverage", "Duration of Protection", "Initial Efficacy")


# ----------------------------------------------------------
# Calculate cutoff criteria
# ----------------------------------------------------------

for(exp in exp_list) {
  
  # Set up
  if (!dir.exists(paste0("./Experiments/", exp, "/Outputs"))) dir.create(paste0("./Experiments/", exp, "/Outputs"))
  
  load(paste0(GROUP_dr, exp, "/param_ranges.RData"))
  param_ranges_cont
  
  SPAQ <- as.matrix(c("Coverage1" = 1, "Coverage2" = 0.95, "Halflife" = 31, "Efficacy" = 1))
  SPAQ <- (SPAQ - param_ranges_cont[, 1]) / (param_ranges_cont[, 2] - param_ranges_cont[, 1])
  
  
  # Import settings
  setting <- c()
  
  for(pred in pred_list) {
    setting <- c(setting, Sys.glob(paste0(GROUP_dr, exp, "/gp/GP_grid_optimization/", pred, "/*")))
  }
  
  setting_id <- sub(".rds", "", sub("opt_", "", basename(setting)))
  
  # Calculate cutoff criteria
  for(i in 1:length(setting_id)){
    print(paste0("Generating table for setting ", i , " of ", length(setting_id), " for experiment ", exp))
  
    # Set up table to store results
    out_int <- matrix(nrow = 2, ncol = 6)
    out_int[1, ] <- c("Scenario", "Target Impact", "Program Coverage", "Round Coverage", "Duration of Protection", "Initial Efficacy")
    out_int[2, 1] <- setting_id[i]
    
    # Read in optimisation result
    df <- readRDS(setting[i])$scenario
    
    # Read in associated GP model
    gp_path <- sub(".rds", "_cv.RData", sub("opt", "seeds", sub("GP_grid_optimization", "trained", setting[i])))
    load(gp_path) # loads file called cv_result
    
    # Identify predicted outcome for SPAQ
    target <- floor(predict(x = t(SPAQ), object = cv_result$GP_model)$mean)
    out_int[2, 2] <- paste0(target, "%")
    
      # Calculate cutoff criteria to meet the predicted impact of SPAQ
    res <- df[df$mean >= target, ] %>%
       summarise(Coverage1 = min(Coverage1),
                 Coverage2 = min(Coverage2),
                 Halflife = min(Halflife),
                 Efficacy = min(Efficacy)) %>%
       as.numeric()
        
    res[c(1, 2, 4)] <- paste0(res[c(1, 2, 4)]*100, "%")
    res[3] <- paste0(res[3], " days")
        
    out_int[2, 3:6] <- res
    
    # Store outputs
    out <- rbind(out, out_int[2, ])
  }
}

# Save resulting matrix against future use
saveRDS(out, paste0("./Experiments/", exp, "/Outputs/tradeoffs_cutoff_criteria.rds"))
#out <- readRDS(paste0("./Experiments/", exp, "/Outputs/tradeoffs_cutoff_criteria.rds"))


# ----------------------------------------------------------
# Format cutoff criteria outputs
# ----------------------------------------------------------

# Add column names
colnames(out) <- out[1, ]
out <- as.data.frame(out[-1, ])

# Add additional columns for each scenario factor
out$Scenario <- sub("_red_int_Tot", "", sub("iTPP3tradeoffs", "", out$Scenario))
out <- out %>%
  separate(col = Scenario,
           into = c("Rounds", "Seasonality", "System", "EIR", "Agegroup", "Decay", "Access", "Timing", "Endpoint"),
           sep= "_")

# Merge in baseline prevalence
prev <- read.csv(paste0("./Experiments/iTPP3_tradeoffs_4rounds/Outputs/Prevalence_prior_to_intervention.csv"))
prev <- prev[, c("Seasonality", "EIR", "Access", "annual_prev")]
prev$annual_prev <- paste0(round(prev$annual_prev*100, 0), "%")
out <- merge(out, prev, by = c("Seasonality", "EIR", "Access"))


# Format text
out$Endpoint <- factor(out$Endpoint, levels = c("inc", "sev", "mor"))
out$Endpoint <- recode(out$Endpoint,
                      "inc" = "Clinical incidence reduction",
                      "sev" = "Severe disease reduction",
                      "mor" = "Mortality reduction")

out$Rounds <- factor(out$Rounds, levels = c("3rounds", "4rounds", "5rounds"))
out$Rounds <- recode(out$Rounds,
                     "3rounds" = "SMC deployed 3 times a year",
                     "4rounds" = "SMC deployed 4 times a year",
                     "5rounds" = "SMC deployed 5 times a year")

out$Seasonality <- factor(out$Seasonality, levels = c("sharpseasonal", "wideseasonal"))
out$Seasonality <- recode(out$Seasonality,
                          "sharpseasonal" = "short seasonal profile",
                          "wideseasonal" = "long seasonal profile")

out$Agegroup <- factor(out$Agegroup, levels = c("5", "10"))
out$Agegroup <- recode(out$Agegroup,
                       "5" = "children aged between 3 and 59 months",
                       "10" = "children aged between 3 and 119 months")

out$Access <- factor(out$Access, levels = c("0.04", "0.24"))
out$Access <- recode(out$Access,
                     "0.04" = "low access to care",
                     "0.24" = "high access to care")

out$Scenario <- paste0(out$Rounds, " to ", out$Agegroup, ", in a setting with a ", out$Seasonality, " and ", out$Access)

head(out)


# Tidy up table and order rows

names(out)[names(out) == "annual_prev"] <- "PfPR2_10"

out$indexOutcome <- sub(" \\(.*", "", out[, "Target Impact"])
out$indexOutcome <- factor(out$indexOutcome, levels = c("High", "Moderate", "Low"))

out <- out[order(out$Seasonality, out$Access, out$Rounds, out$Agegroup, out$PfPR2_10, out$Endpoint, out$indexOutcome), ]

out <- out[, c("Scenario", "PfPR2_10", "Endpoint", "Target Impact", "Program Coverage", "Round Coverage", "Duration of Protection", "Initial Efficacy")]


# Remove duplicated text
out <- out %>%
  mutate(indexScenario = Scenario == lag(Scenario, 1),
         indexPfPR = PfPR2_10 == lag(PfPR2_10, 1))

out[is.na(out)] <- FALSE

out[out$indexScenario, ]$Scenario <- NA
out[out$indexPfPR, ]$PfPR2_10 <- NA

out <- out[, !(names(out) %in% c("indexScenario", "indexPfPR", "indexOutcome"))]

head(out)


# ----------------------------------------------------------
# Write results
# ----------------------------------------------------------

write.csv(out, 
          file = paste0("./Experiments/", "iTPP3_tradeoffs_4rounds", "/Outputs/tradeoffs_cutoff_criteria.csv"),
          na = "",
          row.names = FALSE)


