############################################################
# PS_06_PlotSensitivity
#
# Visualises relationships between emulator input and predictor variables
# Note: This script depends on outputs of the script 4_sensitivityanalysis_workflow.R
#
# Written by Lydia Braunack-Mayer
############################################################

# ----------------------------------------------------------
# Setup
# ----------------------------------------------------------

rm(list = ls())

# !!! Insert your experiment name here as a string, e.g. "MyExperiment" !!!
exp <- c("iTPP3_ChemoBlood_4rounds")

# !!! Insert your predicted parameters here. Note that this must match with one column name in post-processing files !!!
pred_list <- c("inc_red_int_Tot", "sev_red_int_Tot")

library(dplyr)
library(tidyr)

user <- strsplit(getwd(), "/", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]][5]
GROUP_dr <- "/scicore/home/penny/GROUP/M3TPP/"
setwd(paste0("/scicore/home/penny/", user, "/M3TPP/SMC_TPP/"))

if (!dir.exists(paste0("./../Experiments/", exp, "/Outputs"))) dir.create(paste0("./../Experiments/", exp, "/Outputs"))

load(paste0(GROUP_dr, exp, "/param_ranges_manual.RData"))
param_ranges_cont


# ----------------------------------------------------------
# Define data to plot
# ----------------------------------------------------------

# Import settings
setting <- Sys.glob(paste0(GROUP_dr, exp, "/gp/trained/sensitivity/*"))
setting <- setting[setting != paste0(GROUP_dr, exp, "/gp/trained/sensitivity/err")]
setting_id <- sub("_cv_sidx.RData", "", sub(".*agg_", "", sub(".*seeds_", "", setting)))


# Import total effect sizes for each setting
df <- data.frame("S_eff" = c(), "T_eff" = c(), scenario = c())

for (i in 1:length(setting)) {
  
  load(setting[i]) #loads list called sobol_idx_list
  
  sobol_idx_list <- as.data.frame(sobol_idx_list)
  
  sobol_idx_list$S_eff <- sobol_idx_list$S_eff / sum(sobol_idx_list$S_eff) # rescale so total = 1
  sobol_idx_list$T_eff <- sobol_idx_list$T_eff / sum(sobol_idx_list$T_eff) # rescale so total = 1
  
  sobol_idx_list$scenario <- setting_id[i]
  sobol_idx_list$parameter <- rownames(param_ranges_cont)
  
  df <- rbind(df, sobol_idx_list)
  
}

df <- df %>%
  separate(col = scenario, 
           into = c("Experiment", "Seasonality", "System", "EIR", "Agegroup", "Access", "Timing", "IC50", "Outcome", "temp1", "temp2", "temp3"),
           sep = "_",
           remove = FALSE)
df <- df[, !(names(df) %in% c("temp1", "temp2", "temp3"))]


# Import aggregated impact for each setting

setting_id <- unique(sub("_Apr_0.02083134.*", "_Apr_0.02083134", sub("_May_0.02083134.*", "_May_0.02083134", setting_id)))
df_impact <- data.frame()

for (i in 1:length(setting_id)) {
  
  temp <- read.table(paste0(GROUP_dr, exp, "/postprocessing/agg_",setting_id[i], ".txt"), header = TRUE, sep = "")
  temp <- temp %>%
    group_by(Seasonality, Biting_pattern, EIR, MaxAge, Access, Timing) %>%
    summarise(inc = median(inc_red_int_Tot),
              sev = median(sev_red_int_Tot),
              prev = median(prev_red_int_Aug),
              mor = median(mor_red_int_Tot))
  temp$scenario <- setting_id[i]
  
  df_impact <- rbind(df_impact, as.data.frame(temp))
  
}

df_impact <- df_impact[, c("inc", "sev", "prev", "mor", "scenario")] %>%
  separate(col = scenario, 
           into = c("Experiment", "Seasonality", "System", "EIR", "Agegroup", "Access", "Timing", "IC50"),
           sep = "_",
           remove = FALSE)

df_impact <- df_impact %>%
  pivot_longer(cols = c(inc, sev, prev, mor), names_to = "Outcome", values_to = "Median_Reduction")


# Scale total effects by median impact for each setting

df <- merge(df[, !(names(df) == "scenario")],
            df_impact, 
            by = c("Experiment", "Seasonality", "System", "EIR", "Agegroup", "Access", "Timing", "IC50", "Outcome"))

df$T_eff_scaled <- df$T_eff * df$Median_Reduction


# ----------------------------------------------------------
# Import baseline prevalence
# ----------------------------------------------------------

# Import from csv
prev <- read.csv(paste0("./../Experiments/", exp, "/Outputs/Prevalence_prior_to_intervention.csv"))
prev <- prev[, c("Seasonality", "EIR", "Access", "MaxAge", "annual_prev_210_2034")]
names(prev) <- c("Seasonality", "EIR", "Access", "Agegroup", "annual_prev")

# Merge into data
df <- merge(df, prev, by = c("Seasonality", "EIR", "Access", "Agegroup"))


# ----------------------------------------------------------
# Format data for plotting - MUST BE ADJUSTED FOR EACH EXPERIMENT
# ----------------------------------------------------------

df$Seasonality <- factor(df$Seasonality, levels = c("seas3mo", "seas5mo"))
df$Seasonality <- recode(df$Seasonality, "seas3mo" = "3 MONTH SEASON", "seas5mo" = "5 MONTH SEASON")

df$Access <- factor(df$Access, levels = c("0.04", "0.24"))
df$Access <- recode(df$Access, "0.04" = "LOW ACCESS", "0.24" = "HIGH ACCESS")

df$Setting <- paste0(df$Seasonality, "\n", df$Access)
df$Setting <- factor(df$Setting, levels = c("5 MONTH SEASON\nLOW ACCESS", 
                                            "5 MONTH SEASON\nHIGH ACCESS",
                                            "3 MONTH SEASON\nLOW ACCESS",
                                            "3 MONTH SEASON\nHIGH ACCESS"))

df$Agegroup <- factor(df$Agegroup, levels = c("5", "10"))
df$Agegroup <- recode(df$Agegroup,
                      "5" = "CHILDREN 3 TO 59 MONTHS",
                      "10" = "CHILDREN 3 TO 119 MONTHS")

df$parameter <- factor(df$parameter, levels = c("Coverage1", "Halflife", "Coverage2", "MaxKillingRate", "Slope"))
df$parameter <- recode(df$parameter,
                       "Coverage1" = "Program reach [70% - 95%]",
                       "Coverage2" = "Round coverage [70% - 95%]",
                       "MaxKillingRate" = "Emax [2 - 30 units]",
                       "Halflife" = "Elimination half-life [1 - 20 days]",
                       "Slope" = "Slope [6 - 6]")

df$Outcome <- factor(df$Outcome, levels = c("inc", "prev", "sev", "mor"))
df$Outcome <- recode(df$Outcome,
                     "inc" = "CLINICAL INCIDENCE",
                     "prev" = "PREVALENCE",
                     "sev" = "SEVERE DISEASE",
                     "mor" = "MORTALITY")

index <- order(unique(round(df$annual_prev*100, 0)))
df$annual_prev_lab <- paste0(round(df$annual_prev*100, 0), "%")
df$annual_prev_lab <- factor(df$annual_prev_lab, levels = unique(df$annual_prev_lab)[index])

df$label <- ifelse(df$T_eff >= 0.07, paste0(round(df$T_eff*100, 0), "%"), "")
df$label <- ifelse(df$Median_Reduction  >= 20, df$label, "")

# ----------------------------------------------------------
# Write data to file
# ----------------------------------------------------------

saveRDS(df, "./data_and_visualisation/Appendix_Figure21/data_figA21_panelB.rds")
