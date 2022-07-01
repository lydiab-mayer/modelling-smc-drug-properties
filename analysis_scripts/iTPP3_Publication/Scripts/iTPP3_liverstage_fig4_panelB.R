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

library(tidyr)
library(dplyr)
library(hetGP)
library(ggplot2)
library(patchwork)

# !!! Insert your experiment name here as a string, e.g. "MyExperiment" !!!
exp_list <-c("iTPP3_ChemoLiver_TreatLiverBlood_3rounds", "iTPP3_ChemoLiver_TreatLiverBlood_4rounds", "iTPP3_ChemoLiver_TreatLiverBlood_5rounds")
num_rounds <- c("iTPP3_ChemoLiver_TreatLiverBlood_3rounds" = 3, "iTPP3_ChemoLiver_TreatLiverBlood_4rounds" = 4, "iTPP3_ChemoLiver_TreatLiverBlood_5rounds" = 5)

# !!! Insert your predicted parameter here. Note that this must match with one column name in post-processing files !!!
pred_list <- c("inc_red_int_Tot", "sev_red_int_Tot")

# !!! Identify desired targets for optimisation !!!
targets <- seq(50, 90, by = 5)

# !!! Identify desired levels of coverage !!!
coverage1 <- seq(0.75, 0.95, by = 0.1)
coverage2 <- seq(0.75, 0.95, by = 0.1)

user <- strsplit(getwd(), "/", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]][5]
GROUP_dr <- "/scicore/home/penny/GROUP/M3TPP/"
setwd(paste0("/scicore/home/penny/", user, "/M3TPP/"))

varnames <- c("Halflife", "Efficacy")
out <- data.frame()

# Define plot colours
cols <- c("#0f1f2f", "#1d3d5d", "#295784", "#3e6897", "#537bac", "#938998", "#ca9574",
          "#f9a24b", "#fab464", "#fac67c", "#fad694", "#fcdfaa", "#fee8c0")


# ----------------------------------------------------------
# Calculate cutoff criteria
# ----------------------------------------------------------

for(exp in exp_list) {
  
  # Set up
  if (!dir.exists(paste0("./Experiments/", exp, "/Outputs"))) dir.create(paste0("./Experiments/", exp, "/Outputs"))
  
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
    out_int <- expand.grid(setting_id[i], targets, coverage1, coverage2)
    out_int <- cbind(out_int, matrix(nrow = nrow(out_int), ncol = length(varnames)))
    colnames(out_int) <- c("Scenario", "Target", "Coverage1", "Coverage2", varnames)
    
    # Read in optimisation result
    df <- readRDS(setting[i])$scenario
    
    # Add confidence interval bounds
    df$se <- sqrt(df$sd2 + df$nugs)
    df$cl <- qnorm(0.05, df$mean, df$se); df$cl[df$cl <= 0] <- 0
    df$cu <- qnorm(0.95, df$mean, df$se); df$cu[df$cu >= 100] <- 100
    
    # Extract minimum criteria for selected targets based on lower bound of 95% CI
    for (j in 1:nrow(out_int)) {
      for (k in varnames) {
        temp <- df[df$Coverage1 == out_int[j, "Coverage1"] & 
                     df$Coverage2 == out_int[j, "Coverage2"], ]
        out_int[j, k] <- min(temp[temp$cl >= out_int[j, "Target"], k])
      }
    }
    
    # Store outputs
    out <- rbind(out, out_int)
  }

  # Save resulting matrix against future use
  saveRDS(out, paste0("./analysisworkflow/analysis_scripts/iTPP3_Publication/Figures/", exp, "_cutoff_criteria.rds"))
  
}
 
# ----------------------------------------------------------
# Generate and format database of results across all experiments
# ----------------------------------------------------------

# Merge results across experiments

# !!! Note that this must be updated manually should the list of experiments be changed !!!
df3 <- readRDS(paste0("./analysisworkflow/analysis_scripts/iTPP3_Publication/Figures/", exp_list[1], "_cutoff_criteria.rds"))
df4 <- readRDS(paste0("./analysisworkflow/analysis_scripts/iTPP3_Publication/Figures/", exp_list[2], "_cutoff_criteria.rds"))
df5 <- readRDS(paste0("./analysisworkflow/analysis_scripts/iTPP3_Publication/Figures/", exp_list[3], "_cutoff_criteria.rds"))

df <- rbind(df3, df4, df5)
remove(df3, df4, df5)
 
# Replace infs with NAs
df <- data.frame(lapply(df, function(x) gsub("Inf", NA, x)), stringsAsFactors = FALSE)
  
# Format variables as numeric
df$Target <- as.numeric(df$Target)
df$Coverage1 <- as.numeric(df$Coverage1)
df$Coverage2 <- as.numeric(df$Coverage2)
df$Halflife <- as.numeric(df$Halflife)
df$Efficacy <- as.numeric(df$Efficacy)

# Replace minimum characterstics outside of considered range with NA
df$Halflife[df$Halflife >= 60] <- NA
df$Efficacy[df$Efficacy >= 1] <- NA
  
# Add additional columns for each scenario factor
df$Scenario <- sub("_red_int_Tot", "", sub("iTPP3ChemoLiverTreatLiverBlood", "", df$Scenario))
df <- df %>%
  separate(col = Scenario,
           into = c("Rounds", "Seasonality", "System", "EIR", "Agegroup", "Decay", "Access", "Timing", "Endpoint"),
           sep= "_")
# 
# # Transform coverage variables
# df$Rounds <- as.numeric(sub("rounds", "", df$Rounds))
# df$CoverageAllRounds <- df$Coverage1*df$Coverage2^df$Rounds
# df$CoverageOneRound <- 1 - (1 - df$Coverage1) + df$Coverage1*((1 - df$Coverage2)^df$Rounds)
  
# Merge in baseline prevalence
prev <- read.csv(paste0("./Experiments/", exp_list[2], "/Outputs/Prevalence_prior_to_intervention.csv"))
prev <- prev[, c("Seasonality", "EIR", "Access", "MaxAge", "annual_prev_210_2034")]
names(prev) <- c("Seasonality", "EIR", "Access", "Agegroup", "annual_prev_210_2034")

### NOTE that there are small (<1%) differences between baseline prevalence per agegroup, due to stochasticity in OpenMalaria. To have consistent plotting, take the higher value
prev <- prev %>%
  group_by(Seasonality, EIR, Access) %>%
  summarise(annual_prev_210_2034 = max(annual_prev_210_2034))

prev$annual_prev <- paste0(round(prev$annual_prev_210_2034*100, 0), "%")

df <- merge(df, prev, by = c("Seasonality", "EIR", "Access"))
  
# Format text
prev_levels <- paste0(unique(round(prev$annual_prev_210_2034[order(prev$annual_prev_210_2034)]*100, 0)), "%")
df$annual_prev <- factor(df$annual_prev, levels = prev_levels)

df$Endpoint <- factor(df$Endpoint, levels = c("inc", "sev", "prev", "mor"))
df$Endpoint <- recode(df$Endpoint,
                      "inc" = "CLINICAL INCIDENCE REDUCTION",
                      "sev" = "SEVERE DISEASE REDUCTION",
                      "prev" = "PREVALENCE REDUCTION",
                      "mor" = "MORTALITY REDUCTION")

df$Rounds <- sub("iTPP3ChemoBloodTreatLiver", "", df$Rounds)  
df$Rounds <- factor(df$Rounds, levels = c("3rounds", "4rounds", "5rounds"))
df$Rounds <- recode(df$Rounds,
                     "3rounds" = "3 SMC ROUNDS",
                     "4rounds" = "4 SMC ROUNDS",
                     "5rounds" = "5 SMC ROUNDS")
  
df$Seasonality <- factor(df$Seasonality, levels = c("seas3mo", "seas5mo"))
df$Seasonality <- recode(df$Seasonality, "seas3mo" = "3 MONTH SEASON", "seas5mo" = "5 MONTH SEASON")
  
df$Agegroup <- factor(df$Agegroup, levels = c("5", "10"))
df$Agegroup <- recode(df$Agegroup,
                       "5" = "CHILDREN 3 TO 59 MONTHS",
                       "10" = "CHILDREN 3 TO 119 MONTHS")
  
df$Access <- factor(df$Access, levels = c("0.04", "0.24"))
df$Access <- recode(df$Access,
                     "0.04" = "LOW ACCESS",
                     "0.24" = "HIGH ACCESS")
  
head(df)


# ----------------------------------------------------------
# Plot results - figure 4 panel B
# ----------------------------------------------------------

# Subset data
df_plot <- df %>%
  filter(Seasonality == "5 MONTH SEASON",
         Access == "HIGH ACCESS",
         Agegroup == "CHILDREN 3 TO 59 MONTHS",
         Rounds == "4 SMC ROUNDS",
         Endpoint == "CLINICAL INCIDENCE REDUCTION")

# Transform min criteria to factor variables
df_plot <- df_plot %>%
  mutate(Target = paste0(Target, "%"),
         Coverage1 = as.factor(paste0(Coverage1*100, "% PROGRAM REACH")),
         Coverage2 = as.factor(paste0(Coverage2*100, "% ROUND\nCOVERAGE")),
         Halflife = floor(Halflife/5)*5,
         Efficacy = floor(Efficacy/5)*5)

# Correct minimum values and set factor levels
df_plot <- df_plot %>%
  mutate(Coverage2 = factor(Coverage2, levels = rev(levels(Coverage2))),
         Halflife = as.factor(ifelse(Halflife <= 5, 5, Halflife)),
         Efficacy = as.factor(ifelse(Efficacy <= 5, 1, Efficacy)))


# Generate plot
text_size <- 10

p <- ggplot(df_plot, aes(x = annual_prev, y = Target, fill = Halflife))

p <- p + geom_tile(colour = "white")

p <- p + facet_grid(Coverage2 ~ Coverage1)

p <- p + theme(panel.border = element_blank(), 
               panel.background = element_blank(),
               panel.grid = element_blank(),
               text = element_text(family = "Times New Roman", size = text_size),
               strip.background = element_blank(),
               strip.text = element_text(face = "bold"),
               axis.line = element_blank(),
               axis.ticks = element_blank(),
               axis.text.x = element_text(margin = margin(t = 0)),
               axis.text.y = element_text(margin = margin(r = 0)),
               legend.title = element_text(face = "bold", size = text_size),
               legend.title.align = 0.5,
               legend.position = "bottom")

p <- p + scale_fill_manual(values = cols[c(1, 2, 4, 5, 6, 7, 8, 9, 10, 12)], 
                    na.value = "light grey",
                    labels = c(levels(df_plot$Halflife), 
                               "Target not met in\nparameter space")) +
  scale_y_discrete(labels = c("50%", "", "60%", "", "70%", "", "80%", "", "90%"))

p <- p + labs(y = "TARGET  REDUCTION", x = expression(paste("BASELINE ANNUAL ", italic("Pf"), "PR"["2-10"]))) +
  guides(fill = guide_legend(title.position = "top", title = "DURATION  OF  PROTECTION  HALF-LIFE  (DAYS)", nrow = 1))

p + plot_annotation(title = "B.  MINIMUM  DURATION  OF  PROTECTION  CRITERIA") & 
  theme(plot.title = element_text(family = "Times New Roman", face = "bold", size = text_size))

ggsave(filename = paste0("./analysisworkflow/analysis_scripts/iTPP3_Publication/Figures/fig4_panelB_", exp, ".jpg"),
       plot = last_plot(),
       width = 9.1,
       height = 4.5,
       dpi = 400)


