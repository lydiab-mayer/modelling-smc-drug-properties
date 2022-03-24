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
exp_list <-c("iTPP3_bloodstage_4rounds", "iTPP3_bloodstage_5rounds")
num_rounds <- c("iTPP3_bloodstage_4rounds" = 4, "iTPP3_bloodstage_5rounds" = 5)

# !!! Insert your predicted parameter here. Note that this must match with one column name in post-processing files !!!
pred_list <- c("inc_red_int_Tot", "sev_red_int_Tot")

# !!! Identify desired targets for optimisation !!!
targets <- seq(60, 100, by = 5)
library(tidyr)
library(dplyr)
library(hetGP)
library(ggplot2)
library(patchwork)

user <- strsplit(getwd(), "/", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]][5]
GROUP_dr <- "/scicore/home/penny/GROUP/M3TPP/"
setwd(paste0("/scicore/home/penny/", user, "/M3TPP/"))

# 
# load(paste0(GROUP_dr, exp, "/param_ranges.RData"))
# param_ranges_cont

varnames <- c("CoverageAllRounds", "CoverageOneRound", "Halflife", "MaxKillingRate", "Slope")

out <- matrix(nrow = 1, ncol = 2 + length(varnames))
out[1, ] <- c("Scenario", "Target", varnames)

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
    out_int <- matrix(nrow = 1 + length(targets), ncol = 2 + length(varnames))
    out_int[1, ] <- c("Scenario", "Target", varnames)
    out_int[2:(1 + length(targets)), 1] <- setting_id[i]
    out_int[2:(1 + length(targets)), 2] <- targets
    
    # Read in optimisation result
    df <- readRDS(setting[i])$scenario
    df$CoverageAllRounds <- df$Coverage1*df$Coverage2^num_rounds[exp]
    df$CoverageOneRound <- 1 - (1 - df$Coverage1) + df$Coverage1*((1 - df$Coverage2)^num_rounds[exp])
    
    # Add confidence interval bounds
    df$se <- sqrt(df$sd2 + df$nugs)
    df$cl <- qnorm(0.05, df$mean, df$se); df$cl[df$cl <= 0] <- 0
    df$cu <- qnorm(0.95, df$mean, df$se); df$cu[df$cu >= 100] <- 100
    
    # Extract minimum criteria for selected targets based on lower bound of 95% CI
    for (j in 1:length(targets)){
      for (k in 1:length(varnames)) {
        out_int[j + 1, 2 + k] <- min(df[df$cl >= targets[j], varnames[k]])
      }
    }
    
    # Store outputs
    out <- rbind(out, out_int[2:(1 + length(targets)), ])
  }

  # Save resulting matrix against future use
  saveRDS(out, paste0("./analysisworkflow/analysis_scripts/iTPP3_Publication/Figures/", exp, "_cutoff_criteria.rds"))
  #out <- readRDS(paste0("./analysisworkflow/analysis_scripts/iTPP3_Publication/Figures/bloodstage_cutoff_criteria.rds"))
}
 

 
# ----------------------------------------------------------
# Generate and format database of results across all experiments
# ----------------------------------------------------------

# Merge results across experiments

# !!! Note that this must be updated manually should the list of experiments be changed !!!
#df3 <- readRDS(paste0("./analysisworkflow/analysis_scripts/iTPP3_Publication/Figures/iTPP3_bloodstage_3rounds_cutoff_criteria.rds"))
df4 <- readRDS(paste0("./analysisworkflow/analysis_scripts/iTPP3_Publication/Figures/iTPP3_bloodstage_4rounds_cutoff_criteria.rds"))
df5 <- readRDS(paste0("./analysisworkflow/analysis_scripts/iTPP3_Publication/Figures/iTPP3_bloodstage_5rounds_cutoff_criteria.rds"))

#df <- rbind(df3, df4[-1, ], df5[-1, ])
df <- rbind(df4, df5[-1, ])
remove(df4, df5)

# Add column names
colnames(df) <- df[1, ]
df <- as.data.frame(df[-1, ], stringsAsFactors = FALSE)
  
# Replace infs with NAs
df <- data.frame(lapply(df, function(x) gsub("Inf", NA, x)), stringsAsFactors = FALSE)
  
# Format performance characteristics as numeric
df$Target <- as.numeric(df$Target)
df$CoverageAllRounds <- as.numeric(df$CoverageAllRounds)
df$CoverageOneRound <- as.numeric(df$CoverageOneRound)
df$Halflife <- as.numeric(df$Halflife)
df$MaxKillingRate <- as.numeric(df$MaxKillingRate)  
df$Slope <- as.numeric(df$Slope)  

# Replace minimum characterstics outside of considered range with NA
df$Halflife[df$Halflife >= 40] <- NA
df$MaxKillingRate[df$MaxKillingRate >= 30] <- NA
  
# Add additional columns for each scenario factor
df$Scenario <- sub("_red_int_Tot", "", sub("iTPP3bloodstage", "", df$Scenario))
df <- df %>%
  separate(col = Scenario,
           into = c("Rounds", "Seasonality", "System", "EIR", "Agegroup", "Access", "Timing", "IC50", "Endpoint"),
           sep= "_")
  
# Merge in baseline prevalence
prev <- read.csv(paste0("./Experiments/", exp_list[1], "/Outputs/Prevalence_prior_to_intervention.csv"))
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
  
df$Rounds <- sub("iTPP3bloodstage", "", df$Rounds)
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
  
saveRDS(df, paste0("./analysisworkflow/analysis_scripts/iTPP3_Publication/Figures/cutoff_criteria.rds"))
  
  
# ----------------------------------------------------------
# Plot results - annex for clinical incidence reduction
# ----------------------------------------------------------

# Set up plot parameters and data

text_size <- 10
df_plot <- df[df$Endpoint == "CLINICAL INCIDENCE REDUCTION", ]
  

### Generate plot for Halflife ###
  
p <- ggplot(df_plot, aes(x = annual_prev, y = Target, fill = Halflife))
  
p <- p + geom_tile(colour = "white")

p <- p + facet_wrap(Rounds + Agegroup ~ Seasonality + Access, scales = "free", ncol = 4)
  
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
               plot.title = element_text(hjust = 0.5, face = "bold"),
               legend.key.width = unit(1, "cm"),
               legend.title = element_text(face = "bold", size = text_size),
               legend.position = "bottom")
  
p <- p + scale_y_continuous(labels = scales::number_format(accuracy = 1, suffix = "%"), breaks = targets) +
  scale_fill_stepsn(colours = c("#1b4d79", "#5880b1", "#f9a24b", "#f9d48e", "#ffedcb"),
                    na.value = "light grey",
                    breaks = seq(5, 40, by = 5),
                    #labels = function(x) paste0(x, "%"),
                    limits = c(5, 40),
                    show.limits = TRUE)
  
p <- p + labs(title = "MINIMUM CRITERIA FOR CLINICAL INCIDENCE REDUCTION",  y = "TARGET REDUCTION", x = expression(paste("BASELINE ANNUAL ", italic("Pf"), "PR"["2-10"]))) +
  guides(fill = guide_colourbar(title.position = "top", title.hjust = 0.5, title = "CONCENTRATION HALFLIFE"))
  
p

ggsave(filename = paste0("./analysisworkflow/analysis_scripts/iTPP3_Publication/Figures/Annex_iTPP3_cutoff_criteria_INC_halflife.jpg"),
       plot = last_plot(),
       width = 9,
       height = 12,
       dpi = 300)


### Generate plot for MaxKillingRate ###

p <- ggplot(df_plot, aes(x = annual_prev, y = Target, fill = MaxKillingRate))

p <- p + geom_tile(colour = "white")

p <- p + facet_wrap(Rounds + Agegroup ~ Seasonality + Access, scales = "free", ncol = 4)

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
               plot.title = element_text(hjust = 0.5, face = "bold"),
               legend.key.width = unit(1, "cm"),
               legend.title = element_text(face = "bold", size = text_size),
               legend.position = "bottom")

p <- p + scale_y_continuous(labels = scales::number_format(accuracy = 1, suffix = "%"), breaks = targets) +
  scale_fill_stepsn(colours = c("#1b4d79", "#5880b1", "#f9a24b", "#f9d48e", "#ffedcb"),
                    na.value = "light grey",
                    breaks = seq(5, 30, by = 5),
                    #labels = function(x) paste0(x, "%"),
                    limits = c(5, 30),
                    show.limits = TRUE)

p <- p + labs(title = "MINIMUM CRITERIA FOR CLINICAL INCIDENCE REDUCTION",  y = "TARGET REDUCTION", x = expression(paste("BASELINE ANNUAL ", italic("Pf"), "PR"["2-10"]))) +
  guides(fill = guide_colourbar(title.position = "top", title.hjust = 0.5, title = "MAXIMUM PARASITE KILLING EFFECT"))

p

ggsave(filename = paste0("./analysisworkflow/analysis_scripts/iTPP3_Publication/Figures/Annex_iTPP3_cutoff_criteria_INC_maxkilling.jpg"),
       plot = last_plot(),
       width = 9,
       height = 12,
       dpi = 300)


### Generate plot for Coverage of all rounds ###

p <- ggplot(df_plot, aes(x = annual_prev, y = Target, fill = CoverageAllRounds))

p <- p + geom_tile(colour = "white")

p <- p + facet_wrap(Rounds + Agegroup ~ Seasonality + Access, scales = "free", ncol = 4)

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
               plot.title = element_text(hjust = 0.5, face = "bold"),
               legend.key.width = unit(1, "cm"),
               legend.title = element_text(face = "bold", size = text_size),
               legend.position = "bottom")

p <- p + scale_y_continuous(labels = scales::number_format(accuracy = 1, suffix = "%"), breaks = targets) +
  scale_fill_stepsn(colours = c("#1b4d79", "#5880b1", "#f9a24b", "#f9d48e", "#ffedcb"),
                    na.value = "light grey",
                    breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0),
                    labels = function(x) paste0(x*100, "%"),
                    limits = c(0, 1.0),
                    show.limits = TRUE)

p <- p + labs(title = "MINIMUM CRITERIA FOR CLINICAL INCIDENCE REDUCTION",  y = "TARGET REDUCTION", x = expression(paste("BASELINE ANNUAL ", italic("Pf"), "PR"["2-10"]))) +
  guides(fill = guide_colourbar(title.position = "top", title.hjust = 0.5, title = "COVERAGE OF ALL SMC ROUNDS"))

p

ggsave(filename = paste0("./analysisworkflow/analysis_scripts/iTPP3_Publication/Figures/Annex_iTPP3_cutoff_criteria_INC_coverageallrounds.jpg"),
       plot = last_plot(),
       width = 9,
       height = 12,
       dpi = 300)


### Generate plot for Coverage of > one round ###

p <- ggplot(df_plot, aes(x = annual_prev, y = Target, fill = CoverageOneRound))

p <- p + geom_tile(colour = "white")

p <- p + facet_wrap(Rounds + Agegroup ~ Seasonality + Access, scales = "free", ncol = 4)

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
               plot.title = element_text(hjust = 0.5, face = "bold"),
               legend.key.width = unit(1, "cm"),
               legend.title = element_text(face = "bold", size = text_size),
               legend.position = "bottom")

p <- p + scale_y_continuous(labels = scales::number_format(accuracy = 1, suffix = "%"), breaks = targets) +
  scale_fill_stepsn(colours = c("#1b4d79", "#5880b1", "#f9a24b", "#f9d48e", "#ffedcb"),
                    na.value = "light grey",
                    breaks = seq(0.7, 1.0, by = 0.05),
                    labels = function(x) paste0(x*100, "%"),
                    limits = c(0.7, 1.0),
                    show.limits = TRUE)

p <- p + labs(title = "MINIMUM CRITERIA FOR CLINICAL INCIDENCE REDUCTION",  y = "TARGET REDUCTION", x = expression(paste("BASELINE ANNUAL ", italic("Pf"), "PR"["2-10"]))) +
  guides(fill = guide_colourbar(title.position = "top", title.hjust = 0.5, title = "COVERAGE OF AT LEAST ONE SMC ROUND"))

p

ggsave(filename = paste0("./analysisworkflow/analysis_scripts/iTPP3_Publication/Figures/Annex_iTPP3_cutoff_criteria_INC_coverageoneround.jpg"),
       plot = last_plot(),
       width = 9,
       height = 12,
       dpi = 300)



# ----------------------------------------------------------
# Plot results - annex for severe disease reduction
# ----------------------------------------------------------

# Set up plot parameters and data

text_size <- 10
df_plot <- df[df$Endpoint == "SEVERE DISEASE REDUCTION", ]


### Generate plot for Halflife ###

p <- ggplot(df_plot, aes(x = annual_prev, y = Target, fill = Halflife))

p <- p + geom_tile(colour = "white")

p <- p + facet_wrap(Rounds + Agegroup ~ Seasonality + Access, scales = "free", ncol = 4)

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
               plot.title = element_text(hjust = 0.5, face = "bold"),
               legend.key.width = unit(1, "cm"),
               legend.title = element_text(face = "bold", size = text_size),
               legend.position = "bottom")

p <- p + scale_y_continuous(labels = scales::number_format(accuracy = 1, suffix = "%"), breaks = targets) +
  scale_fill_stepsn(colours = c("#1b4d79", "#5880b1", "#f9a24b", "#f9d48e", "#ffedcb"),
                    na.value = "light grey",
                    breaks = seq(5, 40, by = 5),
                    #labels = function(x) paste0(x, "%"),
                    limits = c(5, 40),
                    show.limits = TRUE)

p <- p + labs(title = "MINIMUM CRITERIA FOR SEVERE DISEASE REDUCTION",  y = "TARGET REDUCTION", x = expression(paste("BASELINE ANNUAL ", italic("Pf"), "PR"["2-10"]))) +
  guides(fill = guide_colourbar(title.position = "top", title.hjust = 0.5, title = "CONCENTRATION HALFLIFE"))

p

ggsave(filename = paste0("./analysisworkflow/analysis_scripts/iTPP3_Publication/Figures/Annex_iTPP3_cutoff_criteria_SEV_halflife.jpg"),
       plot = last_plot(),
       width = 9,
       height = 12,
       dpi = 300)


### Generate plot for MaxKillingRate ###

p <- ggplot(df_plot, aes(x = annual_prev, y = Target, fill = MaxKillingRate))

p <- p + geom_tile(colour = "white")

p <- p + facet_wrap(Rounds + Agegroup ~ Seasonality + Access, scales = "free", ncol = 4)

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
               plot.title = element_text(hjust = 0.5, face = "bold"),
               legend.key.width = unit(1, "cm"),
               legend.title = element_text(face = "bold", size = text_size),
               legend.position = "bottom")

p <- p + scale_y_continuous(labels = scales::number_format(accuracy = 1, suffix = "%"), breaks = targets) +
  scale_fill_stepsn(colours = c("#1b4d79", "#5880b1", "#f9a24b", "#f9d48e", "#ffedcb"),
                    na.value = "light grey",
                    breaks = seq(5, 30, by = 5),
                    #labels = function(x) paste0(x, "%"),
                    limits = c(5, 30),
                    show.limits = TRUE)

p <- p + labs(title = "MINIMUM CRITERIA FOR  SEVERE DISEASE REDUCTION",  y = "TARGET REDUCTION", x = expression(paste("BASELINE ANNUAL ", italic("Pf"), "PR"["2-10"]))) +
  guides(fill = guide_colourbar(title.position = "top", title.hjust = 0.5, title = "MAXIMUM PARASITE KILLING EFFECT"))

p

ggsave(filename = paste0("./analysisworkflow/analysis_scripts/iTPP3_Publication/Figures/Annex_iTPP3_cutoff_criteria_SEV_maxkilling.jpg"),
       plot = last_plot(),
       width = 9,
       height = 12,
       dpi = 300)


### Generate plot for Coverage of all rounds ###

p <- ggplot(df_plot, aes(x = annual_prev, y = Target, fill = CoverageAllRounds))

p <- p + geom_tile(colour = "white")

p <- p + facet_wrap(Rounds + Agegroup ~ Seasonality + Access, scales = "free", ncol = 4)

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
               plot.title = element_text(hjust = 0.5, face = "bold"),
               legend.key.width = unit(1, "cm"),
               legend.title = element_text(face = "bold", size = text_size),
               legend.position = "bottom")

p <- p + scale_y_continuous(labels = scales::number_format(accuracy = 1, suffix = "%"), breaks = targets) +
  scale_fill_stepsn(colours = c("#1b4d79", "#5880b1", "#f9a24b", "#f9d48e", "#ffedcb"),
                    na.value = "light grey",
                    breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1.0),
                    labels = function(x) paste0(x*100, "%"),
                    limits = c(0, 1.0),
                    show.limits = TRUE)

p <- p + labs(title = "MINIMUM CRITERIA FOR  SEVERE DISEASE REDUCTION",  y = "TARGET REDUCTION", x = expression(paste("BASELINE ANNUAL ", italic("Pf"), "PR"["2-10"]))) +
  guides(fill = guide_colourbar(title.position = "top", title.hjust = 0.5, title = "COVERAGE OF ALL SMC ROUNDS"))

p

ggsave(filename = paste0("./analysisworkflow/analysis_scripts/iTPP3_Publication/Figures/Annex_iTPP3_cutoff_criteria_SEV_coverageallrounds.jpg"),
       plot = last_plot(),
       width = 9,
       height = 12,
       dpi = 300)


### Generate plot for Coverage of > one round ###

p <- ggplot(df_plot, aes(x = annual_prev, y = Target, fill = CoverageOneRound))

p <- p + geom_tile(colour = "white")

p <- p + facet_wrap(Rounds + Agegroup ~ Seasonality + Access, scales = "free", ncol = 4)

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
               plot.title = element_text(hjust = 0.5, face = "bold"),
               legend.key.width = unit(1, "cm"),
               legend.title = element_text(face = "bold", size = text_size),
               legend.position = "bottom")

p <- p + scale_y_continuous(labels = scales::number_format(accuracy = 1, suffix = "%"), breaks = targets) +
  scale_fill_stepsn(colours = c("#1b4d79", "#5880b1", "#f9a24b", "#f9d48e", "#ffedcb"),
                    na.value = "light grey",
                    breaks = seq(0.7, 1.0, by = 0.05),
                    labels = function(x) paste0(x*100, "%"),
                    limits = c(0.7, 1.0),
                    show.limits = TRUE)

p <- p + labs(title = "MINIMUM CRITERIA FOR  SEVERE DISEASE REDUCTION",  y = "TARGET REDUCTION", x = expression(paste("BASELINE ANNUAL ", italic("Pf"), "PR"["2-10"]))) +
  guides(fill = guide_colourbar(title.position = "top", title.hjust = 0.5, title = "COVERAGE OF AT LEAST ONE SMC ROUND"))

p

ggsave(filename = paste0("./analysisworkflow/analysis_scripts/iTPP3_Publication/Figures/Annex_iTPP3_cutoff_criteria_SEV_coverageoneround.jpg"),
       plot = last_plot(),
       width = 9,
       height = 12,
       dpi = 300)



# ----------------------------------------------------------
# Plot results - executive summary - 5 month season and high access
# ----------------------------------------------------------

text_size <- 10
df_plot <- df

# Calculate most conservative across all scenarios and all outcomes
df_plot <- df_plot %>%
  group_by(Access, Seasonality, Target, annual_prev) %>%
  summarise(sum_na = sum(is.na(CoverageAllRounds)),
            CoverageAllRounds = max(CoverageAllRounds, na.rm = FALSE),
            CoverageOneRound = max(CoverageOneRound, na.rm = FALSE),
            Halflife = max(Halflife, na.rm = FALSE),
            MaxKillingRate = max(MaxKillingRate, na.rm = FALSE))

# Subset data
df_plot <- df_plot[df_plot$Seasonality == "5 MONTH SEASON" & df_plot$Access == "HIGH ACCESS", ]

# Transform min criteria to factor variables
df_plot$CoverageAllRounds_factor <- floor(df_plot$CoverageAllRounds*100 / 10)*10 
df_plot$CoverageAllRounds_factor <- ifelse(is.na(df_plot$CoverageAllRounds_factor), 
                                           NA,
                                           paste0(df_plot$CoverageAllRounds_factor, "%"))
df_plot$CoverageAllRounds_factor <- as.factor(df_plot$CoverageAllRounds_factor)

df_plot$CoverageOneRound_factor <- floor(df_plot$CoverageOneRound*100 / 5)*5 
df_plot$CoverageOneRound_factor <- ifelse(is.na(df_plot$CoverageOneRound_factor), 
                                          NA,
                                          paste0(df_plot$CoverageOneRound_factor, "%"))
df_plot$CoverageOneRound_factor <- as.factor(df_plot$CoverageOneRound_factor)

df_plot$Halflife_factor <- floor(df_plot$Halflife / 5)*5
df_plot$Halflife_factor <- as.factor(ifelse(df_plot$Halflife_factor <= 5, 5, df_plot$Halflife_factor))

df_plot$MaxKillingRate_factor <- floor(df_plot$MaxKillingRate / 5)*5
df_plot$MaxKillingRate_factor <- as.factor(ifelse(df_plot$MaxKillingRate_factor <= 5, 1, df_plot$MaxKillingRate_factor))


### Generate plot for Coverage of All SMC ROUNDS ###

p <- ggplot(df_plot, aes(x = annual_prev, y = Target, fill = CoverageAllRounds_factor))

p <- p + geom_tile(colour = "white")

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
               plot.title = element_text(hjust = 0.5, face = "bold"),
               legend.key.width = unit(1, "cm"),
               legend.title = element_text(face = "bold", size = text_size),
               legend.title.align = 0.5,
               legend.position = "bottom")

p <- p + scale_y_continuous(labels = scales::number_format(accuracy = 1, suffix = "%"), breaks = targets) +
  scale_fill_manual(values = cols[c(2, 3, 4, 5, 6, 7, 8, 10, 12)], 
                    na.value = "light grey",
                    labels = c(levels(df_plot$CoverageAllRounds_factor), 
                               "Target not met within\nparameter space"))

p <- p + labs(title = "COVERAGE OF ALL SMC ROUNDS",  y = "TARGET  REDUCTION", x = expression(paste("BASELINE ANNUAL ", italic("Pf"), "PR"["2-10"]))) +
  guides(fill = guide_legend(title.position = "top", title.hjust = 0.5, title = "COVERAGE  OF ALL  SMC  ROUNDS", nrow = 2))

p1 <- p


### Generate plot for Coverage of one SMC round ###

p <- ggplot(df_plot, aes(x = annual_prev, y = Target, fill = CoverageOneRound_factor))

p <- p + geom_tile(colour = "white")

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
               plot.title = element_text(hjust = 0.5, face = "bold"),
               legend.key.width = unit(1, "cm"),
               legend.title = element_text(face = "bold", size = text_size),
               legend.title.align = 0.5,
               legend.position = "bottom")

p <- p + scale_y_continuous(labels = scales::number_format(accuracy = 1, suffix = "%"), breaks = targets) +
  scale_fill_manual(values = cols[c(2, 4, 6, 8, 10, 12)], 
                    na.value = "light grey",
                    labels = c(levels(df_plot$CoverageOneRound_factor), 
                               "Target not met within\nparameter space"))

p <- p + labs(title = "COVERAGE OF AT LEAST ONE SMC ROUND",  y = "TARGET  REDUCTION", x = expression(paste("BASELINE ANNUAL ", italic("Pf"), "PR"["2-10"]))) +
  guides(fill = guide_legend(title.position = "top", title.hjust = 0.5, title = "COVERAGE  OF AT  LEAST  ONE  SMC  ROUND", nrow = 2))

p2 <- p


### Generate plot for duration ###

p <- ggplot(df_plot, aes(x = annual_prev, y = Target, fill = Halflife_factor))

p <- p + geom_tile(colour = "white")

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
               plot.title = element_text(hjust = 0.5, face = "bold"),
               legend.key.width = unit(1, "cm"),
               legend.title = element_text(face = "bold", size = text_size),
               legend.title.align = 0.5,
               legend.position = "bottom")

p <- p + scale_y_continuous(labels = scales::number_format(accuracy = 1, suffix = "%"), breaks = targets)  +
  scale_fill_manual(values = cols[c(2, 4, 6, 8, 12)], 
                    na.value = "light grey",
                    labels = c(levels(df_plot$Halflife_factor), 
                               "Target not met within\nparameter space"))

p <- p + labs(title = "CONCENTRATION HALF-LIFE",  y = "TARGET  REDUCTION", x = expression(paste("BASELINE ANNUAL ", italic("Pf"), "PR"["2-10"]))) +
  guides(fill = guide_legend(title.position = "top", title.hjust = 0.5, title = "CONCENTRATION HALF-LIFE  (DAYS)", nrow = 1))

p3 <- p
  
  
# Generate plot for Max Killing Rate
  
p <- ggplot(df_plot, aes(x = annual_prev, y = Target, fill = MaxKillingRate_factor))

p <- p + geom_tile(colour = "white")

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
               plot.title = element_text(hjust = 0.5, face = "bold"),
               legend.key.width = unit(1, "cm"),
               legend.title = element_text(face = "bold", size = text_size),
               legend.title.align = 0.5,
               legend.position = "bottom")

p <- p + scale_y_continuous(labels = scales::number_format(accuracy = 1, suffix = "%"), breaks = targets) +
  scale_fill_manual(values = cols[c(2, 8, 12)], 
                    na.value = "light grey",
                    labels = c(levels(df_plot$MaxKillingRate_factor), 
                               "Target not met within\nparameter space"))

p <- p + labs(title = "MAXIMUM EFFECT",  y = "TARGET  REDUCTION", x = expression(paste("BASELINE ANNUAL ", italic("Pf"), "PR"["2-10"]))) +
  guides(fill = guide_legend(title.position = "top", title.hjust = 0.5, title = "MAXIMUM  PARASITE KILLING  EFFECT", nrow = 1))

p4 <- p
  
p1 + p2 + p3 + p4 + plot_layout(nrow = 2, ncol = 2)
  
ggsave(filename = paste0("./analysisworkflow/analysis_scripts/iTPP3_Publication/Figures/iTPP3_cutoff_criteria_SUMMARY_5moHIGHaccess.jpg"),
       plot = last_plot(),
       width = 9.1,
       height = 6,
       dpi = 400)



# ----------------------------------------------------------
# Plot results - executive summary - 3 month season and high access
# ----------------------------------------------------------

text_size <- 10
df_plot <- df

# Calculate most conservative across all scenarios and all outcomes
df_plot <- df_plot %>%
  group_by(Access, Seasonality, Target, annual_prev) %>%
  summarise(sum_na = sum(is.na(CoverageAllRounds)),
            CoverageAllRounds = max(CoverageAllRounds, na.rm = FALSE),
            CoverageOneRound = max(CoverageOneRound, na.rm = FALSE),
            Halflife = max(Halflife, na.rm = FALSE),
            MaxKillingRate = max(MaxKillingRate, na.rm = FALSE))

# Subset data
df_plot <- df_plot[df_plot$Seasonality == "3 MONTH SEASON" & df_plot$Access == "HIGH ACCESS", ]

# Transform min criteria to factor variables
df_plot$CoverageAllRounds_factor <- floor(df_plot$CoverageAllRounds*100 / 10)*10 
df_plot$CoverageAllRounds_factor <- ifelse(is.na(df_plot$CoverageAllRounds_factor), 
                                           NA,
                                           paste0(df_plot$CoverageAllRounds_factor, "%"))
df_plot$CoverageAllRounds_factor <- as.factor(df_plot$CoverageAllRounds_factor)

df_plot$CoverageOneRound_factor <- floor(df_plot$CoverageOneRound*100 / 5)*5 
df_plot$CoverageOneRound_factor <- ifelse(is.na(df_plot$CoverageOneRound_factor), 
                                          NA,
                                          paste0(df_plot$CoverageOneRound_factor, "%"))
df_plot$CoverageOneRound_factor <- as.factor(df_plot$CoverageOneRound_factor)

df_plot$Halflife_factor <- floor(df_plot$Halflife / 5)*5
df_plot$Halflife_factor <- as.factor(ifelse(df_plot$Halflife_factor <= 5, 5, df_plot$Halflife_factor))

df_plot$MaxKillingRate_factor <- floor(df_plot$MaxKillingRate / 5)*5
df_plot$MaxKillingRate_factor <- as.factor(ifelse(df_plot$MaxKillingRate_factor <= 5, 1, df_plot$MaxKillingRate_factor))


### Generate plot for Coverage of All SMC ROUNDS ###

p <- ggplot(df_plot, aes(x = annual_prev, y = Target, fill = CoverageAllRounds_factor))

p <- p + geom_tile(colour = "white")

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
               plot.title = element_text(hjust = 0.5, face = "bold"),
               legend.key.width = unit(1, "cm"),
               legend.title = element_text(face = "bold", size = text_size),
               legend.title.align = 0.5,
               legend.position = "bottom")

p <- p + scale_y_continuous(labels = scales::number_format(accuracy = 1, suffix = "%"), breaks = targets) +
  scale_fill_manual(values = cols[c(2, 3, 4, 5, 6, 7, 8, 10, 12)], 
                    na.value = "light grey",
                    labels = c(levels(df_plot$CoverageAllRounds_factor), 
                               "Target not met within\nparameter space"))

p <- p + labs(title = "COVERAGE OF ALL SMC ROUNDS",  y = "TARGET  REDUCTION", x = expression(paste("BASELINE ANNUAL ", italic("Pf"), "PR"["2-10"]))) +
  guides(fill = guide_legend(title.position = "top", title.hjust = 0.5, title = "COVERAGE  OF ALL  SMC  ROUNDS", nrow = 2))

p1 <- p


### Generate plot for Coverage of one SMC round ###

p <- ggplot(df_plot, aes(x = annual_prev, y = Target, fill = CoverageOneRound_factor))

p <- p + geom_tile(colour = "white")

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
               plot.title = element_text(hjust = 0.5, face = "bold"),
               legend.key.width = unit(1, "cm"),
               legend.title = element_text(face = "bold", size = text_size),
               legend.title.align = 0.5,
               legend.position = "bottom")

p <- p + scale_y_continuous(labels = scales::number_format(accuracy = 1, suffix = "%"), breaks = targets) +
  scale_fill_manual(values = cols[c(2, 4, 6, 8, 10, 12)], 
                    na.value = "light grey",
                    labels = c(levels(df_plot$CoverageOneRound_factor), 
                               "Target not met within\nparameter space"))

p <- p + labs(title = "COVERAGE OF AT LEAST ONE SMC ROUND",  y = "TARGET  REDUCTION", x = expression(paste("BASELINE ANNUAL ", italic("Pf"), "PR"["2-10"]))) +
  guides(fill = guide_legend(title.position = "top", title.hjust = 0.5, title = "COVERAGE  OF AT  LEAST  ONE  SMC  ROUND", nrow = 2))

p2 <- p


### Generate plot for duration ###

p <- ggplot(df_plot, aes(x = annual_prev, y = Target, fill = Halflife_factor))

p <- p + geom_tile(colour = "white")

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
               plot.title = element_text(hjust = 0.5, face = "bold"),
               legend.key.width = unit(1, "cm"),
               legend.title = element_text(face = "bold", size = text_size),
               legend.title.align = 0.5,
               legend.position = "bottom")

p <- p + scale_y_continuous(labels = scales::number_format(accuracy = 1, suffix = "%"), breaks = targets)  +
  scale_fill_manual(values = cols[c(2, 4, 6, 8, 12)], 
                    na.value = "light grey",
                    labels = c(levels(df_plot$Halflife_factor), 
                               "Target not met within\nparameter space"))

p <- p + labs(title = "CONCENTRATION HALF-LIFE",  y = "TARGET  REDUCTION", x = expression(paste("BASELINE ANNUAL ", italic("Pf"), "PR"["2-10"]))) +
  guides(fill = guide_legend(title.position = "top", title.hjust = 0.5, title = "CONCENTRATION HALF-LIFE  (DAYS)", nrow = 1))

p3 <- p


# Generate plot for Max Killing Rate

p <- ggplot(df_plot, aes(x = annual_prev, y = Target, fill = MaxKillingRate_factor))

p <- p + geom_tile(colour = "white")

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
               plot.title = element_text(hjust = 0.5, face = "bold"),
               legend.key.width = unit(1, "cm"),
               legend.title = element_text(face = "bold", size = text_size),
               legend.title.align = 0.5,
               legend.position = "bottom")

p <- p + scale_y_continuous(labels = scales::number_format(accuracy = 1, suffix = "%"), breaks = targets) +
  scale_fill_manual(values = cols[c(2, 8, 12)], 
                    na.value = "light grey",
                    labels = c(levels(df_plot$MaxKillingRate_factor), 
                               "Target not met within\nparameter space"))

p <- p + labs(title = "MAXIMUM EFFECT",  y = "TARGET  REDUCTION", x = expression(paste("BASELINE ANNUAL ", italic("Pf"), "PR"["2-10"]))) +
  guides(fill = guide_legend(title.position = "top", title.hjust = 0.5, title = "MAXIMUM  PARASITE KILLING  EFFECT", nrow = 1))

p4 <- p

p1 + p2 + p3 + p4 + plot_layout(nrow = 2, ncol = 2)

ggsave(filename = paste0("./analysisworkflow/analysis_scripts/iTPP3_Publication/Figures/iTPP3_cutoff_criteria_SUMMARY_3moHIGHaccess.jpg"),
       plot = last_plot(),
       width = 9.1,
       height = 6,
       dpi = 400)


# ----------------------------------------------------------
# Plot results - executive summary - 5 month season and low access
# ----------------------------------------------------------

text_size <- 10
df_plot <- df

# Calculate most conservative across all scenarios and all outcomes
df_plot <- df_plot %>%
  group_by(Access, Seasonality, Target, annual_prev) %>%
  summarise(sum_na = sum(is.na(CoverageAllRounds)),
            CoverageAllRounds = max(CoverageAllRounds, na.rm = FALSE),
            CoverageOneRound = max(CoverageOneRound, na.rm = FALSE),
            Halflife = max(Halflife, na.rm = FALSE),
            MaxKillingRate = max(MaxKillingRate, na.rm = FALSE))

# Subset data
df_plot <- df_plot[df_plot$Seasonality == "5 MONTH SEASON" & df_plot$Access == "LOW ACCESS", ]

# Transform min criteria to factor variables
df_plot$CoverageAllRounds_factor <- floor(df_plot$CoverageAllRounds*100 / 10)*10 
df_plot$CoverageAllRounds_factor <- ifelse(is.na(df_plot$CoverageAllRounds_factor), 
                                           NA,
                                           paste0(df_plot$CoverageAllRounds_factor, "%"))
df_plot$CoverageAllRounds_factor <- as.factor(df_plot$CoverageAllRounds_factor)

df_plot$CoverageOneRound_factor <- floor(df_plot$CoverageOneRound*100 / 5)*5 
df_plot$CoverageOneRound_factor <- ifelse(is.na(df_plot$CoverageOneRound_factor), 
                                          NA,
                                          paste0(df_plot$CoverageOneRound_factor, "%"))
df_plot$CoverageOneRound_factor <- as.factor(df_plot$CoverageOneRound_factor)

df_plot$Halflife_factor <- floor(df_plot$Halflife / 5)*5
df_plot$Halflife_factor <- as.factor(ifelse(df_plot$Halflife_factor <= 5, 5, df_plot$Halflife_factor))

df_plot$MaxKillingRate_factor <- floor(df_plot$MaxKillingRate / 5)*5
df_plot$MaxKillingRate_factor <- as.factor(ifelse(df_plot$MaxKillingRate_factor <= 5, 1, df_plot$MaxKillingRate_factor))


### Generate plot for Coverage of All SMC ROUNDS ###

p <- ggplot(df_plot, aes(x = annual_prev, y = Target, fill = CoverageAllRounds_factor))

p <- p + geom_tile(colour = "white")

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
               plot.title = element_text(hjust = 0.5, face = "bold"),
               legend.key.width = unit(1, "cm"),
               legend.title = element_text(face = "bold", size = text_size),
               legend.title.align = 0.5,
               legend.position = "bottom")

p <- p + scale_y_continuous(labels = scales::number_format(accuracy = 1, suffix = "%"), breaks = targets) +
  scale_fill_manual(values = cols[c(2, 4, 5, 6, 7, 8, 9, 10, 12)], 
                    na.value = "light grey",
                    labels = c(levels(df_plot$CoverageAllRounds_factor), 
                               "Target not met within\nparameter space"))

p <- p + labs(title = "COVERAGE OF ALL SMC ROUNDS",  y = "TARGET  REDUCTION", x = expression(paste("BASELINE ANNUAL ", italic("Pf"), "PR"["2-10"]))) +
  guides(fill = guide_legend(title.position = "top", title.hjust = 0.5, title = "COVERAGE  OF ALL  SMC  ROUNDS", nrow = 2))

p1 <- p


### Generate plot for Coverage of one SMC round ###

p <- ggplot(df_plot, aes(x = annual_prev, y = Target, fill = CoverageOneRound_factor))

p <- p + geom_tile(colour = "white")

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
               plot.title = element_text(hjust = 0.5, face = "bold"),
               legend.key.width = unit(1, "cm"),
               legend.title = element_text(face = "bold", size = text_size),
               legend.title.align = 0.5,
               legend.position = "bottom")

p <- p + scale_y_continuous(labels = scales::number_format(accuracy = 1, suffix = "%"), breaks = targets) +
  scale_fill_manual(values = cols[c(2, 4, 6, 8, 10, 12)], 
                    na.value = "light grey",
                    labels = c(levels(df_plot$CoverageOneRound_factor), 
                               "Target not met within\nparameter space"))

p <- p + labs(title = "COVERAGE OF AT LEAST ONE SMC ROUND",  y = "TARGET  REDUCTION", x = expression(paste("BASELINE ANNUAL ", italic("Pf"), "PR"["2-10"]))) +
  guides(fill = guide_legend(title.position = "top", title.hjust = 0.5, title = "COVERAGE  OF AT  LEAST  ONE  SMC  ROUND", nrow = 2))

p2 <- p


### Generate plot for duration ###

p <- ggplot(df_plot, aes(x = annual_prev, y = Target, fill = Halflife_factor))

p <- p + geom_tile(colour = "white")

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
               plot.title = element_text(hjust = 0.5, face = "bold"),
               legend.key.width = unit(1, "cm"),
               legend.title = element_text(face = "bold", size = text_size),
               legend.title.align = 0.5,
               legend.position = "bottom")

p <- p + scale_y_continuous(labels = scales::number_format(accuracy = 1, suffix = "%"), breaks = targets)  +
  scale_fill_manual(values = cols[c(2, 8, 12)], 
                    na.value = "light grey",
                    labels = c(levels(df_plot$Halflife_factor), 
                               "Target not met within\nparameter space"))

p <- p + labs(title = "CONCENTRATION HALF-LIFE",  y = "TARGET  REDUCTION", x = expression(paste("BASELINE ANNUAL ", italic("Pf"), "PR"["2-10"]))) +
  guides(fill = guide_legend(title.position = "top", title.hjust = 0.5, title = "CONCENTRATION HALF-LIFE  (DAYS)", nrow = 1))

p3 <- p


# Generate plot for Max Killing Rate

p <- ggplot(df_plot, aes(x = annual_prev, y = Target, fill = MaxKillingRate_factor))

p <- p + geom_tile(colour = "white")

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
               plot.title = element_text(hjust = 0.5, face = "bold"),
               legend.key.width = unit(1, "cm"),
               legend.title = element_text(face = "bold", size = text_size),
               legend.title.align = 0.5,
               legend.position = "bottom")

p <- p + scale_y_continuous(labels = scales::number_format(accuracy = 1, suffix = "%"), breaks = targets) +
  scale_fill_manual(values = cols[c(2, 8, 12)], 
                    na.value = "light grey",
                    labels = c(levels(df_plot$MaxKillingRate_factor), 
                               "Target not met within\nparameter space"))

p <- p + labs(title = "MAXIMUM EFFECT",  y = "TARGET  REDUCTION", x = expression(paste("BASELINE ANNUAL ", italic("Pf"), "PR"["2-10"]))) +
  guides(fill = guide_legend(title.position = "top", title.hjust = 0.5, title = "MAXIMUM  PARASITE KILLING  EFFECT", nrow = 1))

p4 <- p

p1 + p2 + p3 + p4 + plot_layout(nrow = 2, ncol = 2)

ggsave(filename = paste0("./analysisworkflow/analysis_scripts/iTPP3_Publication/Figures/iTPP3_cutoff_criteria_SUMMARY_5moLOWaccess.jpg"),
       plot = last_plot(),
       width = 9.1,
       height = 6,
       dpi = 400)


# ----------------------------------------------------------
# Plot results - executive summary - 3 month season and low access
# ----------------------------------------------------------

text_size <- 10
df_plot <- df

# Calculate most conservative across all scenarios and all outcomes
df_plot <- df_plot %>%
  group_by(Access, Seasonality, Target, annual_prev) %>%
  summarise(sum_na = sum(is.na(CoverageAllRounds)),
            CoverageAllRounds = max(CoverageAllRounds, na.rm = FALSE),
            CoverageOneRound = max(CoverageOneRound, na.rm = FALSE),
            Halflife = max(Halflife, na.rm = FALSE),
            MaxKillingRate = max(MaxKillingRate, na.rm = FALSE))

# Subset data
df_plot <- df_plot[df_plot$Seasonality == "3 MONTH SEASON" & df_plot$Access == "LOW ACCESS", ]

# Transform min criteria to factor variables
df_plot$CoverageAllRounds_factor <- floor(df_plot$CoverageAllRounds*100 / 10)*10 
df_plot$CoverageAllRounds_factor <- ifelse(is.na(df_plot$CoverageAllRounds_factor), 
                                           NA,
                                           paste0(df_plot$CoverageAllRounds_factor, "%"))
df_plot$CoverageAllRounds_factor <- as.factor(df_plot$CoverageAllRounds_factor)

df_plot$CoverageOneRound_factor <- floor(df_plot$CoverageOneRound*100 / 5)*5 
df_plot$CoverageOneRound_factor <- ifelse(is.na(df_plot$CoverageOneRound_factor), 
                                          NA,
                                          paste0(df_plot$CoverageOneRound_factor, "%"))
df_plot$CoverageOneRound_factor <- as.factor(df_plot$CoverageOneRound_factor)

df_plot$Halflife_factor <- floor(df_plot$Halflife / 5)*5
df_plot$Halflife_factor <- as.factor(ifelse(df_plot$Halflife_factor <= 5, 5, df_plot$Halflife_factor))

df_plot$MaxKillingRate_factor <- floor(df_plot$MaxKillingRate / 5)*5
df_plot$MaxKillingRate_factor <- as.factor(ifelse(df_plot$MaxKillingRate_factor <= 5, 1, df_plot$MaxKillingRate_factor))


### Generate plot for Coverage of All SMC ROUNDS ###

p <- ggplot(df_plot, aes(x = annual_prev, y = Target, fill = CoverageAllRounds_factor))

p <- p + geom_tile(colour = "white")

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
               plot.title = element_text(hjust = 0.5, face = "bold"),
               legend.key.width = unit(1, "cm"),
               legend.title = element_text(face = "bold", size = text_size),
               legend.title.align = 0.5,
               legend.position = "bottom")

p <- p + scale_y_continuous(labels = scales::number_format(accuracy = 1, suffix = "%"), breaks = targets) +
  scale_fill_manual(values = cols[c(2, 4, 5, 6, 7, 8, 9, 10, 12)], 
                    na.value = "light grey",
                    labels = c(levels(df_plot$CoverageAllRounds_factor), 
                               "Target not met within\nparameter space"))

p <- p + labs(title = "COVERAGE OF ALL SMC ROUNDS",  y = "TARGET  REDUCTION", x = expression(paste("BASELINE ANNUAL ", italic("Pf"), "PR"["2-10"]))) +
  guides(fill = guide_legend(title.position = "top", title.hjust = 0.5, title = "COVERAGE  OF ALL  SMC  ROUNDS", nrow = 2))

p1 <- p


### Generate plot for Coverage of one SMC round ###

p <- ggplot(df_plot, aes(x = annual_prev, y = Target, fill = CoverageOneRound_factor))

p <- p + geom_tile(colour = "white")

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
               plot.title = element_text(hjust = 0.5, face = "bold"),
               legend.key.width = unit(1, "cm"),
               legend.title = element_text(face = "bold", size = text_size),
               legend.title.align = 0.5,
               legend.position = "bottom")

p <- p + scale_y_continuous(labels = scales::number_format(accuracy = 1, suffix = "%"), breaks = targets) +
  scale_fill_manual(values = cols[c(2, 4, 6, 8, 10, 12)], 
                    na.value = "light grey",
                    labels = c(levels(df_plot$CoverageOneRound_factor), 
                               "Target not met within\nparameter space"))

p <- p + labs(title = "COVERAGE OF AT LEAST ONE SMC ROUND",  y = "TARGET  REDUCTION", x = expression(paste("BASELINE ANNUAL ", italic("Pf"), "PR"["2-10"]))) +
  guides(fill = guide_legend(title.position = "top", title.hjust = 0.5, title = "COVERAGE  OF AT  LEAST  ONE  SMC  ROUND", nrow = 2))

p2 <- p


### Generate plot for duration ###

p <- ggplot(df_plot, aes(x = annual_prev, y = Target, fill = Halflife_factor))

p <- p + geom_tile(colour = "white")

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
               plot.title = element_text(hjust = 0.5, face = "bold"),
               legend.key.width = unit(1, "cm"),
               legend.title = element_text(face = "bold", size = text_size),
               legend.title.align = 0.5,
               legend.position = "bottom")

p <- p + scale_y_continuous(labels = scales::number_format(accuracy = 1, suffix = "%"), breaks = targets)  +
  scale_fill_manual(values = cols[c(2, 6, 8, 12)], 
                    na.value = "light grey",
                    labels = c(levels(df_plot$Halflife_factor), 
                               "Target not met within\nparameter space"))

p <- p + labs(title = "CONCENTRATION HALF-LIFE",  y = "TARGET  REDUCTION", x = expression(paste("BASELINE ANNUAL ", italic("Pf"), "PR"["2-10"]))) +
  guides(fill = guide_legend(title.position = "top", title.hjust = 0.5, title = "CONCENTRATION HALF-LIFE  (DAYS)", nrow = 1))

p3 <- p


# Generate plot for Max Killing Rate

p <- ggplot(df_plot, aes(x = annual_prev, y = Target, fill = MaxKillingRate_factor))

p <- p + geom_tile(colour = "white")

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
               plot.title = element_text(hjust = 0.5, face = "bold"),
               legend.key.width = unit(1, "cm"),
               legend.title = element_text(face = "bold", size = text_size),
               legend.title.align = 0.5,
               legend.position = "bottom")

p <- p + scale_y_continuous(labels = scales::number_format(accuracy = 1, suffix = "%"), breaks = targets) +
  scale_fill_manual(values = cols[c(2, 8, 12)], 
                    na.value = "light grey",
                    labels = c(levels(df_plot$MaxKillingRate_factor), 
                               "Target not met within\nparameter space"))

p <- p + labs(title = "MAXIMUM EFFECT",  y = "TARGET  REDUCTION", x = expression(paste("BASELINE ANNUAL ", italic("Pf"), "PR"["2-10"]))) +
  guides(fill = guide_legend(title.position = "top", title.hjust = 0.5, title = "MAXIMUM  PARASITE KILLING  EFFECT", nrow = 1))

p4 <- p

p1 + p2 + p3 + p4 + plot_layout(nrow = 2, ncol = 2)

ggsave(filename = paste0("./analysisworkflow/analysis_scripts/iTPP3_Publication/Figures/iTPP3_cutoff_criteria_SUMMARY_3moLOWaccess.jpg"),
       plot = last_plot(),
       width = 9.1,
       height = 6,
       dpi = 400)
