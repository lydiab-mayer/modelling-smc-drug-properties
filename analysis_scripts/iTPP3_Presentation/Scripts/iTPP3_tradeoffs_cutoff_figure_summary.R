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
num_rounds <- c("iTPP3_tradeoffs_3rounds" = 3, "iTPP3_tradeoffs_4rounds" = 4, "iTPP3_tradeoffs_5rounds" = 5)

# !!! Insert your predicted parameter here. Note that this must match with one column name in post-processing files !!!
pred_list <- c("inc_red_int_Tot", "sev_red_int_Tot")

# !!! Identify desired targets for optimisation !!!
targets <- c(40, 50, 60, 70, 80, 90)

library(tidyr)
library(dplyr)
library(hetGP)
library(ggplot2)
library(patchwork)

user <- strsplit(getwd(), "/", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]][5]
GROUP_dr <- "/scicore/home/penny/GROUP/M3TPP/"
setwd(paste0("/scicore/home/penny/", user, "/M3TPP/"))

out <- matrix(nrow = 1, ncol = 5)
out[1, ] <- c("Scenario", "Target", "Coverage", "Halflife", "Efficacy")


# ----------------------------------------------------------
# Calculate cutoff criteria
# ----------------------------------------------------------

# for(exp in exp_list) {
#   
#   # Set up
#   if (!dir.exists(paste0("./Experiments/", exp, "/Outputs"))) dir.create(paste0("./Experiments/", exp, "/Outputs"))
#   # 
#   # load(paste0(GROUP_dr, exp, "/param_ranges.RData"))
#   # param_ranges_cont
#   
#   # Import settings
#   setting <- c()
#   
#   for(pred in pred_list) {
#     setting <- c(setting, Sys.glob(paste0(GROUP_dr, exp, "/gp/GP_grid_optimization/", pred, "/*")))
#   }
#   
#   setting_id <- sub(".rds", "", sub("opt_", "", basename(setting)))
#   
#   # Calculate cutoff criteria
#   for(i in 1:length(setting_id)){
#     print(paste0("Generating table for setting ", i , " of ", length(setting_id), " for experiment ", exp))
#   
#     # Set up table to store results
#     out_int <- matrix(nrow = 1 + length(targets), ncol = 5)
#     out_int[1, ] <- c("Scenario", "Target", "Coverage", "Halflife", "Efficacy")
#     out_int[2:(1 + length(targets)), 1] <- setting_id[i]
#     out_int[2:(1 + length(targets)), 2] <- targets
#     
#     # Read in optimisation result
#     df <- readRDS(setting[i])$scenario
#     df$CoverageAllRounds <- df$Coverage1*df$Coverage2^num_rounds[exp]
#     
#     # Extract minimum criteria for selected targets
#     for (i in 1:length(targets)){
#       out_int[i + 1, 3] <- min(df[df$mean >= targets[i], "CoverageAllRounds"])
#       out_int[i + 1, 4] <- min(df[df$mean >= targets[i], "Halflife"])
#       out_int[i + 1, 5] <- min(df[df$mean >= targets[i], "Efficacy"])
#     }
#     
#     # Store outputs
#     out <- rbind(out, out_int[2:(1 + length(targets)), ])
#   }
# 
#   # Save resulting matrix against future use
#   saveRDS(out, paste0("./analysisworkflow/analysis_scripts/iTPP3_Figures/", exp, "_tradeoffs_cutoff_criteria.rds"))
#   #out <- readRDS(paste0("./analysisworkflow/analysis_scripts/iTPP3_Figures/tradeoffs_cutoff_criteria.rds"))
# }
#  

 
# ----------------------------------------------------------
# Generate and format database of results across all experiments
# ----------------------------------------------------------

# Merge results across experiments

# !!! Note that this must be updated manually should the list of experiments be changed !!!
df3 <- readRDS(paste0("./analysisworkflow/analysis_scripts/iTPP3_Figures/iTPP3_tradeoffs_3rounds_tradeoffs_cutoff_criteria.rds"))
df4 <- readRDS(paste0("./analysisworkflow/analysis_scripts/iTPP3_Figures/iTPP3_tradeoffs_4rounds_tradeoffs_cutoff_criteria.rds"))
df5 <- readRDS(paste0("./analysisworkflow/analysis_scripts/iTPP3_Figures/iTPP3_tradeoffs_5rounds_tradeoffs_cutoff_criteria.rds"))

df <- rbind(df3, df4[-1, ], df5[-1, ])
remove(df3, df4, df5)

# Add column names
colnames(df) <- df[1, ]
df <- as.data.frame(df[-1, ], stringsAsFactors = FALSE)
  
# Replace infs with NAs
df <- data.frame(lapply(df, function(x) gsub("Inf", NA, x)), stringsAsFactors = FALSE)
  
# Format performance characteristics as numeric
df$Target <- as.numeric(df$Target)
df$Coverage <- as.numeric(df$Coverage)
df$Halflife <- as.numeric(df$Halflife)
df$Efficacy <- as.numeric(df$Efficacy)  
  
# Add additional columns for each scenario factor
df$Scenario <- sub("_red_int_Tot", "", sub("iTPP3tradeoffs", "", df$Scenario))
df <- df %>%
  separate(col = Scenario,
           into = c("Rounds", "Seasonality", "System", "EIR", "Agegroup", "Decay", "Access", "Timing", "Endpoint"),
           sep= "_")
  
# Merge in baseline prevalence
prev <- read.csv(paste0("./Experiments/iTPP3_tradeoffs_4rounds/Outputs/Prevalence_prior_to_intervention.csv"))
prev <- prev[, c("Seasonality", "EIR", "Access", "annual_prev")]
prev$annual_prev <- paste0(round(prev$annual_prev*100, 0), "%")
df <- merge(df, prev, by = c("Seasonality", "EIR", "Access"))
  
  
# Format text
df$Endpoint <- factor(df$Endpoint, levels = c("inc", "sev", "mor"))
df$Endpoint <- recode(df$Endpoint,
                      "inc" = "CLINICAL INCIDENCE REDUCTION",
                      "sev" = "SEVERE DISEASE REDUCTION",
                      "mor" = "MORTALITY REDUCTION")
  
df$Rounds <- factor(df$Rounds, levels = c("3rounds", "4rounds", "5rounds"))
df$Rounds <- recode(df$Rounds,
                     "3rounds" = "3 SMC ROUNDS",
                     "4rounds" = "4 SMC ROUNDS",
                     "5rounds" = "5 SMC ROUNDS")
  
df$Seasonality <- factor(df$Seasonality, levels = c("sharpseasonal", "wideseasonal"))
df$Seasonality <- recode(df$Seasonality,
                          "sharpseasonal" = "SHORT SEASON",
                          "wideseasonal" = "LONG SEASON")
  
df$Agegroup <- factor(df$Agegroup, levels = c("5", "10"))
df$Agegroup <- recode(df$Agegroup,
                       "5" = "CHILDREN 3 TO 59 MONTHS",
                       "10" = "CHILDREN 3 TO 119 MONTHS")
  
df$Access <- factor(df$Access, levels = c("0.04", "0.24"))
df$Access <- recode(df$Access,
                     "0.04" = "LOW ACCESS",
                     "0.24" = "HIGH ACCESS")
  
head(df)
  
#saveRDS(df, paste0("./analysisworkflow/analysis_scripts/iTPP3_Presentation/Figures/iTPP3_tradeoffs_cutoff_criteria.rds"))
  

# ----------------------------------------------------------
# Plot results - executive summary
# ----------------------------------------------------------

text_size <- 6

df_plot <- df

# Calculate max across all scenarios and all outcomes
df_plot <- df_plot %>%
  group_by(Access, Seasonality, Target, annual_prev) %>%
  summarise(sum_na = sum(is.na(Coverage)),
            Coverage = max(Coverage, na.rm = FALSE),
            Halflife = max(Halflife, na.rm = FALSE),
            Efficacy = max(Efficacy, na.rm = FALSE))

#df_plot$sum_na <- 1 - ifelse(df_plot$sum_na == 0, 0, df_plot$sum_na / sum(df_plot$sum_na))


df_plot <- df_plot[df_plot$Seasonality == "SHORT SEASON" & df_plot$Access == "HIGH ACCESS", ]
df_plot <- df_plot[df_plot$Target != 90, ]



# Generate plot for Coverage

p <- ggplot(df_plot, aes(x = annual_prev, y = Target, fill = Coverage))

p <- p + geom_tile(colour = "white")

#p <- p + facet_wrap(.~ Seasonality + Access, scales = "free")

p <- p + theme(panel.border = element_blank(), 
               panel.background = element_blank(),
               panel.grid = element_blank(),
               text = element_text(family = "Arial", size = text_size),
               strip.background = element_blank(),
               strip.text = element_text(face = "bold"),
               axis.line = element_blank(),
               axis.ticks = element_blank(),
               axis.text.x = element_text(margin = margin(t = 0)),
               axis.text.y = element_text(margin = margin(r = 0)),
               axis.title = element_text(face = "bold"),
               plot.title = element_text(hjust = 0.5, face = "bold"),
               legend.key.size = unit(0.25, "cm"),
               legend.text = element_text(size = text_size),
               legend.title = element_text(face = "bold", size = text_size),
               legend.position = "right",
               legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))

p <- p + scale_y_continuous(labels = scales::number_format(accuracy = 1, suffix = "%"), breaks = targets) +
  scale_fill_stepsn(colours = c("#1b4d79", "#5880b1", "#f9a24b", "#f9d48e", "#ffedcb"),
                    na.value = "light grey",
                    breaks = c(0.6, 0.7, 0.8, 0.9, 1.0),
                    labels = function(x) paste0(x*100, "%"),
                    limits = c(0.6, 1.0),
                    show.limits = TRUE)

p <- p + labs(title = "MINIMUM COVERAGE CRITERIA", y = "TARGET REDUCTION", x = expression(paste(bold("BASELINE ANNUAL "), bolditalic("Pf"), bold("PR"["2-10"])))) +
  guides(fill = guide_colourbar(title.position="top", title.hjust = 0.5, title = "COVERAGE OF AT LEAST\nONE SMC ROUND"))

p1 <- p


# Generate plot for duration

p <- ggplot(df_plot, aes(x = annual_prev, y = Target, fill = Halflife))

p <- p + geom_tile(colour = "white")

#p <- p + facet_wrap(.~ Seasonality + Access, scales = "free")

p <- p + theme(panel.border = element_blank(), 
               panel.background = element_blank(),
               panel.grid = element_blank(),
               text = element_text(family = "Arial", size = text_size),
               strip.background = element_blank(),
               strip.text = element_text(face = "bold"),
               axis.line = element_blank(),
               axis.ticks = element_blank(),
               axis.text.x = element_text(margin = margin(t = 0)),
               axis.text.y = element_text(margin = margin(r = 0)),
               axis.title = element_text(face = "bold"),
               plot.title = element_text(hjust = 0.5, face = "bold"),
               legend.key.size = unit(0.25, "cm"),
               legend.text = element_text(size = text_size),
               legend.title = element_text(face = "bold", size = text_size),
               legend.position = "right",
               legend.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))

p <- p + scale_y_continuous(labels = scales::number_format(accuracy = 1, suffix = "%"), breaks = targets) +
  scale_fill_stepsn(colours = c("#1b4d79", "#5880b1", "#f9a24b", "#f9d48e", "#ffedcb"),
                    na.value = "light grey",
                    breaks = c(10, 21, 28, 35, 42, 49, 60),
                    #labels = function(x) paste0(x, "%"),
                    limits = c(10, 60),
                    show.limits = TRUE)

p <- p + labs(title = "MINIMUM DURATION CRITERIA",  y = "TARGET REDUCTION", x = expression(paste(bold("BASELINE ANNUAL "), bolditalic("Pf"), bold("PR"["2-10"])))) +
  guides(fill = guide_colourbar(title.position = "top", title.hjust = 0.5, title = "DURATION OF PROTECTION"))

p2 <- p
  
  
# Generate plot for Efficacy
  
p <- ggplot(df_plot, aes(x = annual_prev, y = Target, fill = Efficacy))
  
p <- p + geom_tile(colour = "white")
  
#p <- p + facet_wrap(.~ Seasonality + Access, scales = "free")
  
p <- p + theme(panel.border = element_blank(), 
               panel.background = element_blank(),
               panel.grid = element_blank(),
               text = element_text(family = "Arial", size = text_size),
               strip.background = element_blank(),
               strip.text = element_text(face = "bold"),
               axis.line = element_blank(),
               axis.ticks = element_blank(),
               axis.text.x = element_text(margin = margin(t = 0)),
               axis.text.y = element_text(margin = margin(r = 0)),
               axis.title = element_text(face = "bold"),
               plot.title = element_text(hjust = 0.5, face = "bold"),
               legend.key.size = unit(0.25, "cm"),
               legend.text = element_text(size = text_size),
               legend.title = element_text(face = "bold", size = text_size),
               legend.position = "right",
               legend.margin = margin(t = 0, r = 30, b = 0, l = 0, unit = "pt"))
  
p <- p + scale_y_continuous(labels = scales::number_format(accuracy = 1, suffix = "%"), breaks = targets) +
  scale_fill_stepsn(colours = c("#1b4d79", "#5880b1", "#f9a24b", "#f9d48e", "#ffedcb"),
                    na.value = "light grey",
                    breaks = c(0.8, 0.85, 0.9, 0.95, 1.0),
                    labels = function(x) paste0(x*100, "%"),
                    limits = c(0.8, 1.0),
                    show.limits = TRUE)
  
p <- p + labs(title = "MINIMUM EFFICACY CRITERIA", y = "TARGET REDUCTION", x = expression(paste(bold("BASELINE ANNUAL "), bolditalic("Pf"), bold("PR"["2-10"])))) +
  guides(fill = guide_colourbar(title.position="top", title.hjust = 0.5, title = "INITIAL EFFICACY"))

p3 <- p
  
p1 + p2 + p3 + plot_layout(nrow = 3, ncol = 1)
  
ggsave(filename = paste0("./analysisworkflow/analysis_scripts/iTPP3_Presentation/Figures/iTPP3_tradeoffs_cutoff_criteria_SUMMARY.jpg"),
       plot = last_plot(),
       width = 5,
       height = 4,
       dpi = 400)
