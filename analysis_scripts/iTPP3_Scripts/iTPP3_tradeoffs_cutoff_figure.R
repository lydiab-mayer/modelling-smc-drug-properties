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
exp_list <-c("iTPP3_tradeoffs_5rounds")
num_rounds <- c("iTPP3_tradeoffs_3rounds" = 3, "iTPP3_tradeoffs_4rounds" = 4, "iTPP3_tradeoffs_5rounds" = 5)

# !!! Insert your predicted parameter here. Note that this must match with one column name in post-processing files !!!
pred_list <- c("inc_red_int_Tot", "sev_red_int_Tot")

# !!! Identify desired targets for optimisation !!!
targets <- c(40, 50, 60, 70, 80)

library(tidyr)
library(dplyr)
library(hetGP)

user <- strsplit(getwd(), "/", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]][5]
GROUP_dr <- "/scicore/home/penny/GROUP/M3TPP/"
setwd(paste0("/scicore/home/penny/", user, "/M3TPP/"))

out <- matrix(nrow = 1, ncol = 5)
out[1, ] <- c("Scenario", "Target", "Coverage", "Halflife", "Efficacy")


# ----------------------------------------------------------
# Calculate cutoff criteria
# ----------------------------------------------------------

for(exp in exp_list) {
  
  # Set up
  if (!dir.exists(paste0("./Experiments/", exp, "/Outputs"))) dir.create(paste0("./Experiments/", exp, "/Outputs"))
  # 
  # load(paste0(GROUP_dr, exp, "/param_ranges.RData"))
  # param_ranges_cont
  
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
    out_int <- matrix(nrow = 1 + length(targets), ncol = 5)
    out_int[1, ] <- c("Scenario", "Target", "Coverage", "Halflife", "Efficacy")
    out_int[2:(1 + length(targets)), 1] <- setting_id[i]
    out_int[2:(1 + length(targets)), 2] <- targets
    
    # Read in optimisation result
    df <- readRDS(setting[i])$scenario
    df$CoverageAllRounds <- df$Coverage1*df$Coverage2^num_rounds[exp]
    
    # Extract minimum criteria for selected targets
    for (i in 1:length(targets)){
      out_int[i + 1, 3] <- min(df[df$mean >= targets[i], "CoverageAllRounds"])
      out_int[i + 1, 4] <- min(df[df$mean >= targets[i], "Halflife"])
      out_int[i + 1, 5] <- min(df[df$mean >= targets[i], "Efficacy"])
    }
    
    # Store outputs
    out <- rbind(out, out_int[2:(1 + length(targets)), ])
  }

  # Save resulting matrix against future use
  saveRDS(out, paste0("./analysisworkflow/analysis_scripts/iTPP3_Figures/", exp, "_tradeoffs_cutoff_criteria.rds"))
  #out <- readRDS(paste0("./analysisworkflow/analysis_scripts/iTPP3_Figures/tradeoffs_cutoff_criteria.rds"))
  
  
  # ----------------------------------------------------------
  # Format cutoff criteria outputs
  # ----------------------------------------------------------
  
  # Add column names
  colnames(out) <- out[1, ]
  out <- as.data.frame(out[-1, ], stringsAsFactors = FALSE)
  
  # Replace infs with NAs
  out <- data.frame(lapply(out, function(x) gsub("Inf", NA, x)), stringsAsFactors = FALSE)
  
  
  # Format performance characteristics as numeric
  out$Target <- as.numeric(out$Target)
  out$Coverage <- as.numeric(out$Coverage)
  out$Halflife <- as.numeric(out$Halflife)
  out$Efficacy <- as.numeric(out$Efficacy)  
  
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
                        "inc" = "CLINICAL INCIDENCE REDUCTION",
                        "sev" = "SEVERE DISEASE REDUCTION",
                        "mor" = "MORTALITY REDUCTION")
  
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
  
  
  # # Tidy up table and order rows
  # 
  # names(out)[names(out) == "annual_prev"] <- "PfPR2_10"
  # 
  # out$indexOutcome <- sub(" \\(.*", "", out[, "Target"])
  # out$indexOutcome <- factor(out$indexOutcome, levels = c("High", "Moderate", "Low"))
  # 
  # out <- out[order(out$Seasonality, out$Access, out$Rounds, out$Agegroup, out$PfPR2_10, out$Endpoint, out$indexOutcome), ]
  # 
  # out <- out[, c("Scenario", "PfPR2_10", "Endpoint", "Target", "Coverage", "Halflife", "Efficacy")]
  # 
  #
  # # Remove duplicated text
  # out <- out %>%
  #   mutate(indexScenario = Scenario == lag(Scenario, 1),
  #          indexPfPR = PfPR2_10 == lag(PfPR2_10, 1))
  # 
  # out[is.na(out)] <- FALSE
  # 
  # out[out$indexScenario, ]$Scenario <- NA
  # out[out$indexPfPR, ]$PfPR2_10 <- NA
  #
  #out <- out[, !(names(out) %in% c("indexScenario", "indexPfPR", "indexOutcome"))]
  #
  #head(out)
  
  
  # ----------------------------------------------------------
  # Plot results
  # ----------------------------------------------------------
  
  text_size <- 10
  
  df_plot <- out[out$Agegroup == "children aged between 3 and 59 months", ]
  df_plot <- df_plot[df_plot$Access == "high access to care", ]
  df_plot <- df_plot[df_plot$Seasonality == "short seasonal profile", ]
  
  # Generate plot for Halflife
  
  p <- ggplot(df_plot, aes(x = annual_prev, y = Target, fill = Halflife))
  
  p <- p + geom_tile(colour = "white")
  
  p <- p + facet_wrap(.~ Endpoint, nrow = 2, ncol = 1)
  
  p <- p + theme(panel.border = element_blank(), 
                 panel.background = element_blank(),
                 panel.grid = element_blank(),
                 text = element_text(family = "Courier", size = text_size),
                 strip.background = element_blank(),
                 strip.text = element_text(face = "bold"),
                 axis.line = element_blank(),
                 axis.ticks = element_blank(),
                 axis.text.x = element_text(margin = margin(t = 0)),
                 axis.text.y = element_text(margin = margin(r = 0)),
                 axis.title = element_text(face = "bold"),
                 plot.title = element_text(hjust = 0.5, face = "bold"),
                 legend.key.width = unit(1, "cm"),
                 legend.title = element_text(face = "bold", size = text_size),
                 legend.position = "bottom")
  
  p <- p + scale_y_continuous(labels = scales::number_format(accuracy = 1, suffix = "%"), breaks = targets) +
    scale_fill_stepsn(colours = c("#1b4d79", "#5880b1", "#f9a24b", "#f9d48e", "#ffedcb"),
                      na.value = "light grey",
                      breaks = c(10, 20, 30, 40, 50, 60),
                      #labels = function(x) paste0(x, "%"),
                      limits = c(10, 60),
                      show.limits = TRUE)
  
  p <- p + labs(title = "MINIMUM DURATION CRITERIA",  y = "TARGET REDUCTION", x = expression(paste(bold("BASELINE ANNUAL "), bolditalic("Pf"), bold("PR"["2-10"])))) +
    guides(fill = guide_colourbar(title.position = "top", title.hjust = 0.5, title = "DURATION OF PROTECTION"))
  
  p1 <- p
  
  
  # Generate plot for Coverage
  
  p <- ggplot(df_plot, aes(x = annual_prev, y = Target, fill = Coverage))
  
  p <- p + geom_tile(colour = "white")
  
  p <- p + facet_wrap(.~ Endpoint, nrow = 2, ncol = 1)
  
  p <- p + theme(panel.border = element_blank(), 
                 panel.background = element_blank(),
                 panel.grid = element_blank(),
                 text = element_text(family = "Courier", size = text_size),
                 strip.background = element_blank(),
                 strip.text = element_text(face = "bold"),
                 axis.line = element_blank(),
                 axis.ticks = element_blank(),
                 axis.text.x = element_text(margin = margin(t = 0)),
                 axis.text.y = element_text(margin = margin(r = 0)),
                 axis.title = element_text(face = "bold"),
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
  
  p <- p + labs(title = "MINIMUM COVERAGE CRITERIA", y = "", x = expression(paste(bold("BASELINE ANNUAL "), bolditalic("Pf"), bold("PR"["2-10"])))) +
    guides(fill = guide_colourbar(title.position="top", title.hjust = 0.5, title = "COVERAGE OF ALL SMC ROUNDS"))
  
  p2 <- p
  
  
  # Generate plot for Efficacy
  
  p <- ggplot(df_plot, aes(x = annual_prev, y = Target, fill = Efficacy))
  
  p <- p + geom_tile(colour = "white")
  
  p <- p + facet_wrap(.~ Endpoint, nrow = 2, ncol = 1)
  
  p <- p + theme(panel.border = element_blank(), 
                 panel.background = element_blank(),
                 panel.grid = element_blank(),
                 text = element_text(family = "Courier", size = text_size),
                 strip.background = element_blank(),
                 strip.text = element_text(face = "bold"),
                 axis.line = element_blank(),
                 axis.ticks = element_blank(),
                 axis.text.x = element_text(margin = margin(t = 0)),
                 axis.text.y = element_text(margin = margin(r = 0)),
                 axis.title = element_text(face = "bold"),
                 plot.title = element_text(hjust = 0.5, face = "bold"),
                 legend.key.width = unit(1, "cm"),
                 legend.title = element_text(face = "bold", size = text_size),
                 legend.position = "bottom")
  
  p <- p + scale_y_continuous(labels = scales::number_format(accuracy = 1, suffix = "%"), breaks = targets) +
    scale_fill_stepsn(colours = c("#1b4d79", "#5880b1", "#f9a24b", "#f9d48e", "#ffedcb"),
                      na.value = "light grey",
                      breaks = c(0.8, 0.85, 0.9, 0.95, 1.0),
                      labels = function(x) paste0(x*100, "%"),
                      limits = c(0.8, 1.0),
                      show.limits = TRUE)
  
  p <- p + labs(title = "MINIMUM EFFICACY CRITERIA", y = "", x = expression(paste(bold("BASELINE ANNUAL "), bolditalic("Pf"), bold("PR"["2-10"])))) +
    guides(fill = guide_colourbar(title.position="top", title.hjust = 0.5, title = "INITIAL EFFICACY"))
  
  p3 <- p
  
  p1 + p2 + p3 + plot_layout(nrow = 1, ncol = 3)
  
  ggsave(filename = paste0("./analysisworkflow/analysis_scripts/iTPP3_Figures/", exp, "_iTPP3_tradeoffs_cutoff_criteria.jpg"),
         plot = last_plot(),
         width = 9,
         height = 5.6,
         dpi = 200)

}