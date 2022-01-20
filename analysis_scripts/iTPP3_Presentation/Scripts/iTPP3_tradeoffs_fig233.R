############################################################
#
# Generates figure 2.3.3
# Note: This script depends on outputs of the script 4_sensitivityanalysis_workflow.R
#
# Written by Lydia Braunack-Mayer
############################################################

# ----------------------------------------------------------
# Setup
# ----------------------------------------------------------

rm(list = ls())

# !!! Insert your experiment name here as a string, e.g. "MyExperiment" !!!
exp <- c("iTPP3_tradeoffs_4rounds")


library(ggplot2)
library(dplyr)
library(tidyr)
library(wesanderson)
library(scales)
library(patchwork)


user <- strsplit(getwd(), "/", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]][5]
GROUP_dr <- "/scicore/home/penny/GROUP/M3TPP/"
setwd(paste0("/scicore/home/penny/", user, "/M3TPP"))


# ----------------------------------------------------------
# GENERATE RESULTS FOR SENSITIVITY ANALYSIS
# ----------------------------------------------------------
  
  if (!dir.exists(paste0("./Experiments/", exp, "/Outputs"))) dir.create(paste0("./Experiments/", exp, "/Outputs"))
  
  load(paste0(GROUP_dr, exp, "/param_ranges.RData"))
  param_ranges_cont
  
  
  # ----------------------------------------------------------
  # Define data to plot
  # ----------------------------------------------------------
  
  # Import settings
  setting <- Sys.glob(paste0(GROUP_dr, exp, "/gp/trained/sensitivity_clinicaltranslation/*"))
  setting <- setting[setting != paste0(GROUP_dr, exp, "/gp/trained/sensitivity_clinicaltranslation/err")]
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
             into = c("Experiment", "Seasonality", "System", "EIR", "Agegroup", "Decay", "Access", "Timing", "Outcome", "temp1", "temp2", "temp3"),
             sep = "_",
             remove = FALSE)
  df <- df[, !(names(df) %in% c("temp1", "temp2", "temp3"))]
    
    
  # Import aggregated impact for each setting
  
  setting_id <- unique(sub("_Apr.*", "_Apr", sub("_May.*", "_May", setting_id)))
  df_impact <- data.frame()
    
  for (i in 1:length(setting_id)) {
      
    temp <- read.table(paste0(GROUP_dr, exp, "/postprocessing/agg_",setting_id[i], ".txt"), header = TRUE, sep = "")
    temp <- temp %>%
      group_by(Seasonality, Biting_pattern, EIR, MaxAge, Decay, Access, Timing) %>%
      summarise(inc = median(inc_red_int_Tot),
                sev = median(sev_red_int_Tot),
                mor = median(mor_red_int_Tot),
                prev = median(prev_red_int_Tot))
    temp$scenario <- setting_id[i]
      
    df_impact <- rbind(df_impact, as.data.frame(temp))
      
  }
    
  df_impact <- df_impact[, c("inc", "sev", "mor", "prev", "scenario")] %>%
    separate(col = scenario, 
             into = c("Experiment", "Seasonality", "System", "EIR", "Agegroup", "Decay", "Access", "Timing"),
             sep = "_",
             remove = FALSE)
    
  df_impact <- df_impact %>%
    pivot_longer(cols = c(inc, sev, mor, prev), names_to = "Outcome", values_to = "Median_Reduction")
    
  
  # Scale total effects by median impact for each setting
  
  df <- merge(df[, !(names(df) == "scenario")],
              df_impact, 
              by = c("Experiment", "Seasonality", "System", "EIR", "Agegroup", "Decay", "Access", "Timing", "Outcome"))
    
  df$T_eff_scaled <- df$T_eff * df$Median_Reduction
    
    
  # ----------------------------------------------------------
  # Import baseline prevalence
  # ----------------------------------------------------------
    
  # Import from csv
  prev <- read.csv(paste0("./Experiments/", exp, "/Outputs/Prevalence_prior_to_intervention.csv"))
  prev <- prev[, c("Seasonality", "EIR", "Access", "annual_prev")]
  prev$annual_prev <- paste0(round(prev$annual_prev*100, 0), "%")
  
  # Merge into data
  df <- merge(df, prev, by = c("Seasonality", "EIR", "Access"))
    
    
  # ----------------------------------------------------------
  # Retain only data relevant for clinical translation
  # ----------------------------------------------------------
  
  df_plot <- df
  df_plot <- df_plot[df_plot$parameter %in% c("Halflife", "Efficacy"), ]
  df_plot <- df_plot[df_plot$Access == 0.24, ]
  df_plot <- df_plot[df_plot$Agegroup == 5, ]
  df_plot <- df_plot[df_plot$Outcome %in% c("inc", "sev"), ]
  
  
  # ----------------------------------------------------------
  # Format data for plotting - MUST BE ADJUSTED FOR EACH EXPERIMENT
  # ----------------------------------------------------------
    
  df_plot$parameter <- as.factor(df_plot$parameter); df_plot$scenario <- as.factor(df_plot$scenario)
  # df <- df %>%
  #   separate(col = scenario, 
  #            into = c("Experiment", "Seasonality", "System", "EIR", "Agegroup", "Decay", "Access", "Timing"),
  #            sep = "_",
  #            remove = FALSE)
  
  df_plot$Seasonality <- ifelse(df_plot$Seasonality == "sharpseasonal", "SHORT SEASON", "LONG SEASON")
  df_plot$Seasonality <- factor(df_plot$Seasonality, levels = c("SHORT SEASON", "LONG SEASON"))
  
  df_plot$Access <- ifelse(df_plot$Access == 0.04, "LOW ACCESS", "HIGH ACCESS")
  df_plot$Access <- factor(df_plot$Access, levels = c("LOW ACCESS", "HIGH ACCESS"))
  
  df_plot$Agegroup <- factor(df_plot$Agegroup, levels = c("5", "10"))
  df_plot$Agegroup <- recode(df_plot$Agegroup,
                        "5" = "CHILDREN 3 TO 59 MONTHS",
                        "10" = "CHILDREN 3 TO 119 MONTHS")
  
  df_plot$parameter <- factor(df_plot$parameter, levels = c("Halflife", "Efficacy"))
  df_plot$parameter <- recode(df_plot$parameter,
                         "Efficacy" = "Initial efficacy [80% - 100%]",
                         "Halflife" = "Duration of protection [10 - 60 days]")

  df_plot$annual_prev <- factor(df_plot$annual_prev)
  
  df_plot$Outcome <- factor(df_plot$Outcome, levels = c("inc", "sev"))
  df_plot$Outcome <- recode(df_plot$Outcome, "inc" = "CLINICAL INCIDENCE", "sev" = "SEVERE DISEASE")
    
    
  # ----------------------------------------------------------
  # Define plot settings
  # ----------------------------------------------------------
    
  # Define name of your plot
  plot_name <- "FIG233"
  
  # Define colours
  cols <- wes_palette("Royal1", n = 4)[c(1, 4)]
  text_cols <- c("#323d42", "#7f4a1b")

  
  # ----------------------------------------------------------
  # Generate plot
  # ----------------------------------------------------------
    
  p <- ggplot(df_plot[df_plot$Seasonality == "SHORT SEASON", ], 
              aes(x = annual_prev, y = T_eff_scaled, fill = parameter, label = paste0(round(T_eff*100, 0), "%")))
    
  p <- p + geom_bar(position = "stack", stat = "identity", colour = "white")
  
  p <- p + geom_text(aes(colour = parameter), 
                     position = position_stack(vjust = 0.5),
                     family = "Arial",
                     fontface = "bold",
                     show.legend = FALSE)
    
  p <- p + facet_wrap(. ~ Outcome, scales = "free", nrow = 2, ncol = 1)
    
  p <- p + theme(panel.border = element_blank(), 
                 panel.background = element_blank(),
                 panel.grid = element_blank(),
                 text = element_text(family = "Arial", size = 10),
                 strip.background = element_blank(),
                 axis.line = element_blank(),
                 axis.ticks = element_blank(),
                 axis.title.x = element_text(face = "bold", margin = margin(t = 10)),
                 axis.title.y = element_text(face = "bold", margin = margin(r = 10)),
                 strip.text = element_text(face = "bold"),
                 legend.key = element_blank(),
                 legend.title = element_text(face = "bold"),
                 legend.text = element_text(size = 8),
                 legend.position = "bottom",
                 legend.margin = margin(t = -5))
  
  p <- p + scale_fill_manual(values = cols) + 
    scale_colour_manual(values = text_cols) +
    scale_y_continuous(breaks = seq(0, 100, 10), labels = paste0(seq(0, 100, 10), "%"))
  
  p <- p + labs(x = expression(paste(bold("BASELINE ANNUAL "), bolditalic("Pf"), bold("PR"["2-10"]))),
                y = "MEDIAN REDUCTION (%)",
                fill = "") +
    guides(fill = guide_legend(nrow = 2))
    
  p
    
  ggsave(filename = paste0("./analysisworkflow/analysis_scripts/iTPP3_Presentation/Figures/iTPP3_", plot_name, ".jpg"),
         plot = last_plot(),
         width = 5,
         height = 4,
         dpi = 400)
  