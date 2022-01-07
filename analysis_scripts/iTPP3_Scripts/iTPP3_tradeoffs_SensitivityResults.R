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
exp_list <- c("iTPP3_tradeoffs_3rounds", "iTPP3_tradeoffs_4rounds", "iTPP3_tradeoffs_5rounds")

# !!! Insert your predicted parameters here. Note that this must match with one column name in post-processing files !!!
pred_list <- c("inc_red_int_Tot", "sev_red_int_Tot")

library(ggplot2)
library(dplyr)
library(tidyr)
library(wesanderson)

user <- strsplit(getwd(), "/", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]][5]
GROUP_dr <- "/scicore/home/penny/GROUP/M3TPP/"
setwd(paste0("/scicore/home/penny/", user, "/M3TPP"))

for(exp in exp_list) {
  
  if (!dir.exists(paste0("./Experiments/", exp, "/Outputs"))) dir.create(paste0("./Experiments/", exp, "/Outputs"))
  
  load(paste0(GROUP_dr, exp, "/param_ranges.RData"))
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
                sev = median(sev_red_int_Tot))
    temp$scenario <- setting_id[i]
      
    df_impact <- rbind(df_impact, as.data.frame(temp))
      
  }
    
  df_impact <- df_impact[, c("inc", "sev", "scenario")] %>%
    separate(col = scenario, 
             into = c("Experiment", "Seasonality", "System", "EIR", "Agegroup", "Decay", "Access", "Timing"),
             sep = "_",
             remove = FALSE)
    
  df_impact <- df_impact %>%
    pivot_longer(cols = c(inc, sev), names_to = "Outcome", values_to = "Median_Reduction")
    
  
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
  # Format data for plotting - MUST BE ADJUSTED FOR EACH EXPERIMENT
  # ----------------------------------------------------------
    
  df$parameter <- as.factor(df$parameter); df$scenario <- as.factor(df$scenario)
  df <- df %>%
    separate(col = scenario, 
             into = c("Experiment", "Seasonality", "System", "EIR", "Agegroup", "Decay", "Access", "Timing"),
             sep = "_",
             remove = FALSE)
  
  df$Seasonality <- ifelse(df$Seasonality == "sharpseasonal", "SHORT SEASON", "LONG SEASON")
  df$Seasonality <- factor(df$Seasonality, levels = c("SHORT SEASON", "LONG SEASON"))
  
  df$Access <- ifelse(df$Access == 0.04, "LOW ACCESS", "HIGH ACCESS")
  df$Access <- factor(df$Access, levels = c("LOW ACCESS", "HIGH ACCESS"))
  
  df$Setting <- paste0(df$Seasonality, "\n", df$Access)
  df$Setting <- factor(df$Setting, levels = c("LONG SEASON\nLOW ACCESS", 
                                              "SHORT SEASON\nLOW ACCESS", 
                                              "LONG SEASON\nHIGH ACCESS",
                                              "SHORT SEASON\nHIGH ACCESS"))
    
  df$Agegroup <- factor(df$Agegroup, levels = c("5", "10"))
  df$Agegroup <- recode(df$Agegroup,
                        "5" = "CHILDREN 3 TO 59 MONTHS",
                        "10" = "CHILDREN 3 TO 119 MONTHS")
  
  df$parameter <- factor(df$parameter, levels = c("Coverage1", "Halflife", "Coverage2", "Efficacy"))
  df$parameter <- recode(df$parameter,
                         "Coverage1" = "SMC program coverage [70% - 100%]",
                         "Coverage2" = "SMC round coverage [70% - 100%]",
                         "Efficacy" = "Initial efficacy [80% - 100%]",
                         "Halflife" = "Duration of protection [10 - 60 days]")
  
  df$Outcome <- factor(df$Outcome, levels = c("inc", "sev"))
  df$Outcome <- recode(df$Outcome,
                       "inc" = "CLINICAL INCIDENCE",
                       "sev" = "SEVERE DISEASE")
    
  df$annual_prev <- factor(df$annual_prev)
    
    
  # ----------------------------------------------------------
  # Define plot settings
  # ----------------------------------------------------------
    
  # Define colours
  cols <- wes_palette("Royal1", n = 4)[c(2, 1, 3, 4)]
  text_cols <- c("#5f1909", "#323d42", "#827d55", "#7f4a1b")
  
  # Define y-axis labels
  y_axis <- c("MEDIAN INCIDENCE REDUCTION", "MEDIAN SEVERE DISEASE REDUCTION")
  names(y_axis) <- pred_list
    
  
  # ----------------------------------------------------------
  # Generate plot - figure 3.3.2
  # ----------------------------------------------------------
    
  p <- ggplot(df[(df$Access == "LOW ACCESS" & df$Seasonality == "LONG SEASON") & df$Agegroup == "CHILDREN 3 TO 59 MONTHS", ], 
              aes(x = annual_prev, y = T_eff_scaled, fill = parameter, label = paste0(round(T_eff*100, 0), "%")))
    
  p <- p + geom_bar(position = "stack", stat = "identity", colour = "white")
  
  p <- p + geom_text(aes(colour = parameter), 
                     position = position_stack(vjust = 0.5),
                     family = "Courier",
                     fontface = "bold",
                     show.legend = FALSE)
    
  p <- p + facet_wrap(.~ Outcome, scales = "free_x")
    
  p <- p + theme(panel.border = element_blank(), 
                 panel.background = element_blank(),
                 panel.grid = element_blank(),
                 text = element_text(family = "Courier", size = 16),
                 strip.background = element_blank(),
                 axis.line = element_blank(),
                 axis.ticks = element_blank(),
                 axis.title.x = element_text(face = "bold", margin = margin(t = 10)),
                 axis.title.y = element_text(face = "bold", margin = margin(r = 10)),
                 strip.text = element_text(face = "bold"),
                 legend.key = element_blank(),
                 legend.title = element_text(face = "bold"),
                 legend.position = "bottom")
  
  p <- p + scale_fill_manual(values = cols) + 
    scale_colour_manual(values = text_cols) +
    scale_y_continuous(breaks = seq(0, 100, 10), labels = paste0(seq(0, 100, 10), "%"))
  
  p <- p + labs(x = expression(paste(bold("BASELINE ANNUAL "), bolditalic("Pf"), bold("PR"["2-10"]))),
                y = "MEDIAN REDUCTION (%)",
                fill = "")
    
  p <- p + guides(fill = guide_legend(nrow = 2)) 
    
  p
    
  ggsave(filename = paste0("./analysisworkflow/analysis_scripts/iTPP3_Figures/iTPP3_", exp, "_FIG332.jpg"),
         plot = last_plot(),
         width = 9.1,
         height = 4,
         dpi = 200)
    
    
  # ----------------------------------------------------------
  # Generate table - table A.2.1
  # ----------------------------------------------------------
    
  tab <- df %>%
    group_by(parameter, Outcome) %>%
    summarise(max = max(S_eff), min = min(S_eff))
    
  tab$minmax <- paste0(round(tab$min*100, 0), "% to ", round(tab$max*100, 0), "%")
  names(tab) <- c("Key performance property", "Max", "Min", "Range of attributable outcome variation")
  write.csv(tab[, c(1, 2, 5)], file = paste0("./analysisworkflow/analysis_scripts/iTPP3_Figures/iTPP3_", exp, "_TABA21.csv"))
    
    
  # ----------------------------------------------------------
  # Generate plot - figure 6.3.1
  # ----------------------------------------------------------
    
  p <- ggplot(df[df$Access == "LOW ACCESS" & df$Seasonality == "LONG SEASON", ], 
              aes(x = annual_prev, y = T_eff_scaled, fill = parameter, label = paste0(round(T_eff*100, 0), "%")))
    
  p <- p + geom_bar(position = "stack", stat = "identity", colour = "white")
    
  p <- p + geom_text(aes(colour = parameter), position = position_stack(vjust = 0.5), family = "Courier", fontface = "bold",
                     show.legend = FALSE)
    
  p <- p + facet_grid(Agegroup ~ Outcome, scales = "free_x")
  
  p <- p + theme(panel.border = element_blank(), 
                 panel.background = element_blank(),
                 panel.grid = element_blank(),
                 text = element_text(family = "Courier", size = 16),
                 strip.background = element_blank(),
                 axis.line = element_blank(),
                 axis.ticks = element_blank(),
                 axis.title.x = element_text(face = "bold", margin = margin(t = 10)),
                 axis.title.y = element_text(face = "bold", margin = margin(r = 10)),
                 strip.text = element_text(face = "bold"),
                 legend.key = element_blank(),
                 legend.title = element_text(face = "bold"),
                 legend.position = "bottom")
  
  p <- p + scale_fill_manual(values = cols) + 
    scale_colour_manual(values = text_cols) +
    scale_y_continuous(breaks = seq(0, 100, 10), labels = paste0(seq(0, 100, 10), "%"))
    
  p <- p + labs(x = expression(paste(bold("BASELINE ANNUAL "), bolditalic("Pf"), bold("PR"["2-10"]))),
                y = "MEDIAN REDUCTION (%)",
                fill = "")
    
  p <- p + guides(fill = guide_legend(nrow = 2)) 
    
  p
    
  ggsave(filename = paste0("./analysisworkflow/analysis_scripts/iTPP3_Figures/iTPP3_", exp,"_FIG631.jpg"),
         plot = last_plot(),
         width = 9.1,
         height = 8,
         dpi = 200)
    
    
    
  # ----------------------------------------------------------
  # Generate plot - figure 7.3.1
  # ----------------------------------------------------------
    
  p <- ggplot(df[df$Agegroup == "CHILDREN 3 TO 59 MONTHS" & df$EIR == 32, ], 
              aes(x = Setting, y = T_eff_scaled, fill = parameter, label = paste0(round(T_eff*100, 0), "%")))
    
  p <- p + geom_bar(position = "stack", stat = "identity", colour = "white")
  
  p <- p + geom_text(aes(colour = parameter), position = position_stack(vjust = 0.5), family = "Courier", fontface = "bold",
                     show.legend = FALSE)
    
  p <- p + facet_wrap(.~ Outcome, scales = "free_x")
    
  p <- p + theme(panel.border = element_blank(), 
                 panel.background = element_blank(),
                 panel.grid = element_blank(),
                 text = element_text(family = "Courier", size = 14),
                 strip.background = element_blank(),
                 axis.line = element_blank(),
                 axis.ticks = element_blank(),
                 axis.text.x = element_text(size = 9),
                 axis.title.x = element_text(face = "bold", margin = margin(t = 10)),
                 axis.title.y = element_text(face = "bold", margin = margin(r = 10)),
                 strip.text = element_text(face = "bold"),
                 legend.key = element_blank(),
                 legend.title = element_text(face = "bold"),
                 legend.position = "bottom")
  
  p <- p + scale_fill_manual(values = cols) + 
      scale_colour_manual(values = text_cols) +
      scale_y_continuous(breaks = seq(0, 100, 10), labels = paste0(seq(0, 100, 10), "%"))
    
  p <- p + labs(x = "SCENARIO",
                y = "MEDIAN REDUCTION (%)",
                fill = "")
    
  p <- p + guides(fill = guide_legend(nrow = 2)) 
    
  p
    
  ggsave(filename = paste0("./analysisworkflow/analysis_scripts/iTPP3_Figures/iTPP3_", exp, "_FIG731.jpg"),
         plot = last_plot(),
         width = 9.1,
         height = 5,
         dpi = 200)

}