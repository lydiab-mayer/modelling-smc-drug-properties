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
exp_list <- c("iTPP3_tradeoffs_4rounds")

# !!! Insert your predicted parameters here. Note that this must match with one column name in post-processing files !!!
pred_list <- c("inc_red_int_Tot", "prev_red_int_Aug", "sev_red_int_Tot")

library(ggplot2)
library(dplyr)
library(tidyr)
library(wesanderson)

user <- strsplit(getwd(), "/", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]][5]
GROUP_dr <- "/scicore/home/penny/GROUP/M3TPP/"
setwd(paste0("/scicore/home/penny/", user, "/M3TPP"))

exp <- exp_list #for(exp in exp_list) {
  
  if (!dir.exists(paste0("./Experiments/", exp, "/Outputs"))) dir.create(paste0("./Experiments/", exp, "/Outputs"))
  
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
             into = c("Experiment", "Seasonality", "System", "EIR", "Agegroup", "Decay", "Access", "Timing"),
             sep = "_",
             remove = FALSE)
    
  df_impact <- df_impact %>%
    pivot_longer(cols = c(inc, sev, prev, mor), names_to = "Outcome", values_to = "Median_Reduction")
    
  
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
  prev <- prev[, c("Seasonality", "EIR", "Access", "MaxAge", "annual_prev_210_2034")]
  names(prev) <- c("Seasonality", "EIR", "Access", "Agegroup", "annual_prev")
  prev$annual_prev <- paste0(round(prev$annual_prev*100, 0), "%")
  
  # Merge into data
  df <- merge(df, prev, by = c("Seasonality", "EIR", "Access", "Agegroup"))
    
    
  # ----------------------------------------------------------
  # Format data for plotting - MUST BE ADJUSTED FOR EACH EXPERIMENT
  # ----------------------------------------------------------
  
  df$Seasonality <- factor(df$Seasonality, levels = c("sharpseasonal", "wideseasonal"))
  df$Seasonality <- recode(df$Seasonality, "sharpseasonal" = "5 MONTH SEASON", "wideseasonal" = "6 MONTH SEASON")
  
  df$Access <- factor(df$Access, levels = c("0.04", "0.24"))
  df$Access <- recode(df$Access, "0.04" = "LOW ACCESS", "0.24" = "HIGH ACCESS")
  
  df$Setting <- paste0(df$Seasonality, "\n", df$Access)
  df$Setting <- factor(df$Setting, levels = c("5 MONTH SEASON\nLOW ACCESS", 
                                              "5 MONTH SEASON\nHIGH ACCESS",
                                              "6 MONTH SEASON\nLOW ACCESS",
                                              "6 MONTH SEASON\nHIGH ACCESS"))
    
  df$Agegroup <- factor(df$Agegroup, levels = c("5", "10"))
  df$Agegroup <- recode(df$Agegroup,
                        "5" = "CHILDREN 3 TO 59 MONTHS",
                        "10" = "CHILDREN 3 TO 119 MONTHS")
  
  df$parameter <- factor(df$parameter, levels = c("Coverage1", "Halflife", "Coverage2", "Efficacy"))
  df$parameter <- recode(df$parameter,
                         "Coverage1" = "Program reach [70% - 95%]",
                         "Coverage2" = "Round coverage [70% - 95%]",
                         "Halflife" = "Duration of protection [10 - 60 days]",
                         "Efficacy" = "Initial efficacy [80% - 100%]")
  
  df$Outcome <- factor(df$Outcome, levels = c("inc", "prev", "sev", "mor"))
  df$Outcome <- recode(df$Outcome,
                       "inc" = "CLINICAL INCIDENCE",
                       "sev" = "SEVERE DISEASE",
                       "prev" = "PREVALENCE",
                       "mor" = "MORTALITY")
    
  df$annual_prev <- factor(df$annual_prev)
    
    
  # ----------------------------------------------------------
  # Define plot settings
  # ----------------------------------------------------------
    
  # Define colours
  cols <- c(wes_palette("Royal1", n = 4)[c(2, 1, 3, 4)], "#425055")
  text_cols <- c("#5f1909", "#323d42", "#827d55", "#7f4a1b", "white")
  
  # Define y-axis labels
  y_axis <- c("MEDIAN INCIDENCE REDUCTION", "MEDIAN PREVALENCE REDUCTION", "MEDIAN SEVERE DISEASE REDUCTION")
  names(y_axis) <- pred_list
  
  fontsize <- 10
  
  df$label <- ifelse(df$T_eff >= 0.05, paste0(round(df$T_eff*100, 0), "%"), "")
    
  
  # ----------------------------------------------------------
  # Generate plot - figure S2.1
  # ----------------------------------------------------------
    
  df_plot <- df %>%
    filter(Access == "HIGH ACCESS",
           Seasonality == "5 MONTH SEASON",
           Agegroup == "CHILDREN 3 TO 59 MONTHS",
           Outcome %in% c("CLINICAL INCIDENCE", "PREVALENCE", "SEVERE DISEASE"))
  
  p <- ggplot(df_plot, aes(x = annual_prev, y = T_eff_scaled, fill = parameter, label = label))
    
  p <- p + geom_bar(position = "stack", stat = "identity", colour = "white")
  
  p <- p + geom_text(aes(colour = parameter), 
                     position = position_stack(vjust = 0.5),
                     family = "Times New Roman",
                     size = fontsize*0.28,
                     show.legend = FALSE)
    
  p <- p + facet_wrap(.~ Outcome, scales = "free_x")
    
  p <- p + theme(panel.border = element_blank(), 
                 panel.background = element_blank(),
                 panel.grid = element_blank(),
                 text = element_text(family = "Times New Roman", size = fontsize),
                 strip.background = element_blank(),
                 axis.line = element_blank(),
                 axis.ticks = element_blank(),
                 axis.title.x = element_text(margin = margin(t = 10)),
                 axis.title.y = element_text(margin = margin(r = 10)),
                 strip.text = element_text(face = "bold"),
                 legend.key = element_blank(),
                 legend.position = "bottom")
  
  p <- p + scale_fill_manual(values = cols) + 
    scale_colour_manual(values = text_cols) +
    scale_y_continuous(breaks = seq(0, 100, 10), labels = paste0(seq(0, 100, 10), "%"))
  
  p <- p + labs(x = expression(paste("BASELINE  ANNUAL  ", italic("Pf"), "PR"["2-10"])),
                y = "MEDIAN  REDUCTION (%)",
                fill = "")
    
  p <- p + guides(fill = guide_legend(nrow = 2)) 
    
  p
    
  ggsave(filename = paste0("./analysisworkflow/analysis_scripts/iTPP3_Publication/Figures/iTPP3_", exp, "_FIGS21.jpg"),
         plot = last_plot(),
         width = 9.1,
         height = 5,
         dpi = 400)
  

    
    
  # ----------------------------------------------------------
  # Generate table - table A.2.1
  # ----------------------------------------------------------
    
  tab <- df %>%
    group_by(parameter, Outcome) %>%
    summarise(max = max(S_eff), min = min(S_eff))
    
  tab$minmax <- paste0(round(tab$min*100, 0), "% to ", round(tab$max*100, 0), "%")
  names(tab) <- c("Key performance property", "Outcome", "Max", "Min", "Range of attributable outcome variation")
  write.csv(tab[, c(1, 2, 5)], file = paste0("./analysisworkflow/analysis_scripts/iTPP3_Publication/Figures/iTPP3_", exp, "_sensitivity_summary.csv"))
    
    
  # ----------------------------------------------------------
  # Generate plot - figure S2.3
  # ----------------------------------------------------------
    
  df_plot <- df %>%
    filter(Access == "HIGH ACCESS",
           Seasonality == "5 MONTH SEASON",
           Outcome %in% c("CLINICAL INCIDENCE", "PREVALENCE", "SEVERE DISEASE"))
  
  p <- ggplot(df_plot, aes(x = annual_prev, y = T_eff_scaled, fill = parameter, label = label))
    
  p <- p + geom_bar(position = "stack", stat = "identity", colour = "white")
    
  p <- p + geom_text(aes(colour = parameter), position = position_stack(vjust = 0.5), 
                     family = "Times New Roman",
                     size = fontsize*0.28,
                     show.legend = FALSE)
    
  p <- p + facet_grid(Outcome ~ Agegroup, scales = "free_x")
  
  p <- p + theme(panel.border = element_blank(), 
                 panel.background = element_blank(),
                 panel.grid = element_blank(),
                 text = element_text(family = "Times New Roman", size = 10),
                 strip.background = element_blank(),
                 axis.line = element_blank(),
                 axis.ticks = element_blank(),
                 axis.title.x = element_text(margin = margin(t = 10)),
                 axis.title.y = element_text(margin = margin(r = 10)),
                 strip.text = element_text(face = "bold"),
                 legend.key = element_blank(),
                 legend.position = "bottom")
  
  p <- p + scale_fill_manual(values = cols) + 
    scale_colour_manual(values = text_cols) +
    scale_y_continuous(breaks = seq(0, 100, 10), labels = paste0(seq(0, 100, 10), "%"))
    
  p <- p + labs(x = expression(paste("BASELINE ANNUAL ", italic("Pf"), "PR"["2-10"])),
                y = "MEDIAN REDUCTION (%)",
                fill = "")
    
  p <- p + guides(fill = guide_legend(nrow = 2)) 
    
  p
    
  ggsave(filename = paste0("./analysisworkflow/analysis_scripts/iTPP3_Publication/Figures/iTPP3_", exp,"_FIGS23.jpg"),
         plot = last_plot(),
         width = 9.1,
         height = 9.1,
         dpi = 400)
    
    
    
  # ----------------------------------------------------------
  # Generate plot - figure S2.5
  # ----------------------------------------------------------
  
  df_plot <- df %>%
    filter(Agegroup == "CHILDREN 3 TO 59 MONTHS",
           EIR == 16,
           Outcome %in% c("CLINICAL INCIDENCE", "PREVALENCE", "SEVERE DISEASE"))
  
  p <- ggplot(df_plot, aes(x = Setting, y = T_eff_scaled, fill = parameter, label = label))
    
  p <- p + geom_bar(position = "stack", stat = "identity", colour = "white")
  
  p <- p + geom_text(aes(colour = parameter), position = position_stack(vjust = 0.5), 
                     family = "Times New Roman",
                     size = fontsize*0.28,
                     show.legend = FALSE)
    
  p <- p + facet_wrap(.~ Outcome, scales = "free_x", ncol = 1, nrow = 3)
    
  p <- p + theme(panel.border = element_blank(), 
                 panel.background = element_blank(),
                 panel.grid = element_blank(),
                 text = element_text(family = "Times New Roman", size = 10),
                 strip.background = element_blank(),
                 axis.line = element_blank(),
                 axis.ticks = element_blank(),
                 axis.title.x = element_text(margin = margin(t = 10)),
                 axis.title.y = element_text(margin = margin(r = 10)),
                 strip.text = element_text(face = "bold"),
                 legend.key = element_blank(),
                 legend.position = "bottom")
  
  p <- p + scale_fill_manual(values = cols) + 
      scale_colour_manual(values = text_cols) +
      scale_y_continuous(breaks = seq(0, 100, 10), labels = paste0(seq(0, 100, 10), "%"))
    
  p <- p + labs(x = "SCENARIO",
                y = "MEDIAN REDUCTION (%)",
                fill = "")
    
  p <- p + guides(fill = guide_legend(nrow = 2)) 
    
  p
    
  ggsave(filename = paste0("./analysisworkflow/analysis_scripts/iTPP3_Publication/Figures/iTPP3_", exp, "_FIGS25.jpg"),
         plot = last_plot(),
         width = 9.1,
         height = 7,
         dpi = 300)

#}