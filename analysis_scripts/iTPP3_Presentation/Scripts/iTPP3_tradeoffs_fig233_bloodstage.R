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
exp <- c("iTPP3_bloodstage")


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
  setting <- Sys.glob(paste0(GROUP_dr, exp, "/gp/trained/sensitivity_CLINICALTRANSLATION/*"))
  setting <- setting[setting != paste0(GROUP_dr, exp, "/gp/trained/sensitivity_CLINICALTRANSLATION/err")]
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
             into = c("Experiment", "Seasonality", "sp", "System", "EIR", "Agegroup", "Access", "Timing", "IC50", "Outcome", "temp1", "temp2", "temp3"),
             sep = "_",
             remove = FALSE)
  df <- df[, !(names(df) %in% c("temp1", "temp2", "temp3"))]
    
    
  # Import aggregated impact for each setting
  
  setting_id <- unique(sub("_0.020831339.*", "_0.020831339", setting_id))
  df_impact <- data.frame()
    
  for (i in 1:length(setting_id)) {
      
    temp <- read.table(paste0(GROUP_dr, exp, "/postprocessing/agg_",setting_id[i], ".txt"), header = TRUE, sep = "")
    temp <- temp %>%
      group_by(Seasonality, Biting_pattern, EIR, MaxAge, Access, Timing) %>%
      summarise(inc = median(inc_red_int_Tot),
                sev = median(sev_red_int_Tot),
                mor = median(mor_red_int_Tot),
                prev = median(prev_red_int_Tot))
    temp$scenario <- setting_id[i]
      
    df_impact <- rbind(df_impact, as.data.frame(temp))
      
  }
    
  df_impact <- df_impact[, c("inc", "sev", "mor", "prev", "scenario")] %>%
    separate(col = scenario, 
             into = c("Experiment", "Seasonality", "sp", "System", "EIR", "Agegroup", "Access", "Timing", "IC50"),
             sep = "_",
             remove = FALSE)
    
  df_impact <- df_impact %>%
    pivot_longer(cols = c(inc, sev, mor, prev), names_to = "Outcome", values_to = "Median_Reduction")
    
  
  # Scale total effects by median impact for each setting
  
  df <- merge(df[, !(names(df) == "scenario")],
              df_impact, 
              by = c("Experiment", "Seasonality", "sp", "System", "EIR", "Agegroup", "Access", "Timing", "IC50", "Outcome"))
    
  df$T_eff_scaled <- df$T_eff * df$Median_Reduction
    
    
  # ----------------------------------------------------------
  # Import baseline prevalence
  # ----------------------------------------------------------
  #   
  # # Import from csv
  # prev <- read.csv(paste0("./Experiments/", exp, "/Outputs/Prevalence_prior_to_intervention.csv"),
  #                  row.names = 1)
  # names(prev)[names(prev) == "annual_prev_210_2034"] <- "annual_prev"
  # 
  # prev <- prev[, c("Seasonality", "EIR", "Access", "annual_prev")]
  # prev$annual_prev <- paste0(round(prev$annual_prev*100, 0), "%")
  # 
  # # Merge into data
  # df <- merge(df, prev, by = c("Seasonality", "EIR", "Access"))
    
    
  # ----------------------------------------------------------
  # Retain only data relevant for clinical translation
  # ----------------------------------------------------------

  df <- df[df$parameter %in% c("Halflife", "MaxKillingRate", "Slope"), ]
  df <- df[df$Access == 0.24, ]
  df <- df[df$Agegroup == 5, ]
  df <- df[df$Outcome %in% c("prev", "inc", "sev"), ]
 # df <- df[!(df$annual_prev %in% c( "9%", "12%")), ]
  
  
  # ----------------------------------------------------------
  # Format data for plotting - MUST BE ADJUSTED FOR EACH EXPERIMENT
  # ----------------------------------------------------------
    
  df$parameter <- as.factor(df$parameter); df$scenario <- as.factor(df$scenario)
  # df <- df %>%
  #   separate(col = scenario, 
  #            into = c("Experiment", "Seasonality", "System", "EIR", "Agegroup", "Decay", "Access", "Timing"),
  #            sep = "_",
  #            remove = FALSE)
  
  df$Seasonality <- factor(df$Seasonality, levels = c("seas3mo", "seas5mo"))
  df$Seasonality <- recode(df$Seasonality,
                           "seas3mo" = "THREE-MONTH SEASON",
                           "seas5mo" = "FIVE-MONTH SEASON")
  
  df$Access <- ifelse(df$Access == 0.04, "LOW ACCESS", "HIGH ACCESS")
  df$Access <- factor(df$Access, levels = c("LOW ACCESS", "HIGH ACCESS"))
  
  df$Agegroup <- factor(df$Agegroup, levels = c("5", "10"))
  df$Agegroup <- recode(df$Agegroup,
                        "5" = "CHILDREN 3 TO 59 MONTHS",
                        "10" = "CHILDREN 3 TO 119 MONTHS")
  
  df$parameter <- factor(df$parameter, levels = c("Halflife", "MaxKillingRate", "Slope"))
  df$parameter <- recode(df$parameter,
                         "Halflife" = "Concentration half-life [5 - 40 days]",
                         "MaxKillingRate" = "Max killing rate [2 - 30 units]",
                         "Slope" = "Slope [1 - 8]")

  #df$annual_prev <- factor(df$annual_prev)
  df$EIR <- factor(df$EIR, levels = c(2, 4, 8, 16, 32, 64, 128, 256))
  
  df$Outcome <- factor(df$Outcome, levels = c("prev", "inc", "sev"))
  df$Outcome <- recode(df$Outcome, 
                            "prev" = "PREVALENCE",
                            "inc" = "CLINICAL INCIDENCE", 
                            "sev" = "SEVERE DISEASE")
    
    
  # ----------------------------------------------------------
  # Define plot settings
  # ----------------------------------------------------------
    
  # Define name of your plot
  plot_name <- "FIG233"
  
  # Define colours
  cols <- c(wes_palette("Royal1", n = 4)[c(1, 4)], "#425055")
  text_cols <- c("#323d42", "#7f4a1b", "white")

  
  # ----------------------------------------------------------
  # Generate plot - incidence and severe disease
  # ----------------------------------------------------------
  
  df_plot <- df[df$Outcome %in% c("PREVALENCE", "CLINICAL INCIDENCE", "SEVERE DISEASE"), ]
    
  p <- ggplot(df_plot, 
              aes(x = EIR, y = T_eff_scaled, fill = parameter, label = paste0(round(T_eff*100, 0), "%")))
    
  p <- p + geom_bar(position = "stack", stat = "identity", colour = "white")
  
  p <- p + geom_text(aes(colour = parameter), 
                     position = position_stack(vjust = 0.5),
                     family = "Courier",
                     fontface = "bold",
                     show.legend = FALSE)
    
  p <- p + facet_grid(Seasonality ~ Outcome, scales = "free")
    
  p <- p + theme(panel.border = element_blank(), 
                 panel.background = element_blank(),
                 panel.grid = element_blank(),
                 text = element_text(family = "Arial", size = 16),
                 strip.background = element_blank(),
                 axis.line = element_blank(),
                 axis.ticks = element_blank(),
                 axis.title.x = element_text(face = "bold", margin = margin(t = 10)),
                 axis.title.y = element_text(face = "bold", margin = margin(r = 10)),
                 strip.text = element_text(face = "bold"),
                 legend.key = element_blank(),
                 legend.title = element_text(face = "bold"),
                 legend.text = element_text(size = 10),
                 legend.position = "bottom")
  
  p <- p + scale_fill_manual(values = cols) + 
    scale_colour_manual(values = text_cols) +
    scale_y_continuous(breaks = seq(0, 100, 10), labels = paste0(seq(0, 100, 10), "%"))
  
  p <- p + labs(x = expression(paste(bold("BASELINE ANNUAL "), bolditalic("Pf"), bold("PR"["2-10"]))),
                y = "MEDIAN REDUCTION (%)",
                fill = "")
    
  p
    
  ggsave(filename = paste0("./analysis_workflow/analysis_scripts/iTPP3_Presentation/Figures/iTPP#_FIG233_Bloodstage.jpg"),
         plot = last_plot(),
         width = 9.1,
         height = 5,
         dpi = 200)
  
  
  # 
  # # ----------------------------------------------------------
  # # Generate plot - prevalence
  # # ----------------------------------------------------------
  # 
  # df_plot <- df[df$Outcome %in% c("PREVALENCE"), ]
  # df_plot$text <- ifelse(df_plot$T_eff <= 0.03, NA, 
  #                        paste0(round(df_plot$T_eff*100, 0), "%"))
  # 
  # p <- ggplot(df_plot, 
  #             aes(x = annual_prev, y = T_eff_scaled, fill = parameter, label = text))
  # 
  # p <- p + geom_bar(position = "stack", stat = "identity", colour = "white")
  # 
  # p <- p + geom_text(aes(colour = parameter), 
  #                    position = position_stack(vjust = 0.5),
  #                    family = "Courier",
  #                    fontface = "bold",
  #                    show.legend = FALSE)
  # 
  # p <- p + facet_grid(. ~ Outcome, scales = "free_x")
  # 
  # p <- p + theme(panel.border = element_blank(), 
  #                panel.background = element_blank(),
  #                panel.grid = element_blank(),
  #                text = element_text(family = "Courier", size = 16),
  #                strip.background = element_blank(),
  #                axis.line = element_blank(),
  #                axis.ticks = element_blank(),
  #                axis.title.x = element_text(face = "bold", margin = margin(t = 10)),
  #                axis.title.y = element_text(face = "bold", margin = margin(r = 10)),
  #                strip.text = element_text(face = "bold"),
  #                legend.key = element_blank(),
  #                legend.title = element_text(face = "bold"),
  #                legend.text = element_text(size = 10),
  #                legend.position = "bottom")
  # 
  # p <- p + scale_fill_manual(values = cols) + 
  #   scale_colour_manual(values = text_cols) +
  #   scale_y_continuous(breaks = seq(0, 100, 10), labels = paste0(seq(0, 100, 10), "%"))
  # 
  # p <- p + labs(x = expression(paste(bold("BASELINE ANNUAL "), bolditalic("Pf"), bold("PR"["2-10"]))),
  #               y = "MEDIAN\nREDUCTION (%)",
  #               fill = "")
  # 
  # p <- p + guides(fill = guide_legend(nrow = 3)) 
  # 
  # p
  # 
  # ggsave(filename = paste0("./Experiments/", exp, "/Outputs/", exp, "_sens_clintranslation_prev.jpg"),
  #        plot = last_plot(),
  #        width = 4.5,
  #        height = 4,
  #        dpi = 200)
  # 