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
pred_list <- c("inc_red_int_Tot", "sev_red_int_Tot", "mor_red_int_Tot")

library(ggplot2)
library(dplyr)
library(tidyr)
library(wesanderson)

user <- strsplit(getwd(), "/", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]][5]
GROUP_dr <- "/scicore/home/penny/GROUP/M3TPP/"
setwd(paste0("/scicore/home/penny/", user, "/M3TPP"))

if (!dir.exists(paste0("./Experiments/", exp_list[1], "/Outputs"))) dir.create(paste0("./Experiments/", exp_list[1], "/Outputs"))

load(paste0(GROUP_dr, exp_list[1], "/param_ranges.RData"))
param_ranges_cont


##########################
# Loop over all outcomes #
##########################

y_axis <- c("MEDIAN\nCLINICAL INCIDENCE\nREDUCTION", "MEDIAN\nSEVERE DISEASE\nREDUCTION", "MEDIAN\nMORTALITY\nREDUCTION")
names(y_axis) <- pred_list

for (pred in pred_list) {

  # ----------------------------------------------------------
  # Define data to plot
  # ----------------------------------------------------------
  
  df <- data.frame("S_eff" = c(), "T_eff" = c(), scenario = c())
  df_impact <- data.frame()
  
  for (experiment in exp_list) {
    
    # Import settings
    setting <- Sys.glob(paste0(GROUP_dr, experiment, "/gp/trained/sensitivity/seeds*"))
    setting <- setting[grepl(pred, setting, fixed = TRUE)]
    setting_id <- sub(paste0("_", pred, ".*"), "", sub(".*seeds_", "", setting))
    
    # Import total effect sizes for each setting
    for (i in 1:length(setting)) {
      
      load(setting[i]) #loads list called sobol_idx_list
      
      sobol_idx_list <- as.data.frame(sobol_idx_list)
      
      sobol_idx_list$S_eff <- sobol_idx_list$S_eff/sum(sobol_idx_list$S_eff) # rescale so total = 1
      sobol_idx_list$T_eff <- sobol_idx_list$T_eff/sum(sobol_idx_list$T_eff) # rescale so total = 1
      
      sobol_idx_list$scenario <- setting_id[i]
      sobol_idx_list$parameter <- rownames(param_ranges_cont)
      
      df <- rbind(df, sobol_idx_list)
      
    }
    
    # Import median impact for each setting
    for (i in 1:length(setting_id)) {
      
      temp <- read.table(paste0(GROUP_dr, experiment, "/postprocessing/agg_", setting_id[i], ".txt"), header = TRUE, sep = "")
      temp <- temp %>%
        group_by(Seasonality, Biting_pattern, EIR, MaxAge, Decay, Access, Timing) %>%
        summarise(inc_red_int_Med = median(inc_red_int_Tot),
                  sev_red_int_Med = median(sev_red_int_Tot),
                  mor_red_int_Med = median(mor_red_int_Tot))
      temp$scenario <- setting_id[i]
      
      df_impact <- rbind(df_impact, as.data.frame(temp))
      
    }
  }
  
  
  # Scale total effects by median impact for each setting
  df <- merge(df, df_impact, by = "scenario")
  if (pred == "inc_red_int_Tot") {
    df$T_eff_scaled <- df$T_eff * df$inc_red_int_Med
  } else if (pred == "sev_red_int_Tot") {
    df$T_eff_scaled <- df$T_eff * df$sev_red_int_Med
  } else {
    df$T_eff_scaled <- df$T_eff * df$mor_red_int_Med
  }
  
  
  # ----------------------------------------------------------
  # Import baseline prevalence
  # ----------------------------------------------------------
  
  # Import from csv
  prev <- read.csv(paste0("./Experiments/", exp_list[1], "/Outputs/Prevalence_prior_to_intervention.csv"))
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
  
  df$Setting <- paste0(df$Seasonality, "-", df$Access)
  df$Setting <- factor(df$Setting, levels = c("LONG SEASON-LOW ACCESS", 
                                              "SHORT SEASON-LOW ACCESS", 
                                              "LONG SEASON-HIGH ACCESS",
                                              "SHORT SEASON-HIGH ACCESS"))
  
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
  
  df$Rounds <- factor(df$Experiment, levels = c("iTPP3tradeoffs3rounds", "iTPP3tradeoffs4rounds", "iTPP3tradeoffs5rounds"))
  df$Rounds <- recode(df$Rounds, 
                      "iTPP3tradeoffs3rounds" = "3 SMC ROUNDS",
                      "iTPP3tradeoffs4rounds" = "4 SMC ROUNDS",
                      "iTPP3tradeoffs5rounds" = "5 SMC ROUNDS")
  
  df$annual_prev <- factor(df$annual_prev)
  
  
  # ----------------------------------------------------------
  # Define plot settings
  # ----------------------------------------------------------
  
  # !!! Define name of your plot !!!
  plot_name <- "FIG334"
  
  cols <- wes_palette("Royal1", n = 4)[c(2, 1, 3, 4)]
  text_cols <- c("#5f1909", "#323d42", "#827d55", "#7f4a1b")
  
  
  # ----------------------------------------------------------
  # Generate plot - figure 3.3.4
  # ----------------------------------------------------------
  
  p <- ggplot(df[df$Access == "HIGH ACCESS" & (df$Agegroup == "CHILDREN 3 TO 59 MONTHS" & df$EIR == 8), ], 
              aes(x = Rounds, y = T_eff_scaled, fill = parameter, label = paste0(round(T_eff*100, 0), "%")))
  
  p <- p + geom_bar(position = "stack", stat = "identity", colour = "white")

  p <- p + geom_text(aes(colour = parameter), position = position_stack(vjust = 0.5), family = "Courier", fontface = "bold",
                     show.legend = FALSE, size = 3)
  
  p <- p + facet_wrap(.~ Seasonality, scales = "free")
  
  p <- p + theme(panel.border = element_blank(), 
                 panel.background = element_blank(),
                 panel.grid = element_blank(),
                 text = element_text(family = "Courier", size = 12),
                 strip.background = element_blank(),
                 axis.line = element_blank(),
                 axis.ticks = element_blank(),
                 axis.title.x = element_text(face = "bold", margin = margin(t = 10)),
                 axis.title.y = element_text(face = "bold", margin = margin(r = 10)),
                 legend.key = element_blank(),
                 legend.title = element_text(face = "bold", size = 10),
                 legend.position = "bottom",
                 legend.margin = margin(t=-1, r=0, b=0, l=0, unit="cm"))
  
  p <- p + scale_fill_manual(values = cols) + 
    scale_colour_manual(values = text_cols) +
    scale_y_continuous(breaks = seq(0, 100, 10), labels = paste0(seq(0, 100, 10), "%"))
  
  p <- p + labs(x = "",
                y = y_axis[pred],
                fill = "")
  
  p <- p + guides(fill = guide_legend(nrow = 2)) 
  
  p
  
  ggsave(filename = paste0("./analysisworkflow/analysis_scripts/iTPP3_Figures/iTPP3_", pred, "_", plot_name, ".jpg"),
         plot = last_plot(),
         width = 9.1,
         height = 3,
         dpi = 200)
  
}
