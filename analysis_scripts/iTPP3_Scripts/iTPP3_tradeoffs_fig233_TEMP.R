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
exp_list <- c("iTPP3_tradeoffs_4rounds")


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
  df_plot <- df_plot[df_plot$Outcome %in% c("prev", "inc", "sev"), ]
  
  
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
  
  df_plot$Outcome <- factor(df_plot$Outcome, levels = c("prev", "inc", "sev"))
  df_plot$Outcome <- recode(df_plot$Outcome,
                            "prev" = "PREVALENCE",
                            "inc" = "CLINICAL INCIDENCE",
                            "sev" = "SEVERE DISEASE")
    
    
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
    
  p <- ggplot(df_plot, 
              aes(x = annual_prev, y = T_eff_scaled, fill = parameter, label = paste0(round(T_eff*100, 0), "%")))
    
  p <- p + geom_bar(position = "stack", stat = "identity", colour = "white")
  
  p <- p + geom_text(aes(colour = parameter), 
                     position = position_stack(vjust = 0.5),
                     family = "Courier",
                     fontface = "bold",
                     show.legend = FALSE)
    
  p <- p + facet_grid(Outcome ~ Seasonality, scales = "free_x")
    
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
    
  p
    
  ggsave(filename = paste0("./analysisworkflow/analysis_scripts/iTPP3_Figures/iTPP3_", plot_name, "_TEMP.jpg"),
         plot = last_plot(),
         width = 9.1,
         height = 7.5,
         dpi = 200)

}

# 
# # ----------------------------------------------------------
# # GENERATE RESULTS FOR OPTIMISATION
# # ----------------------------------------------------------
# 
# 
# for (exp in exp_list) {
#   
#   pred <- "prev_red_int_Tot"
#   
#   setting <- Sys.glob(paste0(GROUP_dr, exp, "/gp/GP_grid_optimization/", pred, "/*"))
#   (setting_id <- sub(".rds", "", sub("opt_", "", basename(setting))))
#   index <- seq(4, 56, by = 4)
#   
#   df <- data.frame()
#   
#   for (i in index) {
#     print(paste0("Preparing dataset ", i/4, " of ", length(index)))
#     temp <- readRDS(setting[i])$scenarios
#     temp$scenario <- sub(paste0("_", pred), "", setting_id[i])
#     
#     temp <- temp %>%
#       separate(col = scenario, 
#                into = c("Experiment", "Seasonality", "System", "EIR", "Agegroup", "Decay", "Access", "Timing"),
#                sep = "_",
#                remove = FALSE)
#     
#     df <- rbind(df, temp)
#     remove(temp)
#   }
#   
#   
#   # ----------------------------------------------------------
#   # Prepare data
#   # ----------------------------------------------------------
#   
#   # Update target_range
#   df$target_range <- floor(df$mean/5)*5
#   
#   # Format data in preparation for plotting
#   df$Seasonality <- ifelse(df$Seasonality == "sharpseasonal", "SHORT SEASON", "LONG SEASON")
#   df$Access <- ifelse(df$Access == 0.04, "LOW ACCESS", "HIGH ACCESS")
#   df$Agegroup <- paste0("CHILDREN 3M TO ", df$Agegroup, "Y")
#   df$Seasonality <- factor(df$Seasonality, levels = c("SHORT SEASON", "LONG SEASON"))
#   df$Access <- factor(df$Access, levels = c("LOW ACCESS", "HIGH ACCESS"))
#   
#   #saveRDS(df, paste0("./Experiments/", exp, "/Outputs/PS_06_optimisation_", exp, "_", pred, "_DATA_clinicaltranslation.rds"))
#   
#   df <- df[df$Coverage1 == 1 & df$Coverage2 == 0.95, ]
#   
#   
#   # ----------------------------------------------------------
#   # Import baseline prevalence
#   # ----------------------------------------------------------
#   
#   # Import from csv
#   prev <- read.csv(paste0("./Experiments/", exp, "/Outputs/Prevalence_prior_to_intervention.csv"))
#   prev <- prev[, c("Seasonality", "EIR", "Access", "annual_prev")]
#   prev$annual_prev <- paste0(round(prev$annual_prev*100, 0), "%")
#   prev$Access <- ifelse(prev$Access == 0.04, "LOW ACCESS", "HIGH ACCESS")
#   prev$Seasonality <- ifelse(prev$Seasonality == "sharpseasonal", "SHORT SEASON", "LONG SEASON")
#   prev <- prev[prev$Access == "HIGH ACCESS", ]
#   
#   # Merge into data
#   df <- merge(df, prev, by = c("Seasonality", "EIR", "Access"))
#   df$annual_prev <- paste0(df$annual_prev, " BASELINE\nPfPR 2-10")
#   
#   # ----------------------------------------------------------
#   # Generate heat map
#   # ----------------------------------------------------------
# 
#   # SHORT SEASON
#   
#   p <- ggplot(df[df$Seasonality == "SHORT SEASON", ], aes(x = Halflife, y = Efficacy, fill = target_range))
#   
#   p <- p + geom_tile()
# 
#   p <- p + facet_wrap(.~ annual_prev, ncol = 7)
#   
#   p <- p + theme(panel.border = element_blank(), 
#                  panel.background = element_blank(),
#                  panel.grid = element_blank(),
#                  text = element_text(family = "Courier", size = 12),
#                  strip.background = element_blank(),
#                  strip.text = element_text(face = "bold"),
#                  axis.line = element_blank(),
#                  axis.ticks = element_blank(),
#                  axis.text.x = element_text(margin = margin(t = 0)),
#                  axis.text.y = element_text(margin = margin(r = 0)),
#                  axis.title.x = element_blank(),
#                  axis.title.y = element_text(face = "bold", margin = margin(r = 10)),
#                  plot.title = element_text(hjust = 0.5, face = "bold"),
#                  legend.title = element_text(face = "bold"))
#   
#   p <- p + scale_y_continuous(labels = scales::number_format(accuracy = 1, scale = 100, suffix = "%")) +
#     scale_fill_stepsn(colours = c("#131f41", "#5880b1", "#ffffff"), #"#1c649e"),
#                       breaks = seq(min(df$target_range), max(df$target_range), 5),
#                       labels = function(x) paste0(x, "%"),
#                       limits = c(min(df$target_range), max(df$target_range)),
#                       show.limits = TRUE) +
#     scale_x_continuous(breaks = seq(10, 60, 20))
#   
#   p <- p + labs(y = "INITIAL EFFICACY (%)",
#                 title = "SHORT SEASON")
#   
#   p <- p + guides(fill = "none")
#   
#   
#   # LONG SEASON
#   
#   q <- ggplot(df[df$Seasonality == "LONG SEASON", ], aes(x = Halflife, y = Efficacy, fill = target_range))
#   
#   q <- q + geom_tile()
#   
#   q <- q + facet_wrap(.~ annual_prev, ncol = 7)
#   
#   q <- q + theme(panel.border = element_blank(), 
#                  panel.background = element_blank(),
#                  panel.grid = element_blank(),
#                  text = element_text(family = "Courier", size = 12),
#                  strip.background = element_blank(),
#                  strip.text = element_text(face = "bold"),
#                  axis.line = element_blank(),
#                  axis.ticks = element_blank(),
#                  axis.text.x = element_text(margin = margin(t = 0)),
#                  axis.text.y = element_text(margin = margin(r = 0)),
#                  axis.title.x = element_text(face = "bold", margin = margin(t = 10)),
#                  axis.title.y = element_text(face = "bold", margin = margin(r = 10)),
#                  plot.title = element_text(hjust = 0.5, face = "bold"),
#                  legend.title = element_text(face = "bold"),
#                  legend.position = "bottom",
#                  legend.key.width = unit(1, "cm"))
#   
#   q <- q + scale_y_continuous(labels = scales::number_format(accuracy = 1, scale = 100, suffix = "%")) +
#     scale_fill_stepsn(colours = c("#131f41", "#5880b1", "#ffffff"), #"#1c649e"),
#                       breaks = seq(min(df$target_range), max(df$target_range), 5),
#                       labels = function(x) paste0(x, "%"),
#                       limits = c(min(df$target_range), max(df$target_range)),
#                       show.limits = TRUE) +
#     scale_x_continuous(breaks = seq(10, 60, 20))
#   
#   q <- q + labs(x = "DURATION OF PROTECTION (DAYS)",
#                 y = "INITIAL EFFICACY (%)",
#                 title = "LONG SEASON")
#   
#   q <- q + guides(fill = guide_colorbar(title = "REDUCTION",
#                                         frame.colour = "white"))
#   
#   
#   p + q + plot_layout(ncol = 1.5)
#   
#   ggsave(filename= paste0("./Experiments/", exp, "/Outputs/PS_06_optimisation_", exp, "_", pred, "_FIGOPT_clinicaltranslation.jpg"),
#          plot = last_plot(),
#          width = 9,
#          height = 6,
#          dpi = 300)
# }
