############################################################
# PS_07_PlotOptimisation
#
# Visualises relationships between emulator input and optimised variables
# Note: This script depends on outputs of the script 5_optimisation_workflow.R
#
# Written by Lydia Braunack-Mayer
############################################################

# ----------------------------------------------------------
# Setup
# ----------------------------------------------------------

rm(list = ls())

# !!! Insert your experiment name here as a string, e.g. "MyExperiment" !!!
exp <- "iTPP3_tradeoffs_4rounds"

# !!! Insert your predicted parameter here. Note that this must match with one column name in post-processing files !!!
pred <- "inc_red_int_Tot"

library(ggplot2)
library(tidyr)
library(dplyr)
library(scales)

user <- strsplit(getwd(), "/", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]][5]
GROUP_dr <- "/scicore/home/penny/GROUP/M3TPP/"
setwd(paste0("/scicore/home/penny/", user, "/M3TPP/"))

if (!dir.exists(paste0("./Experiments/", exp, "/Outputs"))) dir.create(paste0("./Experiments/", exp, "/Outputs"))

# load(paste0(GROUP_dr, exp, "/param_ranges.RData"))
# param_ranges_cont


# ----------------------------------------------------------
# Define data to plot
# ----------------------------------------------------------

setting <- Sys.glob(paste0(GROUP_dr, exp, "/gp/GP_grid_optimization/", pred, "/*"))
setting_id <- sub(".rds", "", sub("opt_", "", basename(setting)))

for (i in 1:length(setting)) {
  
  df <- readRDS(setting[i])$scenarios
  head(df)
  
  # ----------------------------------------------------------
  # Heat map
  # ----------------------------------------------------------
  
  # Filter data based on coverage
  
  Coverage1 <- Coverage2 <- c("low" = 0.75, "moderate" = 0.85, "high" = 0.95)
  
  df <- df[df$Coverage1 %in% Coverage1 & df$Coverage2 %in% Coverage2, ]
  df$target_range <- floor(df$mean/5)*5
  
  # Format data in preparation for plotting
  df$Coverage1 <- factor(df$Coverage1, labels = c("0.75" = "75% COHORT COVERAGE",
                                                  "0.85" = "85% COHORT COVERAGE",
                                                  "0.95" = "95% COHORT COVERAGE"))
  df$Coverage2 <- factor(df$Coverage2, labels = c("0.75" = "75% ROUND COVERAGE",
                                                  "0.85" = "85% ROUND COVERAGE",
                                                  "0.95" = "95% ROUND COVERAGE"))
  df$Coverage2 <- factor(df$Coverage2, levels = rev(levels(df$Coverage2)))
  
  # Plot
  p <- ggplot(df, aes(x = Halflife, y = Efficacy, fill = target_range))
  
  p <- p + geom_tile(colour = "white")
  
  p <- p + facet_grid(Coverage2 ~ Coverage1)
  
  p <- p + theme(panel.border = element_blank(), 
                 panel.background = element_blank(),
                 panel.grid = element_blank(),
                 text = element_text(colour = "#323f4f", family = "Courier", size = 16),
                 axis.line = element_blank(),
                 axis.ticks = element_blank(),
                 axis.title = element_text(face="bold"),
                 strip.background = element_blank(),
                 legend.key = element_blank(),
                 legend.title = element_text(face = "bold"),
                 legend.key.height = unit(1.5,"cm"))
  
  p <- p + scale_y_continuous(labels = scales::number_format(accuracy = 1, scale = 100, suffix = "%")) +
    scale_fill_binned(high = "#1c649e", 
                      low = "#ffffff", 
                      breaks = seq(min(df$target_range), max(df$target_range), 5)) +
    scale_x_continuous(breaks = seq(10, 60, 10))
  
  p <- p + labs(x = "DURATION OF PROTECTION HALFLIFE (DAYS)",
                y = "INITIAL EFFICACY")
  
 p <- p + guides(fill = guide_colorbar(title = "REDUCTION",
                                       frame.colour = "white"))
  
  ggsave(filename= paste0("./Experiments/", exp, "/Outputs/PS_06_optimisation_", setting_id[i], ".jpg"),
         plot = last_plot(),
         width = 15,
         height = 9,
         dpi = 800)
  
}


