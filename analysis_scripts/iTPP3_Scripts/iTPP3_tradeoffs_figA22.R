############################################################
#
# Visualises results of optimisation procedure
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
library(patchwork)

user <- strsplit(getwd(), "/", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]][5]
GROUP_dr <- "/scicore/home/penny/GROUP/M3TPP/"
setwd(paste0("/scicore/home/penny/", user, "/M3TPP/"))

if (!dir.exists(paste0("./Experiments/", exp, "/Outputs"))) dir.create(paste0("./Experiments/", exp, "/Outputs"))


# ----------------------------------------------------------
# Define data to plot
# ----------------------------------------------------------

setting <- Sys.glob(paste0(GROUP_dr, exp, "/gp/GP_grid_optimization/", pred, "/*"))
(setting_id <- sub(".rds", "", sub("opt_", "", basename(setting))))
index <- c(27, 28, 55, 56)
setting_id[index]

df <- data.frame()

for (i in index) {
  temp <- readRDS(setting[i])$scenarios
  temp$scenario <- sub(paste0("_", pred), "", setting_id[i])
  
  temp <- temp %>%
    separate(col = scenario, 
             into = c("Experiment", "Seasonality", "System", "EIR", "Agegroup", "Decay", "Access", "Timing"),
             sep = "_",
             remove = FALSE)
  
  df <- rbind(df, temp)
  remove(temp)
}


# ----------------------------------------------------------
# Prepare data
# ----------------------------------------------------------

# Update target_range
df$target_range <- floor(df$mean/5)*5

# Format data in preparation for plotting
df$Seasonality <- ifelse(df$Seasonality == "sharpseasonal", "SHORT SEASON", "LONG SEASON")
df$Access <- ifelse(df$Access == 0.04, "LOW ACCESS", "HIGH ACCESS")
df$Agegroup <- ifelse(df$Agegroup == 5, "CHILDREN 3 TO 59 MONTHS", "CHILDREN 3 TO 119 MONTHS")
df$Agegroup <- factor(df$Agegroup, levels = c("CHILDREN 3 TO 59 MONTHS", "CHILDREN 3 TO 119 MONTHS"))
df$Seasonality <- factor(df$Seasonality, levels = c("SHORT SEASON", "LONG SEASON"))
df$Access <- factor(df$Access, levels = c("LOW ACCESS", "HIGH ACCESS"))

# Filter data based on coverage
Coverage1 <- 0.95
Coverage2 <- 0.85
df <- df[df$Coverage1 == Coverage1 & df$Coverage2 == Coverage2, ]



# ----------------------------------------------------------
# Generate figure
# ----------------------------------------------------------

p <- ggplot(df, aes(x = Halflife, y = Efficacy, fill = target_range))

p <- p + geom_tile(colour = "white")

p <- p + facet_grid(Seasonality ~ Access)

p <- p + theme(panel.border = element_blank(), 
               panel.background = element_blank(),
               panel.grid = element_blank(),
               text = element_text(family = "Courier", size = 13),
               strip.background = element_blank(),
               axis.line = element_blank(),
               axis.ticks = element_blank(),
               axis.text.x = element_text(margin = margin(t = 0)),
               axis.text.y = element_text(margin = margin(r = 0)),
               axis.title.x = element_text(face = "bold", margin = margin(t = 10)),
               axis.title.y = element_text(face = "bold", margin = margin(r = 10)),
               plot.title = element_text(hjust = 0.5, face = "bold"),
               legend.title = element_text(face = "bold"))

p <- p + scale_y_continuous(labels = scales::number_format(accuracy = 1, scale = 100, suffix = "%")) +
  scale_fill_stepsn(colours = c("#1b4d79", "#5880b1", "#f9a24b", "#f9d48e", "#ffedcb"),
                    breaks = seq(min(df$target_range), max(df$target_range), 5),
                    labels = function(x) paste0(x, "%"),
                    limits = c(min(df$target_range), max(df$target_range)),
                    show.limits = TRUE) +
  scale_x_continuous(breaks = seq(10, 60, 10))

p <- p + labs(x = "DURATION OF PROTECTION (DAYS)",
              y = "INITIAL EFFICACY (%)",
              title = expression(paste(bolditalic("Pf"), bold("PR"["2-10"]), bold(" 44%"))))

p <- p + guides(fill = guide_colorbar(title = "REDUCTION",
                                      frame.colour = "white"))

p

ggsave(filename = paste0("./analysisworkflow/analysis_scripts/iTPP3_Figures/iTPP3_tradeoffs_", exp, "_", pred, "_FIGA22.jpg"),
       plot = last_plot(),
       width = 9,
       height = 4,
       dpi = 200)
