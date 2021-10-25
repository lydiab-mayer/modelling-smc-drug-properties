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
exp <- "iTPP3_tradeoffs"

# !!! Insert your predicted parameter here. Note that this must match with one column name in post-processing files !!!
pred <- "inc_red_int_Avg"

library(ggplot2)
library(tidyr)
library(dplyr)

user <- strsplit(getwd(), "/", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]][5]
GROUP_dr <- "/scicore/home/penny/GROUP/M3TPP/"
setwd(paste0("/scicore/home/penny/", user, "/M3TPP/"))

if (!dir.exists(paste0("./Experiments/", exp, "/Outputs"))) dir.create(paste0("./Experiments/", exp, "/Outputs"))

load(paste0(GROUP_dr, exp, "/param_ranges.RData"))
param_ranges_cont


# ----------------------------------------------------------
# Define data to plot
# ----------------------------------------------------------

# Define parameter to optimise
#param_opt <- "Halflife"

setting <- Sys.glob(paste0(GROUP_dr, exp, "/gp/GP_grid_optimization/", pred, "/*"))
# setting <- setting[grepl(pred, setting, fixed = TRUE)]
# setting_id <- sub(paste0("_", pred, ".*"), "", sub(".*seeds_", "", setting))

df <- read.table(setting[5], sep = " ", header = TRUE)
head(df)

# ----------------------------------------------------------
# Heat map
# ----------------------------------------------------------

# Fix all but two parameters

Coverage1 <- 0.8
Coverage2 <- 0.75

df <- df[df$Coverage1 == Coverage1 & df$Coverage2 == Coverage2, ]
#df <- pivot_longer(data = df, cols = c(Coverage1, Coverage2, Efficacy, optimal_Halflife), names_to = "Parameter")
df$target_range <- floor(df$mean/1)*1


# Plot
p <- ggplot(df, aes(x = Halflife, y = Efficacy, fill = target_range))

p <- p + geom_raster()

p <- p + theme(panel.border = element_blank(), 
               panel.background = element_blank(),
               panel.grid = element_blank(),
               text = element_text(colour = "#323f4f", family = "Simplon Norm (Body)"),
               axis.line = element_blank(),
               axis.ticks = element_blank(),
               axis.title = element_text(face="bold"),
               legend.key = element_blank())

p <- p + scale_y_continuous(labels = scales::percent_format(accurary = 1L)) +
 scale_fill_continuous(high = "#bf3127", 
                       low = "#ffffff", 
                       breaks = seq(60, 80, 5), 
                       limits = c(60, 80),
                       labels = paste0(seq(60, 80, 5), "%")) +
  scale_x_continuous(breaks = seq(10, 60, 5))

p <- p + labs(fill = "Reduction",
              x = "Duration of protection halflife",
              y = "Initial efficacy")

p

ggsave(filename= paste0("./Experiments/", exp, "/Outputs/PS_06_optimisation.jpg"),
       plot = last_plot(),
       width = 15,
       height = 9,
       dpi = 800)

