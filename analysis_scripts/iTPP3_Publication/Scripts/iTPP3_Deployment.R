############################################################
#
# Visualises SMC deployment dynamics
#
# Written by Lydia Braunack-Mayer
############################################################

# ----------------------------------------------------------
# Setup
# ----------------------------------------------------------

# Clear workspace
rm(list = ls())

# !!! Insert your experiment name here as a string, e.g. "MyExperiment" !!!
exp <- "iTPP3_ChemoBlood_4rounds"

# Load required libraries
library(ggplot2)
library(dplyr)
library(patchwork)

# Define user and group
user <- strsplit(getwd(), "/", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]][5]
GROUP_dr <- "/scicore/home/penny/GROUP/M3TPP/"

# Set working directory
setwd(paste0("/scicore/home/penny/",user,"/M3TPP"))

# Source helper functions
source(paste0("./analysisworkflow/analysis_scripts/import_functions.R"))

# Create folder to store outputs
if (!dir.exists(paste0("./Experiments/",exp,"/Outputs"))) dir.create(paste0("./Experiments/",exp,"/Outputs"))

# Import parameter table
param_table <- read.table(paste0(GROUP_dr,exp,"/param_tab.txt"), sep = "\t", header = TRUE, as.is = TRUE, stringsAsFactors = FALSE)


# ----------------------------------------------------------
# Define data to plot
# ----------------------------------------------------------

# !!! Define which setting you want to plot !!! Should be specified as integer (e.g. '1') or integer vector (e.g. 'c(1, 3, 5)')
(setting <- unique(param_table[, !(names(param_table) == "SEED")]))
setting_id <- 7:8
setting[setting_id, ]

# !!! Define which random seeds you want to visualise !!! Should be defined as an integer or integer vector
seeds <- 1:5

# !!! Define which time steps you want to visualise !!! Should be defined as an integer vector
timesteps <- 1:73

# Calculate average input vs. simulated EIR across all seeds for this setting
df <- import_EIRs_cat(exp = exp, 
                      scenario_id = setting[setting_id, "Scenario_Name"], 
                      seeds = seeds, 
                      timesteps = timesteps)

c(paste0("Total input EIR: ", sum(df$Average$input.EIR)), paste0("Total simulated EIR: ", sum(df$Average$simulated.EIR)))

df$Average$scenario_id <- as.factor(df$Average$scenario_id)
df$Average$scenario_id <- recode(df$Average$scenario_id, 
                                 "iTPP3_ChemoBlood_4rounds_7" = "3 MONTH",
                                 "iTPP3_ChemoBlood_4rounds_8" = "5 MONTH")
# ----------------------------------------------------------
# Define plot settings
# ----------------------------------------------------------

month <- 73/12
int <- month*3:7

# ----------------------------------------------------------
# Generate plot
# ----------------------------------------------------------

# Initialise plot
p <- ggplot()

# Add data for average EIR
p <- p + geom_rect(aes(xmin = 3, xmax = 7, ymin = 0, ymax = max(df$Average$simulated.EIR)), fill = "#899DA4", alpha = 0.1) +
geom_line(data = df$Average, aes(x = timestep/(73/12), y = simulated.EIR, colour = scenario_id)) + 
  geom_vline(xintercept = 3:7, colour = "#899DA4", linetype = "dashed")
  
# Add lines marking rounds

# Define plot theme
p <- p + theme(panel.border = element_blank(), 
                 panel.background = element_blank(),
                 panel.grid = element_blank(),
                 text = element_text(family = "Times New Roman", size = 7),
                 strip.background = element_blank(),
                 axis.line = element_blank(),
                 axis.ticks = element_blank(),
               axis.text.y = element_blank(),
                 axis.title.x = element_text(margin = margin(t = 10)),
                 plot.title = element_text(hjust = 0.5, face = "bold"),
                 legend.title = element_text(face = "bold"),
               legend.position = "bottom",
               legend.key = element_blank()) +
  scale_colour_manual(values = c("#DC863B", "#C93312")) +
  scale_x_continuous(breaks = seq(0, 12, by = 2),
                     expand = expansion(mult = .03, add = 0))

# Define plot titles
p <- p + labs(x = "MONTH",
              y = "ENTOMOLOGICAL  INNOCULATION  RATE",
              colour = "SEASONAL  PROFILE")

p

# Save plot
ggsave(filename= paste0("./analysisworkflow/analysis_scripts/iTPP3_Publication/Figures/iTPP3_Deployment.jpg"),
       plot = last_plot(),
       width = 9.1,
       height = 3,
       dpi = 400)
