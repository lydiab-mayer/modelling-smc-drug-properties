############################################################
# PS_00_plotseasonality
#
# Visualises seasonality in input vs. simulated EIR for OpenMalaria simulations from one or more settings
#
# Written by Lydia Burgert
# Adapted by Lydia Braunack-Mayer
############################################################

# ----------------------------------------------------------
# Setup
# ----------------------------------------------------------

# Clear workspace
rm(list = ls())

# !!! Insert your experiment name here as a string, e.g. "MyExperiment" !!!
exp <- "..."

# Load required libraries
library(ggplot2)

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
setting_id <- ...
setting[setting_id, ]

# !!! Define which random seeds you want to visualise !!! Should be defined as an integer or integer vector
seeds <- ...

# !!! Define which time steps you want to visualise !!! Should be defined as an integer vector
timesteps <- ...

# Calculate average input vs. simulated EIR across all seeds for this setting
df <- import_EIRs_cat(exp = exp, 
                      scenario_id = setting[setting_id, "Scenario_Name"], 
                      seeds = seeds, 
                      timesteps = timesteps)

c(paste0("Total input EIR: ", sum(df$Average$input.EIR)), paste0("Total simulated EIR: ", sum(df$Average$simulated.EIR)))


# ----------------------------------------------------------
# Define plot settings
# ----------------------------------------------------------

# Define line types
types <- c("simulated" = "solid", "input" = "dashed")

# !!! Define name of your plot !!!
plot_name <- "..."


# ----------------------------------------------------------
# Generate plot
# ----------------------------------------------------------

# Initialise plot
p <- ggplot()

# Add data for average EIR
p <- p + geom_line(data = df$Average, aes(x = timestep, y = input.EIR, colour = scenario_id, linetype = "input")) + 
  geom_line(data = df$Average, aes(x = timestep, y = simulated.EIR, colour = scenario_id, linetype = "simulated"))

# Add data for simulated EIR for each individual seed
p <- p + geom_line(data = df$All, aes(x = timestep, y = simulated.EIR, colour = scenario_id, linetype = "simulated", group = interaction(scenario_id, seed)), alpha = 0.35, size = 0.25)

# Define plot theme
p <- p + theme(panel.border = element_blank(), 
               panel.background = element_blank(),
               panel.grid = element_blank(),
               axis.line = element_blank(),
               axis.ticks = element_blank(),
               axis.text = element_text(colour = "grey45", margin = margin(t = 5)),
               axis.title = element_text(colour = "grey30", face="bold"), 
               legend.text = element_text(colour = "grey30"),
               legend.title = element_text(colour = "grey30"),
               legend.key = element_blank())

# Define plot titles
p <- p + labs(x = "Time step", 
              y = "EIR", 
              title = plot_name,
              colour = "Setting",
              linetype = "Line type")

# Save plot
ggsave(filename= paste0("./Experiments/", exp, "/Outputs/PS_00_", plot_name, ".jpg"),
       plot = last_plot(),
       width = 9,
       height = 6,
       dpi = 600)
