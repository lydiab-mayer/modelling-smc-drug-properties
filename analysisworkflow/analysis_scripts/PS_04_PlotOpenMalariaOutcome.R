############################################################
# PS_04_PlotOpenMalariaOutcome
#
# Visualises OpenMalaria monitoring outcomes for a single, user-specified survey measure
#
# Written by Lydia Braunack-Mayer
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

# !!! Define which monitoring outcome you want to visualise !!! Should be specified as an integer. A full list of monitoring options is available at https://github.com/SwissTPH/openmalaria/wiki/MonitoringOptions
measure <- ...

# Calculate average OpenMalaria monitoring outcome across all seeds for this setting, by age group
df <- import_monitoring_outcome(exp = exp,
                                scenario_id = setting[setting_id, "Scenario_Name"],
                                measure = measure,
                                seeds = seeds,
                                timesteps = timesteps)


# ----------------------------------------------------------
# Define plot settings
# ----------------------------------------------------------

# Define line types
types <- c("Average" = "solid", "Seed" = "dashed")

# !!! Define name of your plot !!!
plot_name <- "..."


# ----------------------------------------------------------
# Generate plot
# ----------------------------------------------------------

# Initialise plot
p <- ggplot()

# Add data for average survey outcome, faceted by age group
p <- p + geom_line(data = df$Average, aes(x = timestep, y = outcome, colour = scenario_id, linetype = "Average")) +
  facet_wrap(.~ agegroup, scales = "free_y")

# Add data for survey outcome for each individual seed, faceted by age group
p <- p + geom_line(data = df$All, aes(x = timestep, y = outcome, colour = scenario_id, linetype = "Seed", group = interaction(scenario_id, seed)), alpha = 0.4, size = 0.25)

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
              y = paste0("Survey measure ", measure),
              title = plot_name,
              colour = "Setting",
              linetype = "Line type")

# Save plot
ggsave(filename= paste0("./Experiments/", exp, "/Outputs/PS_04_", plot_name, ".jpg"),
       plot = last_plot(),
       width = 15,
       height = 9,
       dpi = 800)
