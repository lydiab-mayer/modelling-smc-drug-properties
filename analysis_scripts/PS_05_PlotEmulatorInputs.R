############################################################
# PS_05_PlotEmulatorInputs
#
# Visualises relationships between emulator input and predictor variables
# Note: This script depends on outputs of the script 2_OM_postprocessing
#
# Written by Lydia Braunack-Mayer
############################################################

# ----------------------------------------------------------
# Setup
# ----------------------------------------------------------

rm(list = ls())

# !!! Insert your experiment name here as a string, e.g. "MyExperiment" !!!
exp <- "..."

# !!! Insert your predicted parameters here. Note that this must match with one or more column names in post-processing files !!!
pred_list <- c("inc_red_int_Avg", "sev_red_int", "mor_red_int")

library(ggplot2)
library(dplyr)

user <- strsplit(getwd(), "/", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]][5]
GROUP_dr <- "/scicore/home/penny/GROUP/M3TPP/"
setwd(paste0("/scicore/home/penny/", user, "/M3TPP"))

if (!dir.exists(paste0("./Experiments/", exp, "/Outputs"))) dir.create(paste0("./Experiments/", exp, "/Outputs"))

load(paste0(GROUP_dr, exp, "/param_ranges.RData"))
param_ranges_cont


# ----------------------------------------------------------
# Define data to plot
# ----------------------------------------------------------

setting <- sub(".txt.*", "", sub(".*agg_", "", Sys.glob(paste0(GROUP_dr, exp, "/postprocessing/agg*"))))
setting

# !!! Define which setting you want to plot !!! Should be specified as integer (e.g. '1')
setting_id <- 14 #14, 16, 22

df <- read.table(paste0(GROUP_dr, exp, "/postprocessing/agg_", setting[setting_id], ".txt"), sep = "\t", header = TRUE, as.is = TRUE, stringsAsFactors = FALSE)
df <- df[, names(df) %in% c(rownames(param_ranges_cont), pred_list)]
df <- pivot_longer(df, cols = all_of(pred_list), names_to = "outcome", values_to = "outcome_value") %>%
  pivot_longer(cols = all_of(rownames(param_ranges_cont)), names_to = "predictor", values_to = "predictor_value")


# ----------------------------------------------------------
# Define plot settings
# ----------------------------------------------------------

# !!! Define name of your plot !!!
plot_name <- paste0("Emulator inputs vs outcomes for: ", setting[setting_id])


# ----------------------------------------------------------
# Generate plot
# ----------------------------------------------------------

p <- ggplot()

p <- p + geom_point(data = df, aes(x = predictor_value, y = outcome_value), colour = "grey55") +
  facet_grid(outcome ~ predictor, scales = "free")

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

p <- p + labs(x = "Value of predictor",
              y = "Value of outcome",
              title = plot_name,
              colour = "Outcome")
p

ggsave(filename= paste0("./Experiments/", exp, "/Outputs/PS_05_", plot_name, ".jpg"),
       plot = last_plot(),
       width = 15,
       height = 9,
       dpi = 800)

