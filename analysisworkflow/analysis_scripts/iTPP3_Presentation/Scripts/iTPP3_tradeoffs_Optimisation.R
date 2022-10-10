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
pred <- "sev_red_int_Tot"

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
index <- c(40, 48, 56, 36, 44)
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
df$Agegroup <- paste0("CHILDREN 3M TO ", df$Agegroup, "Y")
df$Seasonality <- factor(df$Seasonality, levels = c("SHORT SEASON", "LONG SEASON"))
df$Access <- factor(df$Access, levels = c("LOW ACCESS", "HIGH ACCESS"))

saveRDS(df, paste0("./Experiments/", exp, "/Outputs/PS_06_optimisation_", exp, "_", pred, "_DATA333.rds"))

# Filter data based on coverage
Coverage1 <- 0.95
Coverage2 <- 0.85
df <- df[df$Coverage1 == Coverage1 & df$Coverage2 == Coverage2, ]



# ----------------------------------------------------------
# Generate main heat map for figure 3.3.3
# ----------------------------------------------------------

df_plot <- df[df$EIR == 8, ]

p <- ggplot(df_plot, aes(x = Halflife, y = Efficacy, fill = target_range))
  
p <- p + geom_tile(colour = "white")

# Add DHA+PPQ to plot
p <- p + geom_rect(data = df_plot[1, ], 
                   aes(xmin = 27.5, xmax = 56.5, ymin = .935, ymax = 1.005), 
                   fill = "white", alpha = 0.5) +
  geom_label(label = "DHA+PPQ", x = 42, y = .97, fill = "white", label.size = NA, family = "Arial", fontface = "bold")

# Add SP+AQ to plot
p <- p + geom_rect(data = df_plot[1, ],
                   aes(xmin = 20.5, xmax = 31.5, ymin = .825, ymax = 1.005),
                   fill = "#C93312", alpha = 0.5) +
  geom_label(label = "SP+AQ", x = 25.5, y = 0.915, fill = "#C93312", label.size = NA, family = "Arial", fontface = "bold", colour = "white")

p <- p + theme(panel.border = element_blank(), 
               panel.background = element_blank(),
               panel.grid = element_blank(),
               text = element_text(family = "Arial", size = 13),
               strip.background = element_blank(),
               axis.line = element_blank(),
               axis.ticks = element_blank(),
               axis.text.x = element_text(margin = margin(t = 0)),
               axis.text.y = element_text(margin = margin(r = 0)),
               axis.title.x = element_text(face = "bold", margin = margin(t = 10)),
               axis.title.y = element_text(face = "bold", margin = margin(r = 10)),
               plot.title = element_text(hjust = 0.5, face = "bold"),
               # legend.key.width = unit(1, "cm"),
               # legend.key.height = unit(1, "cm"),
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
              title = expression(paste(bold("ANNUAL BASELINE "), bolditalic("Pf"), bold("PR"["2-10"]), bold(" 26%"))))

p <- p + guides(fill = guide_colorbar(title = "REDUCTION",
                                     frame.colour = "white"))


# ----------------------------------------------------------
# Generate panels for figure 3.3.3
# ----------------------------------------------------------

# PANEL 1

df_plot <- df[df$EIR == 2, ]

p1 <- ggplot(df_plot, aes(x = Halflife, y = Efficacy, fill = target_range))

p1 <- p1 + geom_tile()

p1 <- p1 + theme(panel.border = element_blank(), 
               panel.background = element_blank(),
               panel.grid = element_blank(),
               text = element_text(family = "Arial", size = 8),
               strip.background = element_blank(),
               axis.line = element_blank(),
               axis.ticks = element_blank(),
               axis.text.x = element_text(margin = margin(t = 0)),
               axis.text.y = element_text(margin = margin(r = 0)),
               axis.title.x = element_blank(),
               axis.title.y = element_blank(),
               plot.title = element_text(hjust = 0.5, face = "bold"),
               legend.key.width = unit(1, "cm"),
               legend.key.height = unit(1.5, "cm"),
               legend.title = element_text(face = "bold"))

p1 <- p1 + scale_y_continuous(labels = scales::number_format(accuracy = 1, scale = 100, suffix = "%")) +
  scale_fill_stepsn(colours = c("#1b4d79", "#5880b1", "#f9a24b", "#f9d48e", "#ffedcb"),
                    breaks = seq(min(df$target_range), max(df$target_range), 5),
                    labels = function(x) paste0(x, "%"),
                    limits = c(min(df$target_range), max(df$target_range)),
                    show.limits = TRUE) +
  scale_x_continuous(breaks = seq(10, 60, 10))

p1 <- p1 + labs(x = "",
              y = "",
              title = expression(paste(bolditalic("Pf"), bold("PR"["2-10"]), bold(" 13%"))))

p1 <- p1 + guides(fill = "none")


# PANEL 2

df_plot <- df[df$EIR == 4, ]

p2 <- ggplot(df_plot, aes(x = Halflife, y = Efficacy, fill = target_range))

p2 <- p2 + geom_tile()

p2 <- p2 + theme(panel.border = element_blank(), 
                 panel.background = element_blank(),
                 panel.grid = element_blank(),
                 text = element_text(family = "Arial", size = 8),
                 strip.background = element_blank(),
                 axis.line = element_blank(),
                 axis.ticks = element_blank(),
                 axis.text.x = element_text(margin = margin(t = 0)),
                 axis.text.y = element_text(margin = margin(r = 0)),
                 axis.title.x = element_blank(),
                 axis.title.y = element_blank(),
                 plot.title = element_text(hjust = 0.5, face = "bold"),
                 legend.key.width = unit(1, "cm"),
                 legend.key.height = unit(1.5, "cm"),
                 legend.title = element_text(face = "bold"))

p2 <- p2 + scale_y_continuous(labels = scales::number_format(accuracy = 1, scale = 100, suffix = "%")) +
  scale_fill_stepsn(colours = c("#1b4d79", "#5880b1", "#f9a24b", "#f9d48e", "#ffedcb"),
                    breaks = seq(min(df$target_range), max(df$target_range), 5),
                    labels = function(x) paste0(x, "%"),
                    limits = c(min(df$target_range), max(df$target_range)),
                    show.limits = TRUE) +
  scale_x_continuous(breaks = seq(10, 60, 10))

p2 <- p2 + labs(x = "",
                y = "",
                title = expression(paste(bolditalic("Pf"), bold("PR"["2-10"]), bold(" 18%"))))

p2 <- p2 + guides(fill = "none")


# PANEL 3

df_plot <- df[df$EIR == 16, ]

p3 <- ggplot(df_plot, aes(x = Halflife, y = Efficacy, fill = target_range))

p3 <- p3 + geom_tile()

p3 <- p3 + theme(panel.border = element_blank(), 
                 panel.background = element_blank(),
                 panel.grid = element_blank(),
                 text = element_text(family = "Arial", size = 8),
                 strip.background = element_blank(),
                 axis.line = element_blank(),
                 axis.ticks = element_blank(),
                 axis.text.x = element_text(margin = margin(t = 0)),
                 axis.text.y = element_text(margin = margin(r = 0)),
                 axis.title.x = element_blank(),
                 axis.title.y = element_blank(),
                 plot.title = element_text(hjust = 0.5, face = "bold"),
                 legend.key.width = unit(1, "cm"),
                 legend.key.height = unit(1.5, "cm"),
                 legend.title = element_text(face = "bold"))

p3 <- p3 + scale_y_continuous(labels = scales::number_format(accuracy = 1, scale = 100, suffix = "%")) +
  scale_fill_stepsn(colours = c("#1b4d79", "#5880b1", "#f9a24b", "#f9d48e", "#ffedcb"),
                    breaks = seq(min(df$target_range), max(df$target_range), 5),
                    labels = function(x) paste0(x, "%"),
                    limits = c(min(df$target_range), max(df$target_range)),
                    show.limits = TRUE) +
  scale_x_continuous(breaks = seq(10, 60, 10))

p3 <- p3 + labs(title = expression(paste(bolditalic("Pf"), bold("PR"["2-10"]), bold(" 37%"))))

p3 <- p3 + guides(fill = "none")


# PANEL 4

df_plot <- df[df$EIR == 32, ]

p4 <- ggplot(df_plot, aes(x = Halflife, y = Efficacy, fill = target_range))

p4 <- p4 + geom_tile()

p4 <- p4 + theme(panel.border = element_blank(), 
                 panel.background = element_blank(),
                 panel.grid = element_blank(),
                 text = element_text(family = "Arial", size = 8),
                 strip.background = element_blank(),
                 axis.line = element_blank(),
                 axis.ticks = element_blank(),
                 axis.text.x = element_text(margin = margin(t = 0)),
                 axis.text.y = element_text(margin = margin(r = 0)),
                 axis.title.x = element_blank(),
                 axis.title.y = element_blank(),
                 plot.title = element_text(hjust = 0.5, face = "bold"),
                 legend.key.width = unit(1, "cm"),
                 legend.key.height = unit(1.5, "cm"),
                 legend.title = element_text(face = "bold"))

p4 <- p4 + scale_y_continuous(labels = scales::number_format(accuracy = 1, scale = 100, suffix = "%")) +
  scale_fill_stepsn(colours = c("#1b4d79", "#5880b1", "#f9a24b", "#f9d48e", "#ffedcb"),
                    breaks = seq(min(df$target_range), max(df$target_range), 5),
                    labels = function(x) paste0(x, "%"),
                    limits = c(min(df$target_range), max(df$target_range)),
                    show.limits = TRUE) +
  scale_x_continuous(breaks = seq(10, 60, 10))

p4 <- p4 + labs(title = expression(paste(bolditalic("Pf"), bold("PR"["2-10"]), bold(" 49%"))))

p4 <- p4 + guides(fill = "none")



# ----------------------------------------------------------
# Arrange final plot for figure 3.3.3
# ----------------------------------------------------------

p

ggsave(filename = paste0("./analysisworkflow/analysis_scripts/iTPP3_Presentation/Figures/iTPP3_tradeoffs_", exp, "_", pred, "_FIG333_1.jpg"),
       plot = last_plot(),
       width = 6,
       height = 3.25,
       dpi = 300)

p1 | p2 | p3 | p4

ggsave(filename = paste0("./analysisworkflow/analysis_scripts/iTPP3_Presentation/Figures/iTPP3_tradeoffs_", exp, "_", pred, "_FIG333_2.jpg"),
       plot = last_plot(),
       width = 6,
       height = 1.5,
       dpi = 300)

# ----------------------------------------------------------
# Generate figure 5.3.3.2
# ----------------------------------------------------------

# Prepare data
df <- readRDS(paste0("./Experiments/", exp, "/Outputs/PS_06_optimisation_", exp, "_", pred, "_DATA333.rds"))

df <- df[df$Coverage1 %in% c(0.75, 0.85, 0.95) & df$Coverage2 %in%  c(0.75, 0.85, 0.95), ]
df$Coverage1 <- factor(paste0(df$Coverage1*100, "% PROGRAM COVERAGE"), 
                       levels = c("75% PROGRAM COVERAGE", "85% PROGRAM COVERAGE", "95% PROGRAM COVERAGE"))
df$Coverage2 <- factor(paste0(df$Coverage2*100, "% ROUND COVERAGE"),
                       levels = c("95% ROUND COVERAGE", "85% ROUND COVERAGE", "75% ROUND COVERAGE"))

# Plot
df_plot <- df[df$EIR == 8, ]

p <- ggplot(df_plot, aes(x = Halflife, y = Efficacy, fill = target_range))

p <- p + geom_tile()

p <- p + facet_grid(Coverage2 ~ Coverage1)

p <- p + theme(panel.border = element_blank(), 
               panel.background = element_blank(),
               panel.grid = element_blank(),
               text = element_text(family = "Arial", size = 10),
               strip.background = element_blank(),
               strip.text = element_text(face = "bold"),
               axis.line = element_blank(),
               axis.ticks = element_blank(),
               axis.text.x = element_text(margin = margin(t = 0)),
               axis.text.y = element_text(margin = margin(r = 0)),
               axis.title.x = element_text(face = "bold", margin = margin(t = 10)),
               axis.title.y = element_text(face = "bold", margin = margin(r = 10)),
               plot.title = element_text(hjust = 0.5, face = "bold"),
               # legend.key.width = unit(1, "cm"),
               # legend.key.height = unit(1, "cm"),
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
              title = expression(paste(bold("ANNUAL BASELINE "), bolditalic("Pf"), bold("PR"["2-10"]), bold(" 26%"))))

p <- p + guides(fill = guide_colorbar(title = "REDUCTION",
                                      frame.colour = "white"))

p

ggsave(filename = paste0("./analysisworkflow/analysis_scripts/iTPP3_Presentation/Figures/iTPP3_tradeoffs_", exp, "_", pred, "_FIG5332.jpg"),
       plot = last_plot(),
       width = 8,
       height = 5.5,
       dpi = 200)
