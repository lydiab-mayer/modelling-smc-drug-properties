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
exp <- "iTPP3_ChemoLiver_TreatLiverBlood_4rounds"

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

load(paste0(GROUP_dr, exp, "/param_ranges.RData"))
param_ranges_cont

# Define and load scenarios
setting <- Sys.glob(paste0(GROUP_dr, exp, "/gp/GP_grid_optimization/", pred, "/*"))
(setting_id <- sub(".rds", "", sub("opt_", "", basename(setting))))
index <- c(1, 3, 5, 6, 2, 4) #c(12, 20, 28, 8, 16, 24)
setting_id[index]

# Define function to generate predictions
predict.grid <- function(param.ranges, grid.ranges, ngrid, model, scale = TRUE) {
  
  require(hetGP)
  
  # Set up
  D <- nrow(param.ranges)
  scale.params <- t(param.ranges)
  scale.grid <- t(grid.ranges)
  
  # Scale grid to c(0, 1)
  if (scale) {
    for (i in 1:D) {
      scale.grid[, i] <- (scale.grid[, i] - scale.params[1, i]) / (scale.params[2, i] - scale.params[1, i])
    }
  }
  
  # Create grid of scenarios
  scenarios <- list()
  
  for (i in 1:D) {
    scenarios[[i]] <- seq(scale.grid[1, i], scale.grid[2, i], length.out = ngrid[i])
  }
  
  scenarios <- expand.grid(scenarios)
  names(scenarios) <- rownames(param.ranges)
  
  
  # Make predictions using emulator
  preds <- predict(x = as.matrix(scenarios), object = model)
  scenarios$mean <- preds$mean
  scenarios$sd2 <- preds$sd2
  scenarios$nugs <- preds$nugs
  
  # Covert parameter values back to original scale
  for (i in rownames(param.ranges)) {
    scenarios[, i] <- scenarios[, i] * (param.ranges[i, 2] - param.ranges[i, 1]) + param.ranges[i, 1]
  }
  
  # Calculate standard error and 95% confidence interval
  scenarios$se <- sqrt(scenarios$sd2 + scenarios$nugs)
  scenarios$cl <- qnorm(0.05, scenarios$mean, scenarios$se)
  scenarios$cu <- qnorm(0.95, scenarios$mean, scenarios$se)
  
  # Calculate target reduction
  scenarios$target <- floor(scenarios$mean / 5) * 5
  
  return(scenarios)
}

# Define plot colours
cols <- c("25%" = "#010203", "30%" = "#0f1f2f", "35%" = "#1d3d5d", "40%" = "#295784", "45%" = "#3e6897", 
          "50%" = "#537bac", "55%" = "#938998", "60%" = "#ca9574", "65%" = "#f9a24b", "70%" = "#fab464",
          "75%" = "#fac67c", "80%" = "#fad694", "85%" = "#fcdfaa", "90%" = "#fee8c0")

# Define legend titles
leg_title <- c("inc_red_int_Tot" = "CLINICAL\nINCIDENCE\nREDUCTION",
               "sev_red_int_Tot" = "SEVERE\nDISEASE\nREDUCTION",
               "prev_red_int_Aug" = "PREVALENCE\nREDUCTION",
               "mor_red_int_Tot" = "MORTALITY\nREDUCTION")


# ----------------------------------------------------------
# Generate figure 3.3.3
# ----------------------------------------------------------

# Generate grid of predictions
ngrid <- c(1, 1, 51, 21)
grid_ranges_cont <- rbind(Coverage1 = c(0.9, 0.9),
              Coverage2 = c(0.9, 0.9),
              Halflife = c(10, 60),
              Efficacy = c(0.8, 1.0))

df <- data.frame()

for (j in setting_id[index]) {
  # Load GP model
  load(paste0(GROUP_dr, exp, "/gp/trained/", pred, "/seeds_", j, "_cv.RData"))
  
  # Generate model predictions
  temp <- predict.grid(param.ranges = param_ranges_cont, 
                     grid.ranges = grid_ranges_cont, 
                     ngrid = ngrid, 
                     model = cv_result$GP_model)
  
  temp$scenario <- j
  
  temp <- temp %>%
    separate(col = scenario,
             into = c("Experiment", "Seasonality", "System", "EIR", "Agegroup", "Decay", "Access", "Timing", "Outcome", "temp1", "temp2", "temp3"),
             sep = "_",
             remove = FALSE)
  
  df <- rbind(df, temp)
  remove(temp)
}
  
df$target_label <- factor(paste0(df$target, "%"), levels = rev(paste0(unique(df$target)[order(unique(df$target))], "%")))

# ----------------------------------------------------------
# Generate main heat map 
# ----------------------------------------------------------

df_plot <- df[df$EIR == 8, ]
df_cols <- cols[names(cols) %in% unique(df$target_label)]


### Plot with mean predictions ###

p <- ggplot(df_plot, aes(x = Halflife, y = Efficacy, fill = target_label))
  
p <- p + geom_tile(colour = "white")

p <- p + theme(panel.border = element_blank(), 
               panel.background = element_blank(),
               panel.grid = element_blank(),
               text = element_text(family = "Times New Roman", size = 10),
               strip.background = element_blank(),
               axis.line = element_blank(),
               axis.ticks = element_blank(),
               axis.text.x = element_text(margin = margin(t = 0)),
               axis.text.y = element_text(margin = margin(r = 0)),
               axis.title.x = element_text(margin = margin(t = 10)),
               axis.title.y = element_text(margin = margin(r = 10)),
               plot.title = element_text(hjust = 0.5, face = "bold"),
               legend.title = element_text(face = "bold"))
  
p <- p + scale_fill_manual(values = rev(df_cols),
                           limits = levels(df$target_label)) +
  # scale_fill_stepsn(colours = c("#1b4d79", "#5880b1", "#f9a24b", "#f9d48e", "#ffedcb"),
  #                          breaks = seq(min(df$target), max(df$target), 5),
  #                          labels = function(x) paste0(x, "%"),
  #                          limits = c(min(df$target), max(df$target)),
  #                          show.limits = TRUE) +
  scale_x_continuous(breaks = seq(5, 60, 5)) + 
  scale_y_continuous(breaks = seq(0.8, 1.0, 0.05),
                     labels = paste0(seq(80, 100, 5), "%"))
  
p <- p + labs(x = "DURATION  OF  PROTECTION  (DAYS)",
              y = "INITIAL  EFFICACY  (%)",
              title = expression(paste(bold("ANNUAL BASELINE "), bolditalic("Pf"), bold("PR"["2-10"]), bold(" 31%"))))

p <- p + guides(fill = guide_legend(title = leg_title[pred]))

p0 <- p


# ----------------------------------------------------------
# Generate panels for figure 3.3.3
# ----------------------------------------------------------

# PANEL 1

df_plot <- df[df$EIR == 1, ]
df_cols <- cols[names(cols) %in% unique(df_plot$target_label)]

p <- ggplot(df_plot, aes(x = Halflife, y = Efficacy, fill = target_label))

p <- p + geom_tile()

p <- p + theme(panel.border = element_blank(), 
               panel.background = element_blank(),
               panel.grid = element_blank(),
               text = element_text(family = "Times New Roman", size = 10),
               strip.background = element_blank(),
               axis.line = element_blank(),
               axis.ticks = element_blank(),
               axis.text.x = element_text(margin = margin(t = 0)),
               axis.text.y = element_text(margin = margin(r = 0)),
               axis.title = element_blank(),
               plot.title = element_text(hjust = 0.5, face = "bold"),
               legend.position = "none")

p <- p + scale_fill_manual(values = rev(df_cols)) +
  # scale_fill_stepsn(colours = c("#1b4d79", "#5880b1", "#f9a24b", "#f9d48e", "#ffedcb"),
  #                          breaks = seq(min(df$target), max(df$target), 5),
  #                          labels = function(x) paste0(x, "%"),
  #                          limits = c(min(df$target), max(df$target)),
  #                          show.limits = TRUE) +
  scale_x_continuous(breaks = seq(5, 20, 5)) + 
  scale_y_continuous(breaks = seq(2, 30, 4))

p <- p + labs(title = expression(paste(bolditalic("Pf"), bold("PR"["2-10"]), bold(" 11%"))))

p1 <- p


# PANEL 2

df_plot <- df[df$EIR == 2, ]
df_cols <- cols[names(cols) %in% unique(df_plot$target_label)]

p <- ggplot(df_plot, aes(x = Halflife, y = Efficacy, fill = target_label))

p <- p + geom_tile()

p <- p + theme(panel.border = element_blank(), 
               panel.background = element_blank(),
               panel.grid = element_blank(),
               text = element_text(family = "Times New Roman", size = 10),
               strip.background = element_blank(),
               axis.line = element_blank(),
               axis.ticks = element_blank(),
               axis.text.x = element_text(margin = margin(t = 0)),
               axis.text.y = element_text(margin = margin(r = 0)),
               axis.title = element_blank(),
               plot.title = element_text(hjust = 0.5, face = "bold"),
               legend.position = "none")

p <- p + scale_fill_manual(values = rev(df_cols)) +
  # scale_fill_stepsn(colours = c("#1b4d79", "#5880b1", "#f9a24b", "#f9d48e", "#ffedcb"),
  #                          breaks = seq(min(df$target), max(df$target), 5),
  #                          labels = function(x) paste0(x, "%"),
  #                          limits = c(min(df$target), max(df$target)),
  #                          show.limits = TRUE) +
  scale_x_continuous(breaks = seq(5, 20, 5)) + 
  scale_y_continuous(breaks = seq(2, 30, 4))

p <- p + labs(title = expression(paste(bolditalic("Pf"), bold("PR"["2-10"]), bold(" 19%"))))

p2 <- p


# PANEL 3

df_plot <- df[df$EIR == 4, ]
df_cols <- cols[names(cols) %in% unique(df_plot$target_label)]

p <- ggplot(df_plot, aes(x = Halflife, y = Efficacy, fill = target_label))

p <- p + geom_tile()

p <- p + theme(panel.border = element_blank(), 
               panel.background = element_blank(),
               panel.grid = element_blank(),
               text = element_text(family = "Times New Roman", size = 10),
               strip.background = element_blank(),
               axis.line = element_blank(),
               axis.ticks = element_blank(),
               axis.text.x = element_text(margin = margin(t = 0)),
               axis.text.y = element_text(margin = margin(r = 0)),
               axis.title = element_blank(),
               plot.title = element_text(hjust = 0.5, face = "bold"),
               # legend.key.width = unit(1, "cm"),
               # legend.key.height = unit(1, "cm"),
               legend.position = "none")

p <- p + scale_fill_manual(values = rev(df_cols)) +
# + scale_fill_stepsn(colours = c("#1b4d79", "#5880b1", "#f9a24b", "#f9d48e", "#ffedcb"),
#                            breaks = seq(min(df$target), max(df$target), 5),
#                            labels = function(x) paste0(x, "%"),
#                            limits = c(min(df$target), max(df$target)),
#                            show.limits = TRUE) +
  scale_x_continuous(breaks = seq(5, 20, 5)) + 
  scale_y_continuous(breaks = seq(2, 30, 4))

p <- p + labs(title = expression(paste(bolditalic("Pf"), bold("PR"["2-10"]), bold(" 44%"))))

p3 <- p


# PANEL 4

df_plot <- df[df$EIR == 16, ]
df_cols <- cols[names(cols) %in% unique(df_plot$target_label)]

p <- ggplot(df_plot, aes(x = Halflife, y = Efficacy, fill = target_label))

p <- p + geom_tile()

p <- p + theme(panel.border = element_blank(), 
               panel.background = element_blank(),
               panel.grid = element_blank(),
               text = element_text(family = "Times New Roman", size = 10),
               strip.background = element_blank(),
               axis.line = element_blank(),
               axis.ticks = element_blank(),
               axis.text.x = element_text(margin = margin(t = 0)),
               axis.text.y = element_text(margin = margin(r = 0)),
               axis.title = element_blank(),
               plot.title = element_text(hjust = 0.5, face = "bold"),
               # legend.key.width = unit(1, "cm"),
               # legend.key.height = unit(1, "cm"),
               legend.position = "none")

p <- p + scale_fill_manual(values = rev(df_cols)) +
# + scale_fill_stepsn(colours = c("#1b4d79", "#5880b1", "#f9a24b", "#f9d48e", "#ffedcb"),
#                            breaks = seq(min(df$target), max(df$target), 5),
#                            labels = function(x) paste0(x, "%"),
#                            limits = c(min(df$target), max(df$target)),
#                            show.limits = TRUE) +
  scale_x_continuous(breaks = seq(5, 20, 5)) + 
  scale_y_continuous(breaks = seq(2, 30, 4))

p <- p + labs(title = expression(paste(bolditalic("Pf"), bold("PR"["2-10"]), bold(" 57%"))))

p4 <- p

# PANEL 5

df_plot <- df[df$EIR == 32, ]
df_cols <- cols[names(cols) %in% unique(df_plot$target_label)]

p <- ggplot(df_plot, aes(x = Halflife, y = Efficacy, fill = target_label))

p <- p + geom_tile()

p <- p + theme(panel.border = element_blank(), 
               panel.background = element_blank(),
               panel.grid = element_blank(),
               text = element_text(family = "Times New Roman", size = 10),
               strip.background = element_blank(),
               axis.line = element_blank(),
               axis.ticks = element_blank(),
               axis.text.x = element_text(margin = margin(t = 0)),
               axis.text.y = element_text(margin = margin(r = 0)),
               axis.title = element_blank(),
               plot.title = element_text(hjust = 0.5, face = "bold"),
               # legend.key.width = unit(1, "cm"),
               # legend.key.height = unit(1, "cm"),
               legend.position = "none")

p <- p + scale_fill_manual(values = rev(df_cols)) +
  # scale_fill_stepsn(colours = c("#1b4d79", "#5880b1", "#f9a24b", "#f9d48e", "#ffedcb"),
  #                          breaks = seq(min(df$target), max(df$target), 5),
  #                          labels = function(x) paste0(x, "%"),
  #                          limits = c(min(df$target), max(df$target)),
  #                          show.limits = TRUE) +
  scale_x_continuous(breaks = seq(5, 20, 5)) + 
  scale_y_continuous(breaks = seq(2, 30, 4))

p <- p + labs(title = expression(paste(bolditalic("Pf"), bold("PR"["2-10"]), bold(" 70%"))))

p5 <- p


# ----------------------------------------------------------
# Arrange final plot for figure 3.3.3
# ----------------------------------------------------------

p0 / (p1 | p2 | p3 | p4 | p5) + plot_layout(heights = unit(c(7, 1), c('cm', 'null')))

ggsave(filename = paste0("./analysisworkflow/analysis_scripts/iTPP3_Publication/Figures/optim_", exp, "_", pred, ".jpg"),
       plot = last_plot(),
       width = 9.1,
       height = 5,
       dpi = 400)


# ----------------------------------------------------------
# Generate coverage figure
# ----------------------------------------------------------

# Generate grid of predictions
ngrid <- c(3, 3, 51, 21)
grid_ranges_cont <- rbind(Coverage1 = c(0.75, 0.95),
                          Coverage2 = c(0.75, 0.95),
                          Halflife = c(10, 60),
                          Efficacy = c(0.8, 1.0))

df <- data.frame()

for (j in setting_id[index]) {
  # Load GP model
  load(paste0(GROUP_dr, exp, "/gp/trained/", pred, "/seeds_", j, "_cv.RData"))
  
  # Generate model predictions
  temp <- predict.grid(param.ranges = param_ranges_cont, 
                       grid.ranges = grid_ranges_cont, 
                       ngrid = ngrid, 
                       model = cv_result$GP_model)
  
  temp$scenario <- j
  
  temp <- temp %>%
    separate(col = scenario,
             into = c("Experiment", "Seasonality", "System", "EIR", "Agegroup", "Access", "Timing", "IC50", "Outcome", "temp1", "temp2", "temp3"),
             sep = "_",
             remove = FALSE)
  
  df <- rbind(df, temp)
  remove(temp)
}

df$Coverage1 <- factor(paste0(df$Coverage1*100, "% PROGRAM REACH"), 
                       levels = c("75% PROGRAM REACH", "85% PROGRAM REACH", "95% PROGRAM REACH"))
df$Coverage2 <- factor(paste0(df$Coverage2*100, "% ROUND COVERAGE"),
                       levels = c("95% ROUND COVERAGE", "85% ROUND COVERAGE", "75% ROUND COVERAGE"))
df$target_label <- factor(paste0(df$target, "%"), levels = rev(paste0(unique(df$target)[order(unique(df$target))], "%")))

# Plot
df_plot <- df[df$EIR == 16, ]
df_cols <- cols[names(cols) %in% unique(df_plot$target_label)]

p <- ggplot(df_plot, aes(x = Halflife, y = Efficacy, fill = target_label))

p <- p + geom_tile()

p <- p + facet_grid(Coverage2 ~ Coverage1)

p <- p + theme(panel.border = element_blank(), 
               panel.background = element_blank(),
               panel.grid = element_blank(),
               text = element_text(family = "Times New Roman", size = 10),
               strip.background = element_blank(),
               axis.line = element_blank(),
               axis.ticks = element_blank(),
               axis.text.x = element_text(margin = margin(t = 0)),
               axis.text.y = element_text(margin = margin(r = 0)),
               axis.title.x = element_text(margin = margin(t = 10)),
               axis.title.y = element_text(margin = margin(r = 10)),
               plot.title = element_text(hjust = 0.5, face = "bold"),
               legend.title = element_text(face = "bold"))

p <- p + scale_fill_manual(values = rev(df_cols)) +
  scale_x_continuous(breaks = seq(5, 60, 5)) + 
  scale_y_continuous(breaks = seq(0.8, 1.0, 0.05),
                     labels = paste0(seq(80, 100, 5), "%"))

p <- p + labs(x = "DURATION  OF  PROTECTION  (DAYS)",
              y = "INITIAL  EFFICACY  (%)")#,
              #title = expression(paste(bold("ANNUAL BASELINE "), bolditalic("Pf"), bold("PR"["2-10"]), bold(" 34%"))))

p <- p + guides(fill = guide_legend(title = leg_title[pred]))

p

ggsave(filename = paste0("./analysisworkflow/analysis_scripts/iTPP3_Publication/Figures/optim_", exp, "_", pred, "_coverage.jpg"),
       plot = last_plot(),
       width = 9.1,
       height = 5,
       dpi = 400)
