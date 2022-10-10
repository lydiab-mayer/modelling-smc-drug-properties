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
exp <- "iTPP3_ChemoBlood_4rounds"

# !!! Insert your predicted parameter here. Note that this must match with one column name in post-processing files !!!
pred <- "inc_red_int_Tot"

library(hetGP)
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
index <- c(28, 36, 44, 24, 32, 40) # c(27, 35, 43, 23, 31, 39)
setting_id[index]

# Define function to generate predictions
predict.grid <- function(param.ranges, grid.ranges, ngrid, model, scale = TRUE) {
  
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
# cols <- c("0%" = "#010203", "5%" = "#0f1f2f", "10%" = "#152c43", "15%" = "#1d3d5d", "20%" = "#295784",
#           "25%" = "#3e6897", "30%" = "#537bac", "35%" = "#938998", "40%" = "#a6899a", "45%" = "#b98a93",
#           "50%" = "#c78d85", "55%" = "#ca9574", "60%" = "#f9a24b", "65%" = "#fab464", "70%" = "#fac67c",
#           "75%" = "#fad694", "80%" = "#fcdfaa", "85%" = "#fee8c0", "90%" = "#fff4e1")
cols <- c("0%" = "#010203", "10%" = "#0f1f2f", "20%" = "#1d3d5d", "30%" = "#295784", "40%" = "#537bac",
          "50%" = "#938998", "60%" = "#ca9574", "70%" = "#f9a24b", "80%" = "#fad694", "90%" = "#fee8c0",
          "100%" = "#fff4e1")

# Define legend titles
leg_title <- c("inc_red_int_Tot" = "CLINICAL INCIDENCE REDUCTION",
               "sev_red_int_Tot" = "SEVERE\nDISEASE\nREDUCTION",
               "prev_red_int_Aug" = "PREVALENCE\nREDUCTION",
               "mor_red_int_Tot" = "MORTALITY\nREDUCTION")




# ----------------------------------------------------------
# Generate figure 3
# ----------------------------------------------------------

# Generate grid of predictions
ngrid <- c(3, 3, 20, 29, 1)
grid_ranges_cont <- rbind(Coverage1 = c(0.75, 0.95),
                          Coverage2 = c(0.75, 0.95),
                          Halflife = c(1, 20),
                          MaxKillingRate = c(2, 30),
                          Slope = c(6, 6))

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

# Format resulting database
df$Coverage1 <- factor(paste0(df$Coverage1*100, "% PROGRAM REACH"), 
                       levels = c("75% PROGRAM REACH", "85% PROGRAM REACH", "95% PROGRAM REACH"))
df$Coverage2 <- factor(paste0(df$Coverage2*100, "% ROUND COVERAGE"),
                       levels = c("95% ROUND COVERAGE", "85% ROUND COVERAGE", "75% ROUND COVERAGE"))

# For poster, show only 95% round coverage
df <- df[df$Coverage2 == "95% ROUND COVERAGE", ]

# Set target reduction bands
df$target <- round(df$mean / 10, 0)*10 #floor(round(df$mean, 0) / 10) * 10
df$target[df$target < 0] <- 0
df$target_label <- factor(paste0(df$target, "%"), levels = rev(paste0(unique(df$target)[order(unique(df$target))], "%")))

# Generate plot
df_plot <- df[df$EIR == 8, ]
df_cols <- cols[names(cols) %in% unique(df_plot$target_label)]

p <- ggplot(df_plot, aes(x = Halflife, y = MaxKillingRate, fill = target_label))

p <- p + geom_tile()

p <- p + facet_grid(. ~ Coverage1)

p <- p + geom_rect(aes(xmin = 4.74, xmax = 8.81, ymin = 2.28, ymax = 29.96), colour = "white", fill = NA, size = 0.7) + #uncomment for blood stage only
#p <- p + geom_rect(aes(xmin = 5.12, xmax = 8.81, ymin = 2.28, ymax = 29.96), colour = "white", fill = NA, size = 0.7) + #uncomment for dominant blood stage
  annotate(geom = "label", x = 6.9, y = 15, fill = "white", label = "SP-AQ", size = 5, label.size = 0, family = "Arial")

p <- p + theme(panel.border = element_blank(), 
               plot.background = element_rect(fill = "#f1f2f2", colour = "#f1f2f2"),
               panel.background = element_rect(fill = "#f1f2f2"),
               panel.grid = element_blank(),
               text = element_text(family = "Arial", size = 24),
               strip.background = element_blank(),
               axis.line = element_blank(),
               axis.ticks = element_blank(),
               axis.text.x = element_text(margin = margin(t = 0)),
               axis.text.y = element_text(margin = margin(r = 0)),
               axis.title.x = element_text(margin = margin(t = 10)),
               axis.title.y = element_text(margin = margin(r = 10)),
               plot.title = element_text(hjust = 0.5, face = "bold"),
               legend.background = element_rect(fill = "#f1f2f2"),
         #      legend.title = element_text(face = "bold"),
               legend.position = "bottom")

p <- p + scale_fill_manual(values = rev(df_cols)) +
  scale_x_continuous(breaks = seq(2, 20, 4)) + 
  scale_y_continuous(breaks = seq(0, 30, 10))

p <- p + labs(x = "ELIMINATION  HALF-LIFE  (DAYS)",
              y = expression(E["max"]))

p <- p + guides(fill = guide_legend(title = leg_title[pred], nrow = 1, reverse = TRUE))

p

ggsave(filename = paste0("./analysisworkflow/analysis_scripts/iTPP3_Publication/Figures/fig3_POSTER", exp, ".jpg"),
       plot = last_plot(),
       width = 17.7,
       height = 4.5,
       dpi = 400)
