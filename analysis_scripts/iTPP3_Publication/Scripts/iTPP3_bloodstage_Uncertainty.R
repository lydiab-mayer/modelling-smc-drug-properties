############################################################
#
# Visualises emulator fit
# Note: This script depends on outputs of the script 3_GP_train_workflow.R
#
# Written by Lydia Braunack-Mayer
############################################################

# ----------------------------------------------------------
# Setup
# ----------------------------------------------------------

rm(list = ls())

# !!! Insert your experiment name here as a string, e.g. "MyExperiment" !!!
exp <- "iTPP3_bloodstage_4rounds"

# !!! Insert your predicted parameter here. Note that this must match with one column name in post-processing files !!!
pred_list <- c("inc_red_int_Tot")

library(ggplot2)
library(tidyr)
library(dplyr)
library(patchwork)

user <- strsplit(getwd(), "/", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]][5]
GROUP_dr <- "/scicore/home/penny/GROUP/M3TPP/"
setwd(paste0("/scicore/home/penny/", user, "/M3TPP/"))

load(paste0(GROUP_dr, exp, "/param_ranges.RData"))
param_ranges_cont

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

# Define common plotting features
cols <- c("#C93312", "#425055")

# Generate figures for each outcome

#for (pred in pred_list){
pred <- pred_list

  # ----------------------------------------------------------
  # Define and load scenarios
  # ----------------------------------------------------------
  
  setting <- Sys.glob(paste0(GROUP_dr, exp, "/gp/GP_grid_optimization/", pred, "/*"))
  (setting_id <- sub(".rds", "", sub("opt_", "", basename(setting))))
  index <- 50
  setting_id[index]
  
  # ----------------------------------------------------------
  # Generate figure for coverage1
  # ----------------------------------------------------------
  
  # Generate grid of predictions
  ngrid <- c(31, 1, 1, 1, 1)
  grid_ranges_cont <- rbind(Coverage1 = c(0.7, 1,0),
                            Coverage2 = 0.9,
                            Halflife = 15,
                            MaxKillingRate = 3.45,
                            Slope = 6)
  
  
  # Generate predictions
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
    
    temp <- temp[, !(names(temp) %in% c("temp1", "temp2", "temp3"))]
    
    df <- rbind(df, temp)
    remove(temp)
  }
  
  # Plot parameter vs. prediction
  
  p <- ggplot(df, aes(x = Coverage1, y = mean, ymin = cl, ymax = cu)) +
    geom_ribbon(alpha = 0.3) +
    geom_point()
  
  p <- p + theme(panel.border = element_blank(), 
                 panel.background = element_blank(),
                 panel.grid.major = element_line(colour = "grey95"),
                 panel.grid.minor = element_blank(),
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
  
  p <- p + scale_x_continuous(breaks = seq(0.7, 1.0, by= 0.1),
                              limits = c(0.7, 1.0),
                              labels = paste0(seq(0.7, 1.0, by= 0.1)*100, "%"),
                              expand = expansion(mult = .03, add = 0)) + 
    scale_y_continuous(breaks = seq(55, 95, by = 5),
                       limits = c(55, 95),
                       labels = paste0(seq(55, 95, by = 5), "%"),
                       expand = expansion(mult = .03, add = 0)) +
    scale_colour_manual(values = cols) +
    scale_fill_manual(values = cols)
  
  p <- p + labs(x = "PROGRAM REACH", y = "CLINICAL INCIDENCE REDUCTION")
  
  p_coverage1 <- p
  
  
  
  # ----------------------------------------------------------
  # Generate figure for coverage2
  # ----------------------------------------------------------
  
  # Generate grid of predictions
  ngrid <- c(1, 31, 1, 1, 1)
  grid_ranges_cont <- rbind(Coverage1 = 0.9,
                            Coverage2 = c(0.7, 1.0),
                            Halflife = 15,
                            MaxKillingRate = 3.45,
                            Slope = 6)
  
  
  # Generate predictions
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
  
  # Plot parameter vs. prediction
  
  p <- ggplot(df, aes(x = Coverage2, y = mean, ymin = cl, ymax = cu)) +
    geom_ribbon(alpha = 0.3) +
    geom_point()
  
  p <- p + theme(panel.border = element_blank(), 
                 panel.background = element_blank(),
                 panel.grid.major = element_line(colour = "grey95"),
                 panel.grid.minor = element_blank(),
                 text = element_text(family = "Times New Roman", size = 10),
                 strip.background = element_blank(),
                 axis.line = element_blank(),
                 axis.ticks = element_blank(),
                 axis.text.x = element_text(margin = margin(t = 0)),
                 axis.text.y = element_blank(),
                 axis.title.x = element_text(margin = margin(t = 10)),
                 axis.title.y = element_text(margin = margin(r = 10)),
                 plot.title = element_text(hjust = 0.5, face = "bold"),
                 legend.title = element_text(face = "bold"))
  
  p <- p + scale_x_continuous(breaks = seq(0.7, 1.0, by= 0.1),
                              labels = paste0(seq(0.7, 1.0, by= 0.1)*100, "%"),
                              expand = expansion(mult = .03, add = 0)) + 
    scale_y_continuous(breaks = seq(55, 95, by = 5),
                       limits = c(55, 95),
                       expand = expansion(mult = .03, add = 0))
  
  p <- p + labs(x = "ROUND COVERAGE", y = "")
  
  p_coverage2 <- p
  
  
  # ----------------------------------------------------------
  # Generate figure for halflife
  # ----------------------------------------------------------
  
  # Generate grid of predictions
  ngrid <- c(1, 1, 36, 1, 1)
  grid_ranges_cont <- rbind(Coverage1 = 0.9,
                            Coverage2 = 0.9,
                            Halflife = c(5, 40),
                            MaxKillingRate = 3.45,
                            Slope = 6)
  
  
  # Generate predictions
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
  
  # Plot parameter vs. prediction
  
  p <- ggplot(df, aes(x = Halflife, y = mean, ymin = cl, ymax = cu)) +
    geom_ribbon(alpha = 0.3) +
    geom_point()
  
  p <- p + theme(panel.border = element_blank(), 
                 panel.background = element_blank(),
                 panel.grid.major = element_line(colour = "grey95"),
                 panel.grid.minor = element_blank(),
                 text = element_text(family = "Times New Roman", size = 10),
                 strip.background = element_blank(),
                 axis.line = element_blank(),
                 axis.ticks = element_blank(),
                 axis.text.x = element_text(margin = margin(t = 0)),
                 axis.text.y = element_blank(),
                 axis.title.x = element_text(margin = margin(t = 10)),
                 axis.title.y = element_text(margin = margin(r = 10)),
                 plot.title = element_text(hjust = 0.5, face = "bold"),
                 legend.title = element_text(face = "bold"))
  
  p <- p + scale_x_continuous(breaks = seq(0, 40, by = 10),
                              limits = c(0, 40),
                              expand = expansion(mult = .03, add = 0)) + 
    scale_y_continuous(breaks = seq(55, 95, by = 5),
                       limits = c(55, 95),
                       expand = expansion(mult = .03, add = 0))
  
  p <- p + labs(x = "CONCENTRATION\nHALF-LIFE (DAYS)", y = "")
  
  p_halflife <- p
  
  
  # ----------------------------------------------------------
  # Generate figure for max killing rate
  # ----------------------------------------------------------
  
  # Generate grid of predictions
  ngrid <- c(1, 1, 1, 29, 1)
  grid_ranges_cont <- rbind(Coverage1 = 0.9,
                            Coverage2 = 0.9,
                            Halflife = 15,
                            MaxKillingRate = c(2, 30),
                            Slope = 6)
  
  
  # Generate predictions
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
  
  # Plot parameter vs. prediction
  
  p <- ggplot(df, aes(x = MaxKillingRate, y = mean, ymin = cl, ymax = cu)) +
    geom_ribbon(alpha = 0.3) +
    geom_point()
  
  p <- p + theme(panel.border = element_blank(), 
                 panel.background = element_blank(),
                 panel.grid.major = element_line(colour = "grey95"),
                 panel.grid.minor = element_blank(),
                 text = element_text(family = "Times New Roman", size = 10),
                 strip.background = element_blank(),
                 axis.line = element_blank(),
                 axis.ticks = element_blank(),
                 axis.text.x = element_text(margin = margin(t = 0)),
                 axis.text.y = element_blank(),
                 axis.title.x = element_text(margin = margin(t = 10)),
                 axis.title.y = element_text(margin = margin(r = 10)),
                 plot.title = element_text(hjust = 0.5, face = "bold"),
                 legend.title = element_text(face = "bold"))
  
  p <- p + scale_x_continuous(breaks = seq(0, 30, by = 10),
                              limits = c(0, 30),
                              expand = expansion(mult = .03, add = 0)) + 
    scale_y_continuous(breaks = seq(55, 95, by = 5),
                       limits = c(55, 95),
                       expand = expansion(mult = .03, add = 0))
  
  p <- p + labs(x = "MAXIMUM PARASITE\nKILLING RATE", y = "")
  
  p_maxkilling <- p
  
  
  
  # ----------------------------------------------------------
  # Generate figure for slope
  # ----------------------------------------------------------
  
  # Generate grid of predictions
  ngrid <- c(1, 1, 1, 1, 9)
  grid_ranges_cont <- rbind(Coverage1 = 0.9,
                            Coverage2 = 0.9,
                            Halflife = 15,
                            MaxKillingRate = 3.45,
                            Slope = c(1, 8))
  
  
  # Generate predictions
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
  
  # Plot parameter vs. prediction
  
  p <- ggplot(df, aes(x = Slope, y = mean, ymin = cl, ymax = cu)) +
    geom_ribbon(alpha = 0.3) +
    geom_point()
  
  p <- p + theme(panel.border = element_blank(), 
                 panel.background = element_blank(),
                 panel.grid.major = element_line(colour = "grey95"),
                 panel.grid.minor = element_blank(),
                 text = element_text(family = "Times New Roman", size = 10),
                 strip.background = element_blank(),
                 axis.line = element_blank(),
                 axis.ticks = element_blank(),
                 axis.text.x = element_text(margin = margin(t = 0)),
                 axis.text.y = element_blank(),
                 axis.title.x = element_text(margin = margin(t = 10)),
                 axis.title.y = element_text(margin = margin(r = 10)),
                 plot.title = element_text(hjust = 0.5, face = "bold"),
                 legend.title = element_text(face = "bold"))
  
  p <- p + scale_x_continuous(breaks = seq(0, 8, by = 2),
                              limits = c(0, 8),
                              expand = expansion(mult = .03, add = 0)) + 
    scale_y_continuous(breaks = seq(55, 95, by = 5),
                       limits = c(55, 95),
                       expand = expansion(mult = .03, add = 0))
  
  p <- p + labs(x = "PD MODEL SLOPE", y = "")
  
  p_slope <- p
  
  
  # ----------------------------------------------------------
  # Generate final figure
  # ----------------------------------------------------------
  
  p_out <- p_coverage1 | p_coverage2 | p_halflife | p_maxkilling | p_slope
  
  saveRDS(p_out, paste0("./analysisworkflow/analysis_scripts/iTPP3_Publication/Figures/uncertainty_", exp, "_", pred, ".rds"))

#}



ggsave(filename = paste0("./analysisworkflow/analysis_scripts/iTPP3_Publication/Figures/uncertainty_", exp, "_", pred, ".jpg"),
       plot = p_out,
       width = 9,
       height = 3,
       dpi = 500)


