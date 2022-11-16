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
exp <- "iTPP3_ChemoBlood_TreatLiver_4rounds"

# !!! Insert your predicted parameter here. Note that this must match with one column name in post-processing files !!!
pred <- "inc_red_int_Tot"

library(hetGP)
library(tidyr)
library(dplyr)
library(scales)

user <- strsplit(getwd(), "/", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]][5]
GROUP_dr <- "/scicore/home/penny/GROUP/M3TPP/"
setwd(paste0("/scicore/home/penny/", user, "/M3TPP/SMC_TPP/"))

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


# ----------------------------------------------------------
# Generate data
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
df$Coverage1 <- factor(paste0(df$Coverage1*100, "% ROUND COVERAGE"), 
                       levels = c("75% ROUND COVERAGE", "85% ROUND COVERAGE", "95% ROUND COVERAGE"))
df$Coverage2 <- factor(paste0(df$Coverage2*100, "% CYCLE COVERAGE"),
                       levels = c("95% CYCLE COVERAGE", "85% CYCLE COVERAGE", "75% CYCLE COVERAGE"))

# Set target reduction bands
df$target <- round(df$mean / 10, 0)*10
df$target[df$target < 0] <- 0
df$target_label <- factor(paste0(df$target, "%"), levels = rev(paste0(unique(df$target)[order(unique(df$target))], "%")))


# ----------------------------------------------------------
# Write data to file
# ----------------------------------------------------------

saveRDS(df, "./data_and_visualisation/Manuscript_Figure3/data_fig3.rds")

