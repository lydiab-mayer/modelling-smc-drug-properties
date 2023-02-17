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
exp <- "iTPP3_ChemoBlood_TreatLiver_4rounds"

# !!! Insert your predicted parameters here. Note that this must match with one column name in post-processing files !!!
predictors <- c("inc_red_int_Tot", "prev_red_int_Aug", "sev_red_int_Tot")

library(tidyr)
library(dplyr)

user <- strsplit(getwd(), "/", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]][5]
GROUP_dr <- "/scicore/home/penny/GROUP/M3TPP/"
setwd(paste0("/scicore/home/penny/", user, "/M3TPP/SMC_TPP/"))

load(paste0(GROUP_dr, exp, "/param_ranges.RData"))
param_ranges_cont

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

# Define list to store outputs
data <- vector(mode = "list", length = 0L)


# ----------------------------------------------------------
# Define and load scenarios
# ----------------------------------------------------------

setting <- Sys.glob(paste0(GROUP_dr, exp, "/gp/GP_grid_optimization/", predictors[1], "/*"))
(setting_id <- sub(".rds", "", sub("opt_", "", basename(setting))))
index <- 48 #24
setting_id[index]


# ----------------------------------------------------------
# Generate data for coverage1
# ----------------------------------------------------------

# Generate grid of predictions
ngrid <- c(31, 1, 1, 1, 1)
grid_ranges_cont <- rbind(Coverage1 = c(0.7, 1.0),
                          Coverage2 = 0.8,
                          Halflife = 5,
                          MaxKillingRate = 2,
                          Slope = 6)


# Generate predictions
df <- data.frame()

for (pred in predictors) {
  
  setting <- Sys.glob(paste0(GROUP_dr, exp, "/gp/GP_grid_optimization/", pred, "/*"))
  setting_id <- sub(".rds", "", sub("opt_", "", basename(setting)))
  
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
    
    temp$pred <- pred
    
    df <- rbind(df, temp)
    remove(temp)
  }
}


data[["coverage1"]] <- df


# ----------------------------------------------------------
# Generate data for coverage2
# ----------------------------------------------------------

# Generate grid of predictions
ngrid <- c(1, 31, 1, 1, 1)
grid_ranges_cont <- rbind(Coverage1 = 0.8,
                          Coverage2 = c(0.7, 1.0),
                          Halflife = 5,
                          MaxKillingRate = 2,
                          Slope = 6)


# Generate predictions
df <- data.frame()

for (pred in predictors) {
  
  setting <- Sys.glob(paste0(GROUP_dr, exp, "/gp/GP_grid_optimization/", pred, "/*"))
  setting_id <- sub(".rds", "", sub("opt_", "", basename(setting)))
  
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
    
    temp$pred <- pred
    
    df <- rbind(df, temp)
    remove(temp)
  }
}

data[["coverage2"]] <- df


# ----------------------------------------------------------
# Generate data for halflife
# ----------------------------------------------------------

# Generate grid of predictions
ngrid <- c(1, 1, 20, 1, 1)
grid_ranges_cont <- rbind(Coverage1 = 0.8,
                          Coverage2 = 0.8,
                          Halflife = c(1, 20),
                          MaxKillingRate = 2,
                          Slope = 6)


# Generate predictions
df <- data.frame()

for (pred in predictors) {
  
  setting <- Sys.glob(paste0(GROUP_dr, exp, "/gp/GP_grid_optimization/", pred, "/*"))
  setting_id <- sub(".rds", "", sub("opt_", "", basename(setting)))
  
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
    
    temp$pred <- pred
    
    df <- rbind(df, temp)
    remove(temp)
  }
}

data[["halflife"]] <- df


# ----------------------------------------------------------
# Generate figure for max killing rate
# ----------------------------------------------------------

# Generate grid of predictions
ngrid <- c(1, 1, 1, 29, 1)
grid_ranges_cont <- rbind(Coverage1 = 0.8,
                          Coverage2 = 0.8,
                          Halflife = 5,
                          MaxKillingRate = c(2, 30),
                          Slope = 6)


# Generate predictions
df <- data.frame()

for (pred in predictors) {
  
  setting <- Sys.glob(paste0(GROUP_dr, exp, "/gp/GP_grid_optimization/", pred, "/*"))
  setting_id <- sub(".rds", "", sub("opt_", "", basename(setting)))
  
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
    
    temp$pred <- pred
    
    df <- rbind(df, temp)
    remove(temp)
  }
}

data[["MaxKillingRate"]] <- df
  

# ----------------------------------------------------------
# Generate figure for slope
# ----------------------------------------------------------

# Generate grid of predictions
ngrid <- c(1, 1, 1, 1, 9)
grid_ranges_cont <- rbind(Coverage1 = 0.8,
                          Coverage2 = 0.8,
                          Halflife = 5,
                          MaxKillingRate = 2,
                          Slope = c(1, 8))


# Generate predictions
df <- data.frame()

for (pred in predictors) {
  
  setting <- Sys.glob(paste0(GROUP_dr, exp, "/gp/GP_grid_optimization/", pred, "/*"))
  setting_id <- sub(".rds", "", sub("opt_", "", basename(setting)))
  
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
    
    temp$pred <- pred
    
    df <- rbind(df, temp)
    remove(temp)
  }
}

data[["slope"]] <- df


# ----------------------------------------------------------
# Format dataset
# ----------------------------------------------------------

for (param in c("coverage1", "coverage2", "halflife", "MaxKillingRate", "slope")) {
  data[[param]]$pred <- factor(data[[param]]$pred, levels = predictors)
  data[[param]]$pred <- recode_factor(data[[param]]$pred,
                                      "inc_red_int_Tot" = "Clinical incidence",
                                      "prev_red_int_Aug" = "Prevalence",
                                      "sev_red_int_Tot" = "Severe disease")
}


# ----------------------------------------------------------
# Write data to file
# ----------------------------------------------------------

saveRDS(data, "./data_and_visualisation/Appendix_Figure29/data_figA29_panelA.rds")
