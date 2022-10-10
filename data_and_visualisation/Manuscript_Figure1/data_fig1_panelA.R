############################################################
#
# Plotting script to demonstrate model profile for SMC with
# dominant blood stage activity
#
# Written by Lydia Braunack-Mayer
# April 2022
#
############################################################

##########
# SET UP #
##########

# Clear workspace
rm(list = ls())

# !!! Insert your experiment name here as a string, e.g. "MyExperiment" !!!
exp <- "iTPP3_ChemoBlood_replication"

# Load required libraries
library(dplyr)
library(tidyr)
library(patchwork)

# Define user and group
user <- strsplit(getwd(), "/", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]][5]
GROUP_dr <- "/scicore/home/penny/GROUP/M3TPP/"

# Set working directory
setwd(paste0("/scicore/home/penny/",user,"/M3TPP/SMC_TPP"))


##############################
# PREPARE VALIDATION RESULTS #
##############################

# ----------------------------------------------------------
# Import data
# ----------------------------------------------------------

setting <- Sys.glob(paste0(GROUP_dr, exp, "/postprocessing/seeds*"))

df <- data.frame()

for (i in setting) {
  temp <- read.table(i, header = TRUE, sep = "")
  df <- rbind(df, temp)
  remove(temp)
}


# ----------------------------------------------------------
# Calculate protective efficacy
# ----------------------------------------------------------

# Index so that outputs for intervention/control groups for a given experiment are consecutive
df <- df[order(df$EIR, df$IC50, df$Halflife, df$MaxKillingRate, df$Slope, df$seed), ]

# Calculate pe between consecutive rows
df <- df %>%
  group_by(EIR, IC50, Halflife, MaxKillingRate, Slope, seed) %>%
  summarise(across(contains("cpp"), ~ 1 - tail(.x, 1)/head(.x, 1), .names = "pe_{col}"))

# Aggregate over seeds
df <- df %>%
  group_by(EIR, IC50, Halflife, MaxKillingRate, Slope) %>%
  summarise(across(contains("cpp"), mean))

# Set pe to zero if < 0
df[df <= 0] <- 0
# for (i in 1:nrow(df)) {
#   index <- which(df[i, ] <= 0)[1]
#   if (!is.na(index)) {
#     df[i, index] <- 0
#     df[i, (index + 1):ncol(df)] <- NA
#   }
# }

# ----------------------------------------------------------
# Import zongo trial data
# ----------------------------------------------------------

# Import protective efficacy profile from Zongo et al. 2015
zongo <- read.csv("./data_and_visualisation/Manuscript_Figure1/zongo_data_extraction.csv", sep = ";", as.is = TRUE)

# Add dummy columns to match columns in simulation database
df$linetype <- "mean"
df$drug <- "Next-gen SMC"
zongo$EIR <- zongo$IC50 <- zongo$Halflife <- zongo$MaxKillingRate <- zongo$Slope <- 0


# ----------------------------------------------------------
# Find line(s) of best fit
# ----------------------------------------------------------

# Extract SP-AQ data
SPAQ <- zongo[zongo$drug == "SP-AQ" & zongo$linetype == "mean", "y"]

# Calculate RSS between SP-AQ line and simulated protective efficacy
RSS <- function(x, y) {
  # Function to calculate Residual-Sum of Squares between two vectors
  out <- sum((x - y)^2, na.rm = TRUE)
  return(out)
}

df$RSS <- NULL

for (i in 1:nrow(df)) {
  df[i, "RSS"] <- RSS(SPAQ, as.numeric(df[i, 6:17]))
}

# Calculate cutoff for 'close' to SP-AQ
cutoff <- min(df$RSS) + 0.1*sd(df$RSS)

# Identify parameter range for simulations 'close' to SP-AQ
df[df$RSS <= cutoff, ]

nextgen <- df %>%
  filter(RSS <= cutoff) %>%
  group_by() %>%
  summarise(Halflife_max = max(Halflife),
            Halflife_min = min(Halflife),
            MaxKillingRate_max = max(MaxKillingRate),
            MaxKillingRate_min = min(MaxKillingRate),
            Slope_max = max(Slope),
            Slope_min = min(Slope))

# Format resulting dataframe
nextgen <- as.data.frame(t(nextgen))
nextgen$parameter <- gsub("_.*", "", rownames(nextgen))
nextgen$type <- gsub(".*_", "", rownames(nextgen))
names(nextgen)[1] <- "value"
nextgen <- nextgen %>%
  pivot_wider(names_from = type) %>%
  as.data.frame()

# Save results
write.csv(nextgen, paste0("./data_and_visualisation/Manuscript_Figure1/table2_modelA.csv"))


# ----------------------------------------------------------
# Prepare data for plotting
# ----------------------------------------------------------

# Prepare database
df <- df %>%
  pivot_longer(cols = starts_with("pe"), names_to = "x", values_to = "y")
df$x <- as.numeric(sub("pe_cpp_", "", df$x))
df$x <- df$x - min(df$x)
df <- as.data.frame(df)

# Merge with zongo trial data
df_plot_pe <- df %>% select(-RSS)
df_plot_pe <- rbind(zongo, df_plot_pe)
df_plot_pe <- pivot_wider(df_plot_pe, id_cols = c(drug, x, Slope, MaxKillingRate, Halflife, IC50, EIR), names_from = linetype, values_from = y)
df_plot_pe$drug <- factor(df_plot_pe$drug, levels = c("SP-AQ", "DHA-PPQ", "Next-gen SMC"))

# Add column for weeks
df$weeks <- round(df$x*5/7, 1)
df_plot_pe$weeks <- round(df_plot_pe$x*5/7, 1)


####################################
# GENERATE RANGE OF PK/PD PROFILES #
####################################

# ----------------------------------------------------------
# Set up
# ----------------------------------------------------------

# Define model PK/PD components
PK <- function(t, C_0, halflife) {
  # t = time
  # C_0 = concentration at time 0
  # halflife = Time t at which half the initial concentration remains
  
  C_0 * exp(-(log(2)/halflife)*t)
}

PD <- function(C, Emax, IC50, n) {
  # C = Initial concentration
  # Emax = Maximum parasite killing rate
  # IC50 = Drug concentration at which half the maximum killing rate is achieved
  # n = slope
  
  Emax * C^n / (C^n + IC50^n)
}

# Define fixed parameter values
dose <- 18 # dose in mg/kg
weight <- 20 # weight of child receiving drug in kg - set for demonstration purposes
vol <- 173 # volume of distribution in l/kg
C_0 <- (dose*weight) / (vol*weight)
IC50 <- 0.020831339


# ----------------------------------------------------------
# Generate complete database of PK/PD profiles
# ----------------------------------------------------------

# Create grid of scenarios

ngrid <- c(4, 8, 8)
grid_ranges_cont <- rbind(Halflife = c(1, 20),
                          MaxKillingRate = c(2, 30),
                          Slope = c(1, 8))

D <- nrow(grid_ranges_cont)
grid_ranges <- t(grid_ranges_cont)
scenarios <- list()

for (i in 1:D) {
  scenarios[[i]] <- seq(grid_ranges[1, i], grid_ranges[2, i], length.out = ngrid[i])
}

scenarios <- expand.grid(scenarios)
names(scenarios) <- rownames(grid_ranges_cont)
head(scenarios)


# Generate PK/PD values for each scenario

df_PKPD <- data.frame()

for(i in 1:nrow(scenarios)) {
  temp <- data.frame("Halflife" = scenarios[i, "Halflife"],
                     "MaxKillingRate" = scenarios[i, "MaxKillingRate"],
                     "Slope" = scenarios[i, "Slope"],
                     "t" = 0:90)
  temp$PK <- PK(temp$t, C_0 = C_0, halflife = scenarios[i, "Halflife"])
  temp$PD <- PD(temp$PK, Emax = scenarios[i, "MaxKillingRate"], IC50 = IC50, n = scenarios[i, "Slope"])
  
  df_PKPD <- rbind(df_PKPD, temp)  
}

# Calculate max and min concentration and effect across all scenarios

df_PKPD <- df_PKPD %>%
  group_by(t) %>%
  summarise(PK_max = max(PK),
            PK_min = min(PK),
            PD_max = max(PD),
            PD_min = min(PD))

df_PKPD$drug <- "NEXT-GENERATION SMC"


# ----------------------------------------------------------
# Generate PK/PD profiles for SP_AQ
# ----------------------------------------------------------

df_SPAQ <- data.frame()

for(i in c("min", "max")) {
  temp <- data.frame("t" = 0:90, 
                     "id" = i,
                     "drug" = "SP-AQ")
  
  temp$PK <- PK(temp$t, 
                C_0 = C_0, 
                halflife = nextgen[nextgen$parameter == "Halflife", i])
  
  temp$PD <- PD(temp$PK, 
                Emax = nextgen[nextgen$parameter == "MaxKillingRate", i], 
                IC50 = IC50, 
                n = nextgen[nextgen$parameter == "Slope", i])
  
  df_SPAQ <- rbind(df_SPAQ, temp)  
}


df_SPAQ <- data.frame()

for(i in c("min", "max")) {
  temp <- data.frame("t" = 0:90, 
                     "id" = i,
                     "drug" = "SP-AQ")
  
  temp$PK <- PK(temp$t, 
                C_0 = C_0, 
                halflife = nextgen[nextgen$parameter == "Halflife", i])
  
  temp$PD <- PD(temp$PK, 
                Emax = nextgen[nextgen$parameter == "MaxKillingRate", i], 
                IC50 = IC50, 
                n = nextgen[nextgen$parameter == "Slope", i])
  
  df_SPAQ <- rbind(df_SPAQ, temp)  
}

# Combine into data frame for plotting
df_SPAQ <- df_SPAQ %>%
  pivot_wider(id_cols = c(id, drug, t),
              names_from = c(id),
              values_from = c(PK, PD))

df_plot_PKPD <- rbind(df_PKPD, df_SPAQ)
df_plot_PKPD$drug <- factor(df_plot_PKPD$drug, levels = c("SP-AQ", "NEXT-GENERATION SMC"))

# Convert days to weeks
df_plot_PKPD$week <- df_plot_PKPD$t / 7



######################
# WRITE DATA TO FILE #
######################

data <- list("df" = df,
             "df_plot_pe" = df_plot_pe,
             "df_plot_PKPD" = df_plot_PKPD,
             "cutoff" = cutoff)
saveRDS(data, file = "data_and_visualisation/Manuscript_Figure1/data_fig1_panelA.rds")
