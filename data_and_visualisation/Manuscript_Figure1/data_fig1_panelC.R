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
exp <- "iTPP3_ChemoLiver_TreatLiverBlood_replication"

# Load required libraries
library(ggplot2)
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
df <- df[order(df$EIR, df$Halflife, df$Efficacy, df$seed), ]

# Calculate pe between consecutive rows
df <- df %>%
  group_by(EIR, Halflife, Efficacy, seed) %>%
  summarise(across(contains("cpp"), ~ 1 - tail(.x, 1)/head(.x, 1), .names = "pe_{col}"))

# Aggregate over seeds
df <- df %>%
  group_by(EIR, Halflife, Efficacy) %>%
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
zongo$EIR <- zongo$Halflife <- zongo$Efficacy <- 0


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
  df[i, "RSS"] <- RSS(SPAQ, as.numeric(df[i, 4:15]))
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
            Efficacy_max = max(Efficacy),
            Efficacy_min = min(Efficacy))

# Format resulting dataframe
nextgen <- as.data.frame(t(nextgen))
nextgen$parameter <- gsub("_.*", "", rownames(nextgen))
nextgen$type <- gsub(".*_", "", rownames(nextgen))
names(nextgen)[1] <- "value"
nextgen <- nextgen %>%
  pivot_wider(names_from = type) %>%
  as.data.frame()

# Save results
write.csv(nextgen, paste0("./data_and_visualisation/Manuscript_Figure1/table2_modelC.csv"))


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
df_plot_pe <- pivot_wider(df_plot_pe, id_cols = c(drug, x, Efficacy, Halflife, EIR), names_from = linetype, values_from = y)
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
decay <- function(t, halflife, efficacy) {
  efficacy * exp( -(t/halflife)^5.4 * log(2) )
}


# ----------------------------------------------------------
# Generate complete database of PK/PD profiles
# ----------------------------------------------------------

df_PKPD <- data.frame("t" = 0:90)
df_PKPD$max <- decay(df_PKPD$t, 60, 1)
df_PKPD$min <- decay(df_PKPD$t, 10, 0.8)
df_PKPD$drug <- "NEXT-GENERATION SMC"


# ----------------------------------------------------------
# Generate database of PK/PD profiles that replicate Zongo et al trial results
# ----------------------------------------------------------

# Generate PK/PD profiles for SP-AQ
df_SPAQ <- data.frame()

for(i in c("min", "max")) {
  temp <- data.frame("t" = 0:90, 
                     "id" = i,
                     "drug" = "SP-AQ")
  
  temp$decay <- decay(t = temp$t, 
                      halflife = nextgen[nextgen$parameter == "Halflife", i],
                      efficacy = nextgen[nextgen$parameter == "Efficacy", i])

  df_SPAQ <- rbind(df_SPAQ, temp)  
}

# Combine into data frame for plotting
df_SPAQ <- df_SPAQ %>%
  pivot_wider(id_cols = c(id, drug, t),
              names_from = id,
              values_from = decay)

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
saveRDS(data, file = "data_and_visualisation/Manuscript_Figure1/data_fig1_panelC.rds")