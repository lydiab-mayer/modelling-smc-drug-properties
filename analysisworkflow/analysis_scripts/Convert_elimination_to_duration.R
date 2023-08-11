############################################################
#
# Sscript to convert parameter values for elimination 
# half-life to duration
#
# Written by Lydia Braunack-Mayer
# July 2023
#
############################################################

# ----------------------------------------------------------
# Set up
# ----------------------------------------------------------

library("dplyr")

# Define model PK/PD components
PK <- function(t, C_0, halflife) {
  # t = time
  # C_0 = concentration at time 0
  # halflife = Time t at which half the initial concentration remains
  
  C_0 * exp(-(log(2)/halflife)*t)
}

PD <- function(C, Emax, IC50, n, NegConcentration) {
  # C = Concentration
  # Emax = Maximum parasite killing rate
  # IC50 = Drug concentration at which half the maximum killing rate is achieved
  # n = slope
  # NegConcentration = concentration below which effect is 0
  
  PD <- Emax * C^n / (C^n + IC50^n)
  PD[C <= NegConcentration] <- 0

  return(PD)
}

# Define fixed parameter values
dose <- 18 # dose in mg/kg
weight <- 20 # weight of child receiving drug in kg - set for demonstration purposes
vol <- 173 # volume of distribution in l/kg
C_0 <- (dose*weight) / (vol*weight)
IC50 <- 0.020831339
NegConcentration <- 0.005 # concentration below which zero effect, in mg/l


# ----------------------------------------------------------
# Generate complete database of PK/PD profiles
# ----------------------------------------------------------

# Create grid of scenarios

ngrid <- c(19, 8, 1)
grid_ranges_cont <- rbind(Halflife = c(2, 20),
                          MaxKillingRate = c(2, 30),
                          Slope = 6)

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
  temp$PD <- PD(temp$PK, Emax = scenarios[i, "MaxKillingRate"], IC50 = IC50, n = scenarios[i, "Slope"], NegConcentration = NegConcentration)
  
  df_PKPD <- rbind(df_PKPD, temp)  
}


# ----------------------------------------------------------
# Calculate duration
# ----------------------------------------------------------

# Add calcuations for duration and duration half-life
duration <- df_PKPD %>%
  group_by(Halflife, MaxKillingRate, Slope) %>%
  mutate(Duration = max(t[PD >= MaxKillingRate*0.01]),
         DurationHalflife = max(t[PD >= MaxKillingRate/2]))
duration <- unique(duration[, c("Halflife", "Duration", "DurationHalflife")])

# Write to file
write.csv(duration, "./analysisworkflow/analysis_scripts/elimination_to_duration.csv")
