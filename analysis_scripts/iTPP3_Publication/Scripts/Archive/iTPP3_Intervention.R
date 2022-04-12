############################################################
#
# Visualises intervention dynamics for blood and liver stage SMC
#
# Written by Lydia Braunack-Mayer
############################################################

rm(list = ls())

library(ggplot2)
library(dplyr)
library(patchwork)

###################### BLOOD STAGE SMC######################

# ----------------------------------------------------------
# PK/PD components
# ----------------------------------------------------------

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


# ----------------------------------------------------------
# PK/PD Set up
# ----------------------------------------------------------

dose <- 18 # dose in mg/kg
weight <- 20 # weight of individual receiving drug in kg
vol <- 173 # volume of distribution in l/kg
C_0 <- (dose*weight) / (vol*weight)
IC50 <- 0.020831339

# halflife <- 5
# Emax <- 3.45
# slope <- 8
# 
# nby <- 0.1
# nto <- 90
# 
# PKPD <- data.frame("t" = seq(0, nto, by = nby))
# PKPD$PK <- PK(PKPD$t, C_0, halflife) # first dose
# # PKPD$PK <- PKPD$PK + c(rep(0, 1/nby), PK(PKPD$t[1:(nto/nby - 1/nby + 1)], C_0, halflife)) # second dose
# # PKPD$PK <- PKPD$PK + c(rep(0, 1/nby*2), PK(PKPD$t[1:(nto/nby - 1/nby*2 + 1)], C_0, halflife)) # third dose
# PKPD$PD <- PD(PKPD$PK, Emax, IC50, n = slope)
# 
# par(mfrow = c(1, 3))
# plot(PKPD$t, PKPD$PK, type = "l", xlab = "Time (days)", main = "PK", ylab = "Concentration")
# plot(PKPD$PK, PKPD$PD, type = "l", xlab = "Concentration", main = "PD", ylab = "Effect")
# plot(PKPD$t, PKPD$PD, type = "l", xlab = "Time (days)", main = "PK/PD", ylab = "Effect")
# par(mfrow = c(1, 1))


# ----------------------------------------------------------
# Generate database of PK/PD profiles
# ----------------------------------------------------------

ngrid <- c(4, 8, 8)
grid_ranges_cont <- rbind(Halflife = c(5, 20),
                          MaxKillingRate = c(2, 30),
                          Slope = c(1, 8))

# Create grid of scenarios
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

df <- data.frame()

for(i in 1:nrow(scenarios)) {
  temp <- data.frame("Halflife" = scenarios[i, "Halflife"],
                     "MaxKillingRate" = scenarios[i, "MaxKillingRate"],
                     "Slope" = scenarios[i, "Slope"],
                     "t" = 0:90)
  temp$PK <- PK(temp$t, C_0 = C_0, halflife = scenarios[i, "Halflife"])
  temp$PD <- PD(temp$PK, Emax = scenarios[i, "MaxKillingRate"], IC50 = IC50, n = scenarios[i, "Slope"])
  
  df <- rbind(df, temp)  
}

# Calculate max and min concentration and effect across all scenarios

df_plot <- df %>%
  group_by(t) %>%
  summarise(maxPK = max(PK),
            minPK = min(PK),
            maxPD = max(PD),
            minPD = min(PD))

# Add dummy PK and PD columns for purposing of plotting sample lines
df_plot$PK <- seq(0, max(df_plot$maxPK), length = 91)
df_plot$PD <- seq(min(df_plot$minPD), max(df_plot$maxPD), length = 91)

# Prepare sample PK/PD profiles
df_sample1 <- df[(df$Halflife == 20 & df$Slope == 8) & df$MaxKillingRate == 30, ]
df_sample2 <- df[(df$Halflife == 5 & df$Slope == 1) & df$MaxKillingRate == 30, ]
df_sample3 <- df[(df$Halflife == 10 & df$Slope == 4) & df$MaxKillingRate == 18, ]
df_sample4 <- df[(df$Halflife == 15 & df$Slope == 6) & df$MaxKillingRate == 2, ]

# ----------------------------------------------------------
# Generate plot of range of PK/PD profiles
# ----------------------------------------------------------

# PK profile
p <- ggplot(df_sample1, aes(x = t, y = PK)) +
  geom_ribbon(data = df_plot, aes(x = t, ymin = minPK, ymax = maxPK), alpha = 0.1) +
  geom_line(colour = "#C93312") +
  geom_line(data = df_sample2, colour = "#899DA4") +
  geom_line(data = df_sample3, colour = "#DC863B") +
  geom_line(data = df_sample4, colour = "#425055") +
  theme(panel.border = element_blank(), 
                 panel.background = element_blank(),
                 panel.grid.major = element_line(colour = "grey95"),
                 panel.grid.minor = element_blank(),
                 text = element_text(family = "Times New Roman", size = 10),
                 strip.background = element_blank(),
                 axis.line = element_blank(),
                 axis.ticks = element_blank(),
                 axis.title.x = element_text(margin = margin(t = 10)),
                 axis.title.y = element_text(margin = margin(r = 10)),
                 plot.title = element_text(hjust = 0.5, face = "bold"),
                 legend.title = element_text(face = "bold")) +
  scale_x_continuous(breaks = seq(0, 90, by = 10),
                     expand = expansion(mult = .03, add = 0)) +
  scale_y_continuous(expand = expansion(mult = .03, add = 0)) +
  labs(x = "TIME  (DAYS)", y = "CONCENTRATION\n(MG/L)", title = "PK")

# PD profile
q <- ggplot(df_sample1, aes(x = PK, y = PD)) +
  geom_ribbon(data = df_plot, aes(x = maxPK, ymin = minPD, ymax = maxPD), alpha = 0.1) +
  geom_line(colour = "#C93312") +
  geom_line(data = df_sample2, colour = "#899DA4") +
  geom_line(data = df_sample3, colour = "#DC863B") +
  geom_line(data = df_sample4, colour = "#425055") +
  theme(panel.border = element_blank(), 
        panel.background = element_blank(),
        panel.grid.major = element_line(colour = "grey95"),
        panel.grid.minor = element_blank(),
        text = element_text(family = "Times New Roman", size = 10),
        strip.background = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 10)),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.title = element_text(face = "bold")) +
  scale_x_continuous(expand = expansion(mult = .03, add = 0)) +
  scale_y_continuous(breaks = seq(0, 30, by = 5),
                     expand = expansion(mult = .03, add = 0)) +
  labs(x = "CONCENTRATION  (MG/L)", y = "PARASITE\nKILLING  RATE", title = "PD")

# PK/PD profile
r <- ggplot(df_sample1, aes(x = t, y = PD)) +
  geom_ribbon(data = df_plot, aes(x = t, ymin = minPD, ymax = maxPD), alpha = 0.1) +
  geom_line(colour = "#C93312") +
  geom_line(data = df_sample2, colour = "#899DA4") +
  geom_line(data = df_sample3, colour = "#DC863B") +
  geom_line(data = df_sample4, colour = "#425055") +
  theme(panel.border = element_blank(), 
        panel.background = element_blank(),
        panel.grid.major = element_line(colour = "grey95"),
        panel.grid.minor = element_blank(),
        text = element_text(family = "Times New Roman", size = 10),
        strip.background = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 10)),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.title = element_text(face = "bold")) +
  scale_x_continuous(breaks = seq(0, 90, by = 30),
                     expand = expansion(mult = .03, add = 0)) +
  scale_y_continuous(breaks = seq(0, 30, by = 5),
                     expand = expansion(mult = .03, add = 0)) +
  labs(x = "TIME  (DAYS)", y = "PARASITE\nKILLING  RATE", title = "PK/PD")

p / q / r + plot_annotation(title = "A. Next-generation SMC with blood stage activity") &
  theme(text = element_text(family = "Times New Roman", face = "bold", size = 10))

ggsave(filename = paste0("./analysisworkflow/analysis_scripts/iTPP3_Publication/Figures/iTPP3_PKPD_BS.jpg"),
       plot = last_plot(),
       width = 4,
       height = 5,
       dpi = 400)



###################### LIVER STAGE SMC######################

# ----------------------------------------------------------
# Liver stage set up
# ----------------------------------------------------------

decay <- function(t, halflife, efficacy) {
  efficacy * exp( -(t/halflife)^5.4 * log(2) )
}

# ----------------------------------------------------------
# Generate database of PK/PD profiles
# ----------------------------------------------------------

df_decay <- data.frame("t" = 0:90)
df_decay$decay <- decay(df_decay$t, 31, 1)
df_decay$decay1 <- decay(df_decay$t, 14, .85)
df_decay$decay2 <- decay(df_decay$t, 21, .9)
df_decay$decay3 <- decay(df_decay$t, 28, .95)
df_decay$decay4 <- decay(df_decay$t, 35, 1)
df_decay$decay5 <- decay(df_decay$t, 42, 1)
df_decay$decay6 <- decay(df_decay$t, 49, 1)
df_decay$decay7 <- decay(df_decay$t, 56, 1)
df_decay$decay_max <- decay(df_decay$t, 60, 1)
df_decay$decay_min <- decay(df_decay$t, 10, 0.8)

# ----------------------------------------------------------
# Generate plot of range of PK/PD profiles
# ----------------------------------------------------------

s <- ggplot(df_decay, aes(x = t, y = decay)) +
  theme(panel.border = element_blank(), 
        panel.background = element_blank(),
        panel.grid.major = element_line(colour = "grey95"),
        panel.grid.minor = element_blank(),
        text = element_text(family = "Times New Roman", size = 10),
        strip.background = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 10)),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.title = element_text(face = "bold")) +
  labs(x = "TIME  (DAYS)",
       y = "CONCENTRATION\n(MG/L)",
       title  = "PK")

t <- ggplot(df_decay, aes(x = t, y = decay)) +
  theme(panel.border = element_blank(), 
        panel.background = element_blank(),
        panel.grid.major = element_line(colour = "grey95"),
        panel.grid.minor = element_blank(),
        text = element_text(family = "Times New Roman", size = 10),
        strip.background = element_blank(),
        axis.line = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        axis.title.x = element_text(margin = margin(t = 10)),
        axis.title.y = element_text(margin = margin(r = 10)),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        legend.title = element_text(face = "bold")) +
  labs(x = "CONCENTRATION\n(MG/L)",
       y = "PROB  PREVENT\nINFECTION",
       title  = "PD")

# PK/PD profile
  
u <- ggplot(df_decay, aes(x = t, y = decay)) +
  geom_ribbon(aes(ymin = decay_min, ymax = decay_max), alpha = 0.1) +
  geom_line() +
  geom_label(label = "SP+AQ", x = 31, y = .5, fill = "black", label.size = NA, family = "Times New Roman", fontface = "bold", colour = "white", size = 2)

u <- u + theme(panel.border = element_blank(), 
               panel.background = element_blank(),
               panel.grid.major = element_line(colour = "grey95"),
               panel.grid.minor = element_blank(),
               text = element_text(family = "Times New Roman", size = 10),
               strip.background = element_blank(),
               axis.line = element_blank(),
               axis.ticks = element_blank(),
               axis.title.x = element_text(margin = margin(t = 10)),
               axis.title.y = element_text(margin = margin(r = 10)),
               plot.title = element_text(hjust = 0.5, face = "bold"),
               legend.title = element_text(face = "bold")) +
  scale_x_continuous(breaks = seq(0, 90, by = 30),
                     expand = expansion(mult = .03, add = 0)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2), 
                     labels = paste0(seq(0, 100, by = 20), "%"),
                     expand = expansion(mult = .03, add = 0)) +
  labs(x = "TIME  (DAYS)",
       y = "PROB  PREVENT\nINFECTION",
       title  = "PK/PD")

s / t / u + plot_annotation(title = "B. Next-generation SMC with liver stage activity") &
  theme(text = element_text(family = "Times New Roman", face = "bold", size = 10))

ggsave(filename = paste0("./analysisworkflow/analysis_scripts/iTPP3_Publication/Figures/iTPP3_PKPD_LS.jpg"),
       plot = last_plot(),
       width = 4,
       height = 5,
       dpi = 400)
