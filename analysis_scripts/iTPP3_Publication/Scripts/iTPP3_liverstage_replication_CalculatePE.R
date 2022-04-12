############################################################
#
# Plotting script to demonstrate model profile for SMC with
# dominant blood stage activity
#
# Written by Lydia Braunack-Mayer
# April 2022
#
############################################################

# ----------------------------------------------------------
# Setup
# ----------------------------------------------------------

# Clear workspace
rm(list = ls())

# !!! Insert your experiment name here as a string, e.g. "MyExperiment" !!!
exp <- "iTPP3_liverstage_replication"

# Load required libraries
library(ggplot2)
library(dplyr)
library(patchwork)

# Define user and group
user <- strsplit(getwd(), "/", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]][5]
GROUP_dr <- "/scicore/home/penny/GROUP/M3TPP/"

# Set working directory
setwd(paste0("/scicore/home/penny/",user,"/M3TPP"))

# Create folder to store outputs
if (!dir.exists(paste0("./Experiments/",exp,"/Outputs"))) dir.create(paste0("./Experiments/",exp,"/Outputs"))

# # Import parameter table
# param_table <- read.table(paste0(GROUP_dr,exp,"/param_tab.txt"), sep = "\t", header = TRUE, as.is = TRUE, stringsAsFactors = FALSE)

cols <- c("#C93312", "#1B4D79", "#899DA4", "#DC863B")


###########################
# PLOT VALIDATION RESULTS #
###########################

# ----------------------------------------------------------
# Import data
# ----------------------------------------------------------

# setting <- Sys.glob(paste0(GROUP_dr, exp, "/postprocessing/seeds*"))
# 
# df <- data.frame()
# 
# for (i in setting) {
#   temp <- read.table(i, header = TRUE, sep = "")
#   df <- rbind(df, temp)
#   remove(temp)
# }
# 
# 
# # ----------------------------------------------------------
# # Calculate protective efficacy
# # ----------------------------------------------------------
# 
# # Index so that outputs for intervention/control groups for a given experiment are consecutive
# df <- df[order(df$EIR, df$IC50, df$Halflife, df$MaxKillingRate, df$Slope, df$seed), ]
# 
# # Calculate pe between consecutive rows
# df <- df %>%
#   group_by(EIR, IC50, Halflife, MaxKillingRate, Slope, seed) %>%
#   summarise(across(contains("cpp"), ~ 1 - tail(.x, 1)/head(.x, 1), .names = "pe_{col}"))
# 
# # Aggregate over seeds
# df <- df %>%
#   group_by(EIR, IC50, Halflife, MaxKillingRate, Slope) %>%
#   summarise(across(contains("cpp"), mean))
# 
# 
# # ----------------------------------------------------------
# # Import zongo trial data
# # ----------------------------------------------------------
# 
# # Import protective efficacy profile from Zongo et al. 2015
# zongo <- read.csv(paste0("./Experiments/",exp,"/analysis_scripts/zongo_data_extraction.csv"), sep = ";")
# 
# # Add dummy columns to match columns in simulation database
# df$linetype <- "mean"
# df$drug <- "Next-gen SMC"
# df$EIR <- as.numeric(df$EIR)
# zongo$EIR <- zongo$IC50 <- zongo$Halflife <- zongo$MaxKillingRate <- zongo$Slope <- 0
# zongo$drug <- as.character(zongo$drug); zongo$linetype <- as.character(zongo$linetype)
# 
# # Prepare database
# df_plot <- pivot_longer(df,
#                         cols = starts_with("pe"),
#                         names_to = "x",
#                         values_to = "y")
# df_plot$x <- as.numeric(sub("pe_cpp_", "", df_plot$x))
# df_plot$x <- df_plot$x - min(df_plot$x)
# 
# # Merge with zongo trial data
# df_plot <- rbind(zongo, as.data.frame(df_plot))
# df_plot <- pivot_wider(df_plot,
#                        id_cols = c(drug, x, Slope, MaxKillingRate, Halflife, IC50, EIR),
#                        names_from = linetype,
#                        values_from = y)
# df_plot$drug <- factor(df_plot$drug, levels = c("SP-AQ", "DHA-PPQ", "Next-gen SMC"))
# 
# # Add column for weeks
# df_plot$weeks <- round(df_plot$x*5/7, 1)

df_plot <- data.frame("drug" = c("SP-AQ", "DHA-PPQ", "NEXT-GENERATION SMC"),
                      "weeks" = rep(NA, 3),
                      "mean" = rep(NA, 3),
                      "Halflife" = rep(NA, 3),
                      "cl" = rep(NA, 3),
                      "cu" = rep(NA, 3))

# ----------------------------------------------------------
# Create plot
# ----------------------------------------------------------


p <- ggplot(data = df_plot[df_plot$drug == "Next-gen SMC", ],
            aes(x = weeks, y = mean, group = interaction(Halflife, drug))) #+
 # geom_line(colour = "#899DA4", linetype = "solid", alpha = 0.4, size = 0.1) 

#p <- p + geom_ribbon(data = df_plot[df_plot$drug %in% c("SP-AQ", "DHA-PPQ"), ],
#                     aes(ymin = cl, ymax = cu, fill = drug, colour = drug), alpha = 0.4)

p <- p + theme(panel.border = element_blank(), 
               panel.background = element_blank(),
               panel.grid = element_blank(),
               text = element_text(family = "Times New Roman", size = 10),
               strip.background = element_blank(),
               axis.line = element_blank(),
               axis.ticks = element_blank(),
               axis.title.x = element_text(margin = margin(t = 10)),
               axis.title.y = element_text(margin = margin(r = 10)),
               plot.title = element_text(hjust = 0.5, face = "bold"),
               legend.position = "none") +
  scale_colour_manual(values = cols) +
  scale_fill_manual(values = cols)

p <- p + scale_x_continuous(breaks = 0:8,
                            expand = expansion(mult = .03, add = 0)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2),
                     labels = paste0(seq(0, 100, by = 20), "%"),
                     expand = expansion(mult = .03, add = 0)) +
  coord_cartesian(ylim = c(0, 1), xlim = c(0, 8))
  
p <- p + labs(x = "WEEKS  AFTER  FINAL  SMC  ROUND", 
              y = "PROTECTIVE  EFFICACY", 
              title = "MODEL  VALIDATION  TO  EXISTING  SMC")

p



################################
# PLOT RANGE OF PK/PD PROFILES #
################################

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

df <- data.frame("t" = 0:90)
df$max <- decay(df$t, 60, 1)
df$min <- decay(df$t, 10, 0.8)
df$drug <- "NEXT-GENERATION SMC"


# ----------------------------------------------------------
# Generate database of PK/PD profiles that replicate Zongo et al trial results
# ----------------------------------------------------------

# Identify max and min PK/PD parameters
Halflife_SPAQ <- c("min" = 30, "max" = 35)
Efficacy_SPAQ <- c("min" = 0.95, "max" = 1)

Halflife_DHAPPQ <- c("min" = 25, "max" = 30)
Efficacy_DHAPPQ <- c("min" = 0.9, "max" = 1)

df_zongo <- rbind(Halflife_SPAQ, Halflife_DHAPPQ,
                  Efficacy_SPAQ, Efficacy_DHAPPQ)
df_zongo <- as.data.frame(df_zongo)

df_zongo$parameter = rep(c("Halflife", "Efficacy"), each = 2)
df_zongo$drug <- rep(c("SP-AQ", "DHA-PPQ"), 2)
rownames(df_zongo) <- 1:nrow(df_zongo)

# Generate PK/PD profiles for SP-AQ
df_SPAQ <- data.frame()

for(i in c("min", "max")) {
  temp <- data.frame("t" = 0:90, 
                     "id" = i,
                     "drug" = "SP-AQ")
  
  temp$decay <- decay(t = temp$t, 
                      halflife = df_zongo[df_zongo$drug == "SP-AQ" & df_zongo$parameter == "Halflife", i],
                      efficacy = df_zongo[df_zongo$drug == "SP-AQ" & df_zongo$parameter == "Efficacy", i])

  df_SPAQ <- rbind(df_SPAQ, temp)  
}

# Generate PK/PD profiles for DHA-PPQ

df_DHAPPQ <- data.frame()

for(i in c("min", "max")) {
  temp <- data.frame("t" = 0:90, 
                     "id" = i,
                     "drug" = "DHA-PPQ")
  
  temp$decay <- decay(t = temp$t, 
                      halflife = df_zongo[df_zongo$drug == "DHA-PPQ" & df_zongo$parameter == "Halflife", i],
                      efficacy = df_zongo[df_zongo$drug == "DHA-PPQ" & df_zongo$parameter == "Efficacy", i])
  
  df_DHAPPQ <- rbind(df_DHAPPQ, temp)  
}

# Combine into data frame for plotting
df_zongo_plot <- rbind(df_SPAQ, df_DHAPPQ)
df_zongo_plot <- df_zongo_plot %>%
  pivot_wider(id_cols = c(id, drug, t),
              names_from = id,
              values_from = decay)

df_plot <- rbind(df, df_zongo_plot)
df_plot$drug <- factor(df_plot$drug, levels = c("SP-AQ", "DHA-PPQ", "NEXT-GENERATION SMC"))

# ----------------------------------------------------------
# Generate plot of range of PK/PD profiles
# ----------------------------------------------------------

q <- ggplot(df_plot, aes(x = t, ymin = min, ymax = max, fill = drug, colour = drug)) +
  geom_ribbon(alpha = 0.4)

q <- q + theme(panel.border = element_blank(), 
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
              legend.title = element_blank(),
              legend.position = "bottom")

q <- q + scale_x_continuous(breaks = seq(0, 90, by = 30),
                            expand = expansion(mult = .03, add = 0)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2),
                     labels = paste0(seq(0, 100, by = 20), "%"),
                     expand = expansion(mult = .03, add = 0)) +
  scale_fill_manual(values = cols) +
  scale_colour_manual(values = cols)
  
q <- q + labs(x = "TIME  (DAYS)", 
              y = "PROBABILITY  OF\nPREVENTING  INFECTION", 
              title = "MODELLED  RANGE  OF  PK/PD  PROFILES")

q

p / q + plot_annotation(title = "B. Next-generation SMC with dominant liver stage activity") &
  theme(plot.title = element_text(family = "Times New Roman", face = "bold", size = 10))

ggsave(filename = paste0("./analysisworkflow/analysis_scripts/iTPP3_Publication/Figures/iTPP3_fig1_panelB.jpg"),
       plot = last_plot(),
       width = 4.5,
       height = 5,
       dpi = 400)
