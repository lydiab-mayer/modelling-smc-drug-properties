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
setwd(paste0("/scicore/home/penny/",user,"/M3TPP"))

# Define colours
cols <- c("#C93312", "#1B4D79", "#899DA4", "#DC863B")


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
zongo <- read.csv("./analysisworkflow/analysis_scripts/iTPP3_Publication/zongo_data_extraction.csv", sep = ";", as.is = TRUE)

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
write.csv(nextgen, paste0("./analysisworkflow/analysis_scripts/iTPP3_Publication/Figures/", exp, "_zongo_best_fit.csv"))


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


##################
# GENERATE PLOTS #
##################

# ----------------------------------------------------------
# Plot trial validation
# ----------------------------------------------------------

p <- ggplot(data = df_plot_pe[df_plot_pe$drug == "Next-gen SMC", ],
            aes(x = weeks, y = mean, group = interaction(Halflife, Efficacy, drug))) +
  geom_line(colour = "#899DA4", linetype = "solid", alpha = 0.4, size = 0.2) 

p <- p + geom_line(data = df[df$RSS <= cutoff, ], 
                   aes(x = weeks, y = y, group = interaction(Halflife, Efficacy, drug)),
                   colour = cols[1], alpha = 0.8, size = 0.2, linetype = "dashed")
#p <- p + geom_ribbon(data = df_plot_pe[df_plot_pe$drug %in% c("SP-AQ"), ],
#                     aes(ymin = cl, ymax = cu, fill = drug, colour = drug), 
#                     alpha = 0.4, linetype = "dashed")

p <- p + geom_line(data = df_plot_pe[df_plot_pe$drug %in% c("SP-AQ"), ],
                   aes(x = weeks, y = mean), colour = "#781e0b", size = 1)

p <- p + theme(panel.border = element_blank(), 
               panel.background = element_blank(),
               panel.grid = element_blank(),
               text = element_text(family = "Times New Roman", size = 10),
               strip.background = element_blank(),
               axis.line = element_blank(),
               axis.ticks = element_blank(),
               axis.title.x = element_text(margin = margin(t = 10)),
               axis.title.y = element_text(margin = margin(r = 10)),
               plot.title = element_blank(),
               legend.position = "none") +
  scale_colour_manual(values = cols) +
  scale_fill_manual(values = cols)

p <- p + scale_x_continuous(breaks = 0:8,
                            expand = expansion(mult = .03, add = 0),
                            limits = c(0, 8)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2),
                     labels = paste0(seq(0, 100, by = 20), "%"),
                     expand = expansion(mult = .03, add = 0))

p <- p + labs(x = "WEEKS  AFTER  FINAL  SMC  ROUND", 
              y = "PROTECTIVE  EFFICACY")

p


# ----------------------------------------------------------
# Plot PK/PD profiles
# ----------------------------------------------------------

q <- ggplot(df_plot_PKPD, aes(x = week, ymin = min, ymax = max, fill = drug, colour = drug)) +
  geom_ribbon(alpha = 0.4)
# q <- ggplot(df_plot[df_plot$drug == "NEXT-GENERATION SMC", ], aes(x = t, ymin = min, ymax = max)) +
#   geom_ribbon(alpha = 0.4, fill = cols[3], colour = cols[3]) +
#   geom_ribbon(data = df_plot[df_plot$drug == "SP-AQ", ], alpha = 0.5, fill = cols[1], colour = cols[1])

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
              plot.title = element_blank(),
              legend.title = element_blank(),
              legend.position = "bottom")

q <- q + scale_x_continuous(breaks = seq(0, 12, by = 1),
                            limits = c(0, 8),
                            expand = expansion(mult = .03, add = 0)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2),
                     labels = paste0(seq(0, 100, by = 20), "%"),
                     expand = expansion(mult = .03, add = 0)) +
  scale_colour_manual(values = cols[c(1, 3)]) +
  scale_fill_manual(values = cols[c(1, 3)])
  
q <- q + labs(x = "WEEKS  AFTER  ONE  SMC  ROUND",  
              y = "PROBABILITY  OF\nPREVENTING  INFECTION")

q

p + q + plot_annotation(title = "C. Next-generation SMC with dominant liver stage activity") +
  plot_layout(guides = "collect") &
  theme(plot.title = element_text(family = "Times New Roman", face = "bold", size = 10),
        legend.position = "bottom")

ggsave(filename = paste0("./analysisworkflow/analysis_scripts/iTPP3_Publication/Figures/fig1_panelC.jpg"),
       plot = last_plot(),
       width = 9.1,
       height = 3,
       dpi = 400)
