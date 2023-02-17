############################################################
#
# Visualises SMC deployment dynamics
#
# Written by Lydia Braunack-Mayer
# October 2022
#
############################################################

# ----------------------------------------------------------
# Setup
# ----------------------------------------------------------

# Load required packages
library(ggplot2)
library(patchwork)

# Set working directory
setwd(paste0("/scicore/home/penny/brauna0000/M3TPP/SMC_TPP/"))

# Read data
df <- readRDS("./data_and_visualisation/Appendix_Figure12/data_figA12.rds")

# Define plot settings
month <- 73/12
int <- month*3:7


# ----------------------------------------------------------
# Generate plot
# ----------------------------------------------------------

# Initialise plot
p <- ggplot()

# Add data for average EIR
p <- p + geom_rect(aes(xmin = 3, xmax = 8, ymin = 0, ymax = max(df$Average$simulated.EIR)), fill = "#899DA4", alpha = 0.1) +
geom_line(data = df$Average, aes(x = timestep/(73/12), y = simulated.EIR, colour = scenario_id)) + 
  geom_vline(xintercept = 3:7, colour = "#899DA4", linetype = "dashed") 

# Add labels
p <- p + geom_text(aes(x = 3, y = max(df$Average$simulated.EIR)*1.05, label = "Cycle 5"), family = "Times New Roman") +
  geom_text(aes(x = 4, y = max(df$Average$simulated.EIR)*1.05, label = "Cycle 1"), family = "Times New Roman") +
  geom_text(aes(x = 5, y = max(df$Average$simulated.EIR)*1.05, label = "Cycle 2"), family = "Times New Roman") +
  geom_text(aes(x = 6, y = max(df$Average$simulated.EIR)*1.05, label = "Cycle 3"), family = "Times New Roman") +
  geom_text(aes(x = 7, y = max(df$Average$simulated.EIR)*1.05, label = "Cycle 4"), family = "Times New Roman") +
  geom_text(aes(x = 5.5, y = 0.02, label = "Follow up period"), family = "Times New Roman")

# Add lines marking rounds

# Define plot theme
p <- p + theme(panel.border = element_blank(), 
                 panel.background = element_blank(),
                 panel.grid = element_blank(),
                 text = element_text(family = "Times New Roman", size = 12),
                 strip.background = element_blank(),
                 axis.line = element_blank(),
                 axis.ticks = element_blank(),
               axis.text.y = element_blank(),
                 axis.title.x = element_text(margin = margin(t = 12)),
                 plot.title = element_text(hjust = 0.5, face = "bold"),
                 legend.title = element_text(face = "bold"),
               legend.position = "bottom",
               legend.key = element_blank()) +
  scale_colour_manual(values = c("#DC863B", "#C93312")) +
  scale_x_continuous(breaks = seq(0, 12, by = 2),
                     expand = expansion(mult = .03, add = 0))

# Define plot titles
p <- p + labs(x = "MONTH",
              y = "ENTOMOLOGICAL\nINNOCULATION  RATE",
              colour = "SEASONAL  PROFILE")

p

# Save plot
ggsave(filename= paste0("./data_and_visualisation/Appendix_Figure12/plot_figA12.jpg"),
       plot = last_plot(),
       width = 9.1,
       height = 4,
       dpi = 400)
