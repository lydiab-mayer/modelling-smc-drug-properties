############################################################
#
# Plotting script to demonstrate model profile for SMC with
# dominant blood stage activity
#
# Written by Lydia Braunack-Mayer
# October 2022
#
############################################################

# ----------------------------------------------------------
# Set up
# ----------------------------------------------------------

# Load required packages
require(ggplot2)
require(patchwork)

# Set working directory
setwd(paste0("/scicore/home/penny/brauna0000/M3TPP/SMC_TPP/"))

# Read data
data <- readRDS("data_and_visualisation/Manuscript_Figure1/data_fig1_panelB.rds")
df <- data[["df"]]
df_plot_pe <- data[["df_plot_pe"]]
df_plot_PKPD <- data[["df_plot_PKPD"]]
cutoff <- data[["cutoff"]]

# Define plot settings
cols <- c("#C93312", "#1B4D79", "#899DA4", "#DC863B")

# ----------------------------------------------------------
# Plot trial validation
# ----------------------------------------------------------

p <- ggplot(data = df_plot_pe[df_plot_pe$drug == "Next-gen SMC", ],
            aes(x = weeks, y = mean, group = interaction(Halflife, drug))) +
  geom_line(colour = "#899DA4", linetype = "solid", alpha = 0.4, size = 0.2) 

p <- p + geom_line(data = df[df$RSS <= cutoff, ], 
                   aes(x = weeks, y = y, group = interaction(Halflife, drug)),
                   colour = cols[1], alpha = 0.8, size = 0.2, linetype = "dashed")
#p <- p + geom_ribbon(data = df_plot_pe[df_plot_pe$drug %in% c("SP-AQ"), ],
#                     aes(ymin = cl, ymax = cu, fill = drug, colour = drug), 
#                     alpha = 0.4, linetype = "dashed")

p <- p + geom_line(data = df_plot_pe[df_plot_pe$drug %in% c("SP-AQ"), ],
                   aes(x = weeks, y = mean), size = 1)

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
                            expand = expansion(mult = .03, add = 0),
                            limits = c(0, 8)) +
  scale_y_continuous(breaks = seq(0, 1, by = 0.2),
                     labels = paste0(seq(0, 100, by = 20), "%"),
                     expand = expansion(mult = .03, add = 0))

p <- p + labs(x = "WEEKS  AFTER  FINAL  SMC  CYCLE", 
              y = "PROTECTIVE  EFFICACY")


# ----------------------------------------------------------
# Plot PK/PD profiles
# ----------------------------------------------------------


q <- ggplot(df_plot_PKPD, aes(x = week, ymin = PD_min, ymax = PD_max, fill = drug, colour = drug)) +
  geom_ribbon(alpha = 0.4)
# q <- ggplot(df_plot[df_plot$drug == "NEXT-GENERATION SMC", ], aes(x = t, ymin = PD_min, ymax = PD_max)) +
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
               plot.title = element_text(hjust = 0.5, face = "bold"),
               legend.title = element_blank(),
               legend.position = "bottom")

q <- q + scale_x_continuous(breaks = seq(0, 12, by = 1),
                            limits = c(0, 8),
                            expand = expansion(mult = .03, add = 0)) +
  scale_y_continuous(breaks = seq(0, 30, by = 5),
                     expand = expansion(mult = .03, add = 0)) +
  scale_colour_manual(values = cols[c(1, 3)]) +
  scale_fill_manual(values = cols[c(1, 3)])

q <- q + labs(x = "WEEKS  AFTER  ONE  SMC  CYCLE",  
              y = "KILLING  RATE")

p + q + plot_annotation(title = "B. Next-generation SMC with dominant blood stage activity and initial, complete liver stage clearance") +
  plot_layout(guides = "collect") &
  theme(plot.title = element_text(family = "Times New Roman", face = "bold", size = 10),
        legend.position = "none")

ggsave(filename = paste0("./data_and_visualisation/Manuscript_Figure1/fig1_panelB.jpg"),
       plot = last_plot(),
       width = 9.1,
       height = 2.5,
       dpi = 400)
