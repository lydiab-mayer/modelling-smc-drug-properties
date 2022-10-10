############################################################
#
# Visualises emulator fit
#
# Written by Lydia Braunack-Mayer
############################################################

# ----------------------------------------------------------
# Setup
# ----------------------------------------------------------

# Load required packages
require(ggplot2)
require(patchwork)

# Set working directory
setwd(paste0("/scicore/home/penny/brauna0000/M3TPP/SMC_TPP/"))

# Load data
data <- readRDS("./data_and_visualisation/Manuscript_Figure2/data_fig2_panelA.rds")


# ----------------------------------------------------------
# Generate plot for coverage1
# ----------------------------------------------------------

col <- "#C93312"

p <- ggplot(data[["coverage1"]], aes(x = Coverage1, y = mean, ymin = cl, ymax = cu)) +
  geom_ribbon(alpha = 0.3, fill = col) +
  geom_point(colour = col)

p <- p + theme(panel.border = element_blank(), 
               panel.background = element_blank(),
               panel.grid.major = element_line(colour = "grey95"),
               panel.grid.minor = element_blank(),
               text = element_text(family = "Times New Roman", size = 10),
               strip.background = element_blank(),
               axis.line = element_blank(),
               axis.ticks = element_blank(),
               axis.text.x = element_text(margin = margin(t = 0)),
               axis.text.y = element_text(margin = margin(r = 0)),
               axis.title.x = element_text(margin = margin(t = 10)),
               axis.title.y = element_text(margin = margin(r = 10)),
               plot.title = element_text(hjust = 0.5, face = "bold"),
               legend.title = element_text(face = "bold"))

p <- p + scale_x_continuous(breaks = seq(0.7, 1.0, by= 0.1),
                            limits = c(0.7, 1.0),
                            labels = paste0(seq(0.7, 1.0, by= 0.1)*100, "%"),
                            expand = expansion(mult = .03, add = 0)) + 
  scale_y_continuous(breaks = seq(0, 100, by = 10),
                     limits = c(0, 100),
                     labels = paste0(seq(0, 100, by = 10), "%"),
                     expand = expansion(mult = .03, add = 0))

p <- p + labs(x = "PROGRAM  REACH", y = "CLINICAL  INCIDENCE\nREDUCTION")

p_coverage1 <- p
  

# ----------------------------------------------------------
# Generate figure for coverage2
# ----------------------------------------------------------

col_fill <- "#FAEFD1"
col_colour <- "#827d55"

p <- ggplot(data[["coverage2"]], aes(x = Coverage2, y = mean, ymin = cl, ymax = cu)) +
  geom_ribbon(alpha = 0.8, fill = col_fill) +
  geom_point(colour = col_colour)

p <- p + theme(panel.border = element_blank(), 
               panel.background = element_blank(),
               panel.grid.major = element_line(colour = "grey95"),
               panel.grid.minor = element_blank(),
               text = element_text(family = "Times New Roman", size = 10),
               strip.background = element_blank(),
               axis.line = element_blank(),
               axis.ticks = element_blank(),
               axis.text.x = element_text(margin = margin(t = 0)),
               axis.text.y = element_blank(),
               axis.title.x = element_text(margin = margin(t = 10)),
               axis.title.y = element_text(margin = margin(r = 10)),
               plot.title = element_text(hjust = 0.5, face = "bold"),
               legend.title = element_text(face = "bold"))

p <- p + scale_x_continuous(breaks = seq(0.7, 1.0, by= 0.1),
                            labels = paste0(seq(0.7, 1.0, by= 0.1)*100, "%"),
                            expand = expansion(mult = .03, add = 0)) + 
  scale_y_continuous(breaks = seq(0, 100, by = 10),
                     limits = c(0, 100),
                     expand = expansion(mult = .03, add = 0))

p <- p + labs(x = "ROUND  COVERAGE", y = "")

p_coverage2 <- p
  
  
# ----------------------------------------------------------
# Generate figure for halflife
# ----------------------------------------------------------

col <- "#899DA4"

p <- ggplot(data[["halflife"]], aes(x = Halflife, y = mean, ymin = cl, ymax = cu)) +
  geom_ribbon(alpha = 0.3, fill = col) +
  geom_point(colour = col)

p <- p + theme(panel.border = element_blank(), 
               panel.background = element_blank(),
               panel.grid.major = element_line(colour = "grey95"),
               panel.grid.minor = element_blank(),
               text = element_text(family = "Times New Roman", size = 10),
               strip.background = element_blank(),
               axis.line = element_blank(),
               axis.ticks = element_blank(),
               axis.text.x = element_text(margin = margin(t = 0)),
               axis.text.y = element_text(margin = margin(r = 0)),
               axis.title.x = element_text(margin = margin(t = 10)),
               axis.title.y = element_text(margin = margin(r = 10)),
               plot.title = element_text(hjust = 0.5, face = "bold"),
               legend.title = element_text(face = "bold"))

p <- p + scale_x_continuous(breaks = seq(0, 20, by = 5),
                            limits = c(0, 20),
                            expand = expansion(mult = .03, add = 0)) + 
  scale_y_continuous(breaks = seq(0, 100, by = 10),
                     limits = c(0, 100),
                     labels = paste0(seq(0, 100, by = 10), "%"),
                     expand = expansion(mult = .03, add = 0))

p <- p + labs(x = "ELIMINATION\nHALF-LIFE  (DAYS)", y = "CLINICAL  INCIDENCE\nREDUCTION")

p_halflife <- p
  

# ----------------------------------------------------------
# Generate figure for max killing rate
# ----------------------------------------------------------

col <- "#DC863B"

p <- ggplot(data[["MaxKillingRate"]], aes(x = MaxKillingRate, y = mean, ymin = cl, ymax = cu)) +
  geom_ribbon(alpha = 0.3, fill = col) +
  geom_point(colour = col)

p <- p + theme(panel.border = element_blank(), 
               panel.background = element_blank(),
               panel.grid.major = element_line(colour = "grey95"),
               panel.grid.minor = element_blank(),
               text = element_text(family = "Times New Roman", size = 10),
               strip.background = element_blank(),
               axis.line = element_blank(),
               axis.ticks = element_blank(),
               axis.text.x = element_text(margin = margin(t = 0)),
               axis.text.y = element_blank(),
               axis.title.x = element_text(margin = margin(t = 10)),
               axis.title.y = element_text(margin = margin(r = 10)),
               plot.title = element_text(hjust = 0.5, face = "bold"),
               legend.title = element_text(face = "bold"))

p <- p + scale_x_continuous(breaks = seq(0, 30, by = 5),
                            limits = c(0, 30),
                            expand = expansion(mult = .03, add = 0)) + 
  scale_y_continuous(breaks = seq(0, 100, by = 10),
                     limits = c(0, 100),
                     expand = expansion(mult = .03, add = 0))

p <- p + labs(x = expression("E"["max"]), y = "")

p_maxkilling <- p
  
  
# ----------------------------------------------------------
# Generate figure for slope
# ----------------------------------------------------------

col <- "#425055"

p <- ggplot(data[["slope"]], aes(x = Slope, y = mean, ymin = cl, ymax = cu)) +
  geom_ribbon(alpha = 0.3, fill = col) +
  geom_point(colour = col)

p <- p + theme(panel.border = element_blank(), 
               panel.background = element_blank(),
               panel.grid.major = element_line(colour = "grey95"),
               panel.grid.minor = element_blank(),
               text = element_text(family = "Times New Roman", size = 10),
               strip.background = element_blank(),
               axis.line = element_blank(),
               axis.ticks = element_blank(),
               axis.text.x = element_text(margin = margin(t = 0)),
               axis.text.y = element_blank(),
               axis.title.x = element_text(margin = margin(t = 10)),
               axis.title.y = element_text(margin = margin(r = 10)),
               plot.title = element_text(hjust = 0.5, face = "bold"),
               legend.title = element_text(face = "bold"))

p <- p + scale_x_continuous(breaks = seq(0, 8, by = 2),
                            limits = c(0, 8),
                            expand = expansion(mult = .03, add = 0)) + 
  scale_y_continuous(breaks = seq(0, 100, by = 10),
                     limits = c(0, 100),
                     expand = expansion(mult = .03, add = 0))

p <- p + labs(x = "PD  MODEL  SLOPE", y = "")

p_slope <- p
  

# ----------------------------------------------------------
# Generate final figure
# ----------------------------------------------------------

p_out <- (p_coverage1 | p_coverage2) / (p_halflife | p_maxkilling) + plot_annotation(title = "A") & theme(plot.title = element_text(family = "Times New Roman", face = "bold"))

ggsave(filename = paste0("./data_and_visualisation/Manuscript_Figure2/fig2_panelA.jpg"),
       plot = p_out,
       width = 4.5,
       height = 5,
       dpi = 400)
