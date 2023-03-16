############################################################
#
# Visualises results of optimisation procedure
#
# Written by Lydia Braunack-Mayer
# October 2022
#
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
df <- readRDS("./data_and_visualisation/Manuscript_Figure3/data_fig3.rds")

# ----------------------------------------------------------
# Define plot settings
# ----------------------------------------------------------

# Define plot colours
# cols <- c("0%" = "#010203", "5%" = "#0f1f2f", "10%" = "#152c43", "15%" = "#1d3d5d", "20%" = "#295784",
#           "25%" = "#3e6897", "30%" = "#537bac", "35%" = "#938998", "40%" = "#a6899a", "45%" = "#b98a93",
#           "50%" = "#c78d85", "55%" = "#ca9574", "60%" = "#f9a24b", "65%" = "#fab464", "70%" = "#fac67c",
#           "75%" = "#fad694", "80%" = "#fcdfaa", "85%" = "#fee8c0", "90%" = "#fff4e1")
cols <- c("0%" = "#010203", "10%" = "#0f1f2f", "20%" = "#1d3d5d", "30%" = "#295784", "40%" = "#537bac",
          "50%" = "#938998", "60%" = "#ca9574", "70%" = "#f9a24b", "80%" = "#fad694", "90%" = "#fee8c0",
          "100%" = "#fff4e1")

# Define legend titles
leg_title <- c("inc_red_int_Tot" = "CLINICAL\nINCIDENCE\nREDUCTION",
               "sev_red_int_Tot" = "SEVERE\nDISEASE\nREDUCTION",
               "prev_red_int_Aug" = "PREVALENCE\nREDUCTION",
               "mor_red_int_Tot" = "MORTALITY\nREDUCTION")


# ----------------------------------------------------------
# Generate figure
# ----------------------------------------------------------

# Define data
df_plot <- df[df$EIR == 8, ]
df_cols <- cols[names(cols) %in% unique(df_plot$target_label)]

# Generate plot
p <- ggplot(df_plot, aes(x = Halflife, y = MaxKillingRate, fill = target_label))

p <- p + geom_tile()

p <- p + facet_grid(Coverage2 ~ Coverage1)

p <- p + geom_rect(aes(xmin = 4.74, xmax = 8.81, ymin = 2.28, ymax = 29.96), colour = "white", fill = NA, size = 0.7) + #uncomment for blood stage only
#p <- p + geom_rect(aes(xmin = 5.12, xmax = 8.81, ymin = 2.28, ymax = 29.96), colour = "white", fill = NA, size = 0.7) + #uncomment for dominant blood stage
  annotate(geom = "label", x = 6.9, y = 15, fill = "white", label = "SP-AQ", size = 2.5, label.size = 0, family = "serif")

p <- p + theme(panel.border = element_blank(), 
               panel.background = element_blank(),
               panel.grid = element_blank(),
               text = element_text(family = "serif", size = 10),
               strip.background = element_blank(),
               axis.line = element_blank(),
               axis.ticks = element_blank(),
               axis.text.x = element_text(margin = margin(t = 0)),
               axis.text.y = element_text(margin = margin(r = 0)),
               axis.title.x = element_text(margin = margin(t = 10)),
               axis.title.y = element_text(margin = margin(r = 10)),
               plot.title = element_text(hjust = 0.5, face = "bold"),
               legend.title = element_text(face = "bold"),
               legend.position = "bottom")

p <- p + scale_fill_manual(values = rev(df_cols)) +
  scale_x_continuous(breaks = seq(2, 20, 2)) + 
  scale_y_continuous(breaks = seq(2, 30, 4))

p <- p + labs(x = "ELIMINATION  HALF-LIFE  (DAYS)",
              y = expression("E"["max"]))

p <- p + guides(fill = guide_legend(title = leg_title["inc_red_int_Tot"], nrow = 1, reverse = TRUE))

p

ggsave(filename = paste0("./data_and_visualisation/Manuscript_Figure3/fig3.jpg"),
       plot = last_plot(),
       width = 9.1,
       height = 7,
       dpi = 300)

ggsave(filename = paste0("./data_and_visualisation/Manuscript_Figure3/fig3.pdf"),
       plot = last_plot(),
       width = 9.1,
       height = 7,
       dpi = 300)


# Identify point predictions reported in main body of paper

df_point <- df_plot[df_plot$Halflife == 10 & df_plot$MaxKillingRate == 10, ]
df_point[df_point$Coverage1 == "85% ROUND COVERAGE" & df_point$Coverage2 == "85% CYCLE COVERAGE", ]
df_point[df_point$Coverage1 == "95% ROUND COVERAGE" & df_point$Coverage2 == "85% CYCLE COVERAGE", ]

