############################################################
#
# Generates table of cut-off criteria for all experiments
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

# Load data
df <- readRDS("./data_and_visualisation/Manuscript_Figure4/data_fig4_panelA.rds")


# ----------------------------------------------------------
# Define plot settings
# ----------------------------------------------------------

# Define plot colours
cols <- c("#0f1f2f", "#1d3d5d", "#295784", "#3e6897", "#537bac", "#938998", "#ca9574",
          "#f9a24b", "#fab464", "#fac67c", "#fad694", "#fcdfaa", "#fee8c0")
text_size <- 10

  
# ----------------------------------------------------------
# Generate plot
# ----------------------------------------------------------

# Generate plot for program reach

p <- ggplot(df, aes(x = annual_prev, y = Target, fill = CoverageAllRounds))

p <- p + geom_tile(colour = "white", size = 1.5)

p <- p + theme(panel.border = element_blank(), 
               panel.background = element_blank(),
               panel.grid = element_blank(),
               text = element_text(family = "Times New Roman", size = text_size),
               axis.line = element_blank(),
               axis.ticks = element_blank(),
               axis.text.x = element_text(margin = margin(t = 0)),
               axis.text.y = element_text(margin = margin(r = 0)),
               legend.title = element_text(face = "bold", size = text_size),
               legend.title.align = 0.5,
               legend.position = "bottom",
               legend.key = element_blank())

p <- p + scale_fill_manual(values = cols[c(1, 2, 4, 6, 7, 8, 10, 12)],
                           na.value = "light grey",
                           labels = c(levels(df$CoverageAllRounds), "Target not met in\nparameter space")) +
  scale_y_discrete(labels = c("50%", "", "60%", "", "70%", "", "80%", "", "90%"))

p <- p + labs(y = "TARGET  REDUCTION", x = expression(paste("BASELINE  ANNUAL  ", italic("Pf"), "PR"["2-10"]))) +
  guides(fill = guide_legend(title.position = "top", title.hjust = 0.5, title = "COVERAGE  OF  ALL  SMC  CYCLES", nrow = 2))


# Generate plot for round coverage

q <- ggplot(df, aes(x = annual_prev, y = Target, fill = CoverageOneRound))

q <- q + geom_tile(colour = "white", size = 1.5)

q <- q + theme(panel.border = element_blank(), 
               panel.background = element_blank(),
               panel.grid = element_blank(),
               text = element_text(family = "Times New Roman", size = text_size),
               axis.line = element_blank(),
               axis.ticks = element_blank(),
               axis.text.x = element_text(margin = margin(t = 0)),
               axis.text.y = element_text(margin = margin(r = 0)),
               legend.title = element_text(face = "bold", size = text_size),
               legend.title.align = 0.5,
               legend.position = "bottom",
               legend.key = element_blank())

q <- q + scale_fill_manual(values = cols[c(1, 2, 4, 6, 8, 12)],
                           na.value = "light grey",
                           labels = c(levels(df$CoverageOneRound), "Target not met in\nparameter space")) +
  scale_y_discrete(labels = c("50%", "", "60%", "", "70%", "", "80%", "", "90%"))

q <- q + labs(y = "", x = expression(paste("BASELINE  ANNUAL  ", italic("Pf"), "PR"["2-10"]))) +
  guides(fill = guide_legend(title.position = "top", title.hjust = 0.5, title = "COVERAGE OF AT LEAST ONE SMC CYCLE", nrow = 2))

# Arrange figure panels

p + q + plot_annotation(title = "A.  MINIMUM  COVERAGE  CRITERIA") & 
                 theme(plot.title = element_text(family = "Times New Roman", face = "bold", size = text_size))

ggsave(filename = paste0("./data_and_visualisation/Manuscript_Figure4/fig4_panelA.jpg"),
       plot = last_plot(),
       width = 9.1,
       height = 4,
       dpi = 300)
