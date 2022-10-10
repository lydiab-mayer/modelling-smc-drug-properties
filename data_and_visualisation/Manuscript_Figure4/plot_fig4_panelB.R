############################################################
#
# Generates figures cut-off criteria by coverage
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
setwd(paste0("/scicore/home/penny/brauna0000/M3TPP/SMC_TPP"))

# Load data
df <- readRDS("./data_and_visualisation/Manuscript_Figure4/data_fig4_panelB.rds")


# ----------------------------------------------------------
# Define plot settings
# ----------------------------------------------------------

# Define plot colours
cols <- c("#0f1f2f", "#1d3d5d", "#295784", "#3e6897", "#537bac", "#938998", "#ca9574",
          "#f9a24b", "#fab464", "#fac67c", "#fad694", "#fcdfaa", "#fee8c0")

#Set text size
text_size <- 10


# ----------------------------------------------------------
# Generate plot
# ----------------------------------------------------------

p <- ggplot(df, aes(x = annual_prev, y = Target, fill = Halflife))

p <- p + geom_tile(colour = "white", size = 1.5)

p <- p + facet_grid(Coverage2 ~ Coverage1)

p <- p + theme(panel.border = element_blank(), 
               panel.background = element_blank(),
               panel.grid = element_blank(),
               text = element_text(family = "Times New Roman", size = text_size),
               strip.background = element_blank(),
               strip.text = element_text(face = "bold"),
               axis.line = element_blank(),
               axis.ticks = element_blank(),
               axis.text.x = element_text(margin = margin(t = 0)),
               axis.text.y = element_text(margin = margin(r = 0)),
               legend.title = element_text(face = "bold", size = text_size),
               legend.title.align = 0.5,
               legend.position = "bottom",
               legend.key = element_blank())

p <- p + scale_fill_manual(values = cols[c(4, 6, 8, 12)], 
                    na.value = "light grey",
                    labels = c(levels(df$Halflife), 
                               "Target not met in\nparameter space")) +
  scale_y_discrete(labels = c("50%", "", "60%", "", "70%", "", "80%", "", "90%"))

p <- p + labs(y = "TARGET  REDUCTION", x = expression(paste("BASELINE ANNUAL ", italic("Pf"), "PR"["2-10"]))) +
  guides(fill = guide_legend(title.position = "top", title = "ELIMINATION  HALF-LIFE  (DAYS)", nrow = 1))

p + plot_annotation(title = "B.  MINIMUM  ELIMINATION  HALF-LIFE  CRITERIA") & 
  theme(plot.title = element_text(family = "Times New Roman", face = "bold", size = text_size))

ggsave(filename = paste0("./data_and_visualisation/Manuscript_Figure4/plot_fig4_panelB.jpg"),
       plot = last_plot(),
       width = 9.1,
       height = 6,
       dpi = 300)


