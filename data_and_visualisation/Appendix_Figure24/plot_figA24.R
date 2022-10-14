############################################################
#
# Visualises outputs of sensitivity analysis for SMC with 
# dominant blood stage activity
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
require(wesanderson)

# Set working directory
library(dplyr)
library(tidyr)


setwd(paste0("/scicore/home/penny/brauna0000/M3TPP/SMC_TPP/"))

# Load data
df <- readRDS("./data_and_visualisation/Appendix_Figure24/data_figA24.rds")


# ----------------------------------------------------------
# Define plot settings
# ----------------------------------------------------------

# Define colours
cols <- c(wes_palette("Royal1", n = 4)[c(2, 1, 3, 4)], "#425055")
text_cols <- c("#5f1909", "#323d42", "#827d55", "#7f4a1b", "white")

# Define y-axis labels
y_axis <- c("MEDIAN INCIDENCE REDUCTION", "MEDIAN SEVERE DISEASE REDUCTION")
names(y_axis) <-  c("inc_red_int_Tot", "sev_red_int_Tot")

fontsize <- 10


# ----------------------------------------------------------
# Generate plot
# ----------------------------------------------------------

df_plot <- df[df$Access == "HIGH ACCESS" & df$Seasonality == "5 MONTH SEASON", ]
df_plot <- df_plot[df_plot$Outcome %in% c("CLINICAL INCIDENCE", "SEVERE DISEASE", "PREVALENCE"), ]
df_plot <- df_plot[df_plot$parameter != "Slope [6 - 6]", ]

p <- ggplot(df_plot, aes(x = annual_prev_lab, y = T_eff_scaled, fill = parameter, label = label))

p <- p + geom_bar(position = "stack", stat = "identity", colour = "white")

p <- p + geom_text(aes(colour = parameter), position = position_stack(vjust = 0.5), 
                   family = "Times New Roman",
                   size = fontsize*0.28,
                   show.legend = FALSE)

p <- p + facet_grid(Outcome ~ Agegroup, scales = "free_x")

p <- p + theme(panel.border = element_blank(), 
               panel.background = element_blank(),
               panel.grid = element_blank(),
               text = element_text(family = "Times New Roman", size = 10),
               strip.background = element_blank(),
               axis.line = element_blank(),
               axis.ticks = element_blank(),
               axis.title.x = element_text(margin = margin(t = 10)),
               axis.title.y = element_text(margin = margin(r = 10)),
               strip.text = element_text(face = "bold"),
               legend.key = element_blank(),
               legend.position = "bottom")

p <- p + scale_fill_manual(values = cols) + 
  scale_colour_manual(values = text_cols) +
  scale_y_continuous(breaks = seq(0, 100, 10), labels = paste0(seq(0, 100, 10), "%"))

p <- p + labs(x = expression(paste("BASELINE  ANNUAL  ", italic("Pf"), "PR"["2-10"])),
              y = "MEDIAN  REDUCTION  (%)",
              fill = "")

p <- p + guides(fill = guide_legend(nrow = 2)) 

p

ggsave(filename = paste0("./data_and_visualisation/Appendix_Figure24/figA24.jpg"),
       plot = last_plot(),
       width = 9.1,
       height = 9.1,
       dpi = 400)
