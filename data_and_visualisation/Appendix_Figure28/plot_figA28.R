############################################################
#
# Visualises outputs of sensitivity analysis for SMC with 
# dominant liver stage activity
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
df <- readRDS("./data_and_visualisation/Appendix_Figure28/data_figA28.rds")


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

df_plot <- df[df$Agegroup == "CHILDREN 3 TO 59 MONTHS" & df$EIR == 8, ]
df_plot <- df_plot[df_plot$Outcome %in% c("CLINICAL INCIDENCE", "PREVALENCE", "SEVERE DISEASE"), ]
df_plot <- df_plot[df_plot$parameter != "Slope [6 - 6]", ]

p <- ggplot(df_plot, aes(x = Setting, y = T_eff_scaled, fill = parameter, label = label))

p <- p + geom_bar(position = "stack", stat = "identity", colour = "white")

p <- p + geom_text(aes(colour = parameter), position = position_stack(vjust = 0.5), 
                   family = "Times New Roman",
                   size = fontsize*0.28,
                   show.legend = FALSE)

p <- p + facet_wrap(.~ Outcome, scales = "free_x", ncol = 1)

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

p <- p + scale_fill_manual(breaks =,
                           values = cols) + 
  scale_colour_manual(values = text_cols) +
  scale_y_continuous(breaks = seq(0, 100, 10), labels = paste0(seq(0, 100, 10), "%"))

p <- p + labs(x = "SCENARIO",
              y = "MEDIAN  REDUCTION  (%)",
              fill = "")

p <- p + guides(fill = guide_legend(nrow = 2)) 

p

ggsave(filename = paste0("./data_and_visualisation/Appendix_Figure28/figA28.jpg"),
       plot = last_plot(),
       width = 6,
       height = 7,
       dpi = 400)
