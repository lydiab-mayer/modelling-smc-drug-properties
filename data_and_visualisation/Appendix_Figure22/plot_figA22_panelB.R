############################################################
#
# Visualises relationships between emulator input and predictor variables
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
df <- readRDS("./data_and_visualisation/Appendix_Figure22/data_figA22_panelB.rds")


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

df_plot <- df[(df$Access == "HIGH ACCESS" & df$Seasonality == "5 MONTH SEASON") & df$Agegroup == "CHILDREN 3 TO 59 MONTHS", ]
df_plot <- df_plot[df_plot$Outcome %in% c("CLINICAL INCIDENCE", "SEVERE DISEASE", "PREVALENCE"), ]

p <- ggplot(df_plot, aes(x = annual_prev_lab, y = T_eff_scaled, fill = parameter, label = label))

p <- p + geom_bar(position = "stack", stat = "identity", colour = "white")

p <- p + geom_text(aes(colour = parameter), 
                   position = position_stack(vjust = 0.5),
                   family = "Times New Roman",
                   size = fontsize*0.28,
                   show.legend = FALSE)

p <- p + facet_wrap(Outcome ~ ., scales = "free_x", ncol = 1)

p <- p + theme(panel.border = element_blank(), 
               panel.background = element_blank(),
               panel.grid = element_blank(),
               text = element_text(family = "Times New Roman", size = fontsize),
               title = element_text(face = "bold"),
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
  scale_y_continuous(breaks = seq(0, 100, 20), labels = paste0(seq(0, 100, 20), "%"))

p <- p + labs(x = expression(paste("BASELINE  ANNUAL  ", italic("Pf"), "PR"["2-10"])),
              y = "MEDIAN  REDUCTION  (%)",
              fill = "",
              title = "B")

p <- p + guides(fill = guide_legend(nrow = 2)) 

p

ggsave(filename = paste0("./data_and_visualisation/Appendix_Figure22/figA22_panelB.jpg"),
       plot = last_plot(),
       width = 4.5,
       height = 5,
       dpi = 400)


# ----------------------------------------------------------
# Generate supporting table with sensitivity results across scenarios
# ----------------------------------------------------------

tab <- df %>%
  select(-label)

temp <- tab %>%
  filter(parameter %in% c("Round coverage [70% - 95%]", "Cycle coverage [70% - 95%]")) %>%
  group_by(Seasonality, EIR, Access, Agegroup, Outcome) %>%
  mutate(S_eff = sum(S_eff),
         T_eff = sum(T_eff),
         T_eff_scaled = sum(T_eff_scaled),
         parameter = "Combined coverage")

temp <- as.data.frame(unique(temp))

tab$parameter <- as.character(tab$parameter)

tab <- rbind(tab, temp)

tab <- tab %>%
  group_by(parameter, Outcome) %>%
  summarise(max = max(S_eff), min = min(S_eff))

tab$minmax <- paste0(round(tab$min*100, 0), "% to ", round(tab$max*100, 0), "%")
names(tab) <- c("Key performance property", "Outcome", "Max", "Min", "Range of attributable outcome variation")

write.csv(tab[, c(1, 2, 5)], file = paste0("./data_and_visualisation/Appendix_Figure22/tab.csv"))