############################################################
#
# Compares outcomes across endpoints for SP+AQ
# Note: This script depends on outputs of the script 5_optimisation_workflow.R
#
# Written by Lydia Braunack-Mayer
############################################################

# ----------------------------------------------------------
# Setup
# ----------------------------------------------------------

rm(list = ls())

# !!! Insert your experiment name here as a string, e.g. "MyExperiment" !!!
exp_list <- c("iTPP3_tradeoffs_3rounds", "iTPP3_tradeoffs_4rounds", "iTPP3_tradeoffs_5rounds")

# !!! Insert your predicted parameter here. Note that this must match with one column name in post-processing files !!!
pred_list <- c("inc_red_int_Tot", "sev_red_int_Tot", "mor_red_int_Tot")

library(ggplot2)
library(tidyr)
library(dplyr)
library(scales)
library(patchwork)
library(hetGP)

user <- strsplit(getwd(), "/", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]][5]
GROUP_dr <- "/scicore/home/penny/GROUP/M3TPP/"
setwd(paste0("/scicore/home/penny/", user, "/M3TPP/"))

if (!dir.exists(paste0("./Experiments/", exp_list[2], "/Outputs"))) dir.create(paste0("./Experiments/", exp_list[2], "/Outputs"))

load(paste0(GROUP_dr, exp_list[2], "/param_ranges.RData"))
param_ranges_cont

SPAQ <- as.matrix(c("Coverage1" = 1, "Coverage2" = 0.95, "Halflife" = 31, "Efficacy" = 1))
SPAQ <- (SPAQ - param_ranges_cont[, 1]) / (param_ranges_cont[, 2] - param_ranges_cont[, 1])


# ----------------------------------------------------------
# For each outcome and each experiment, generate SPAQ predictions
# ----------------------------------------------------------

out <- data.frame()

for(exp in exp_list) {
  for (pred in pred_list) {
    
    setting <- basename(Sys.glob(paste0(GROUP_dr, exp, "/gp/trained/", pred, "/*")))
    
    df <- data.frame("setting" = setting) %>%
      separate(col = setting, 
               into = c("Seeds", "Experiment", "Seasonality", "System", "EIR", "Agegroup", "Decay", "Access", "Timing",
                        "Outcome", "Extra1", "Extra2", "Extra3", "Type"),
               sep = "_",
               remove = FALSE)
    
    df <- df[df$Type == "cv.RData", ]
    df <- df[, names(df) %in% c("setting", "Experiment", "Seasonality", "EIR", "Agegroup", "Access", "Outcome")]
    head(df)
    
    df$SPAQ <- NA
    
    for(j in 1:nrow(df)) {
      
      # Load trained emulator
      load(paste0(GROUP_dr, exp, "/gp/trained/", pred, "/", df$setting[j])) # loads object called cv_result
      GP_model <- cv_result$GP_model
      
      # Use emulator to predict likely impact for SPAQ
      SPAQ_pred <- predict(x = t(SPAQ), object = cv_result$GP_model)$mean
      
      df$SPAQ[j] <- SPAQ_pred
      
    }
    
    out <- rbind(out, df)
    
  }
}

remove(cv_result, GP_model, df)


# ----------------------------------------------------------
# Format data
# ----------------------------------------------------------

# Merge in baseline prevalence
prev <- read.csv(paste0("./Experiments/iTPP3_tradeoffs_4rounds/Outputs/Prevalence_prior_to_intervention.csv"))
prev <- prev[, c("Seasonality", "EIR", "Access", "annual_prev")]
prev$annual_prev <- paste0(round(prev$annual_prev*100, 0), "%")
out <- merge(out, prev, by = c("Seasonality", "EIR", "Access"))

# Format labels
out$Seasonality <- ifelse(out$Seasonality == "sharpseasonal", "SHORT SEASON", "LONG SEASON")
out$Access <- ifelse(out$Access == 0.04, "LOW ACCESS", "HIGH ACCESS")
out$Agegroup <- paste0("CHILDREN 3M TO ", out$Agegroup, "Y")
out$Seasonality <- factor(out$Seasonality, levels = c("SHORT SEASON", "LONG SEASON"))
out$Access <- factor(out$Access, levels = c("LOW ACCESS", "HIGH ACCESS"))
out$Outcome <- ifelse(out$Outcome == "inc", "CLINICAL INCIDENCE",
                      ifelse(out$Outcome == "sev", "SEVERE DISEASE",
                             "MORTALITY"))
out$Outcome <- factor(out$Outcome, levels = c("CLINICAL INCIDENCE",
                                              "SEVERE DISEASE",
                                              "MORTALITY"))
out$Experiment <- ifelse(out$Experiment == "iTPP3tradeoffs3rounds", "3 ROUNDS",
                         ifelse(out$Experiment == "iTPP3tradeoffs4rounds", "4 ROUNDS",
                                "5 ROUNDS"))
out$Experiment <- factor(out$Experiment, levels = c("5 ROUNDS", "4 ROUNDS", "3 ROUNDS"))

head(out)


# ----------------------------------------------------------
# Generate plot
# ----------------------------------------------------------


cols <- c("#FAD2CA", "#F2B8AB", "#EA9D8D", "#E2836E", "#D9684F", "#D14E31", "#C93312")
shapes <- c(16, 15, 17)

# Short season

df_plot <- out[out$Seasonality == "SHORT SEASON" & (out$Access == "HIGH ACCESS" & out$Agegroup == "CHILDREN 3M TO 5Y"), ]

p1 <- ggplot(df_plot, aes(x = Outcome, y = SPAQ, colour = annual_prev, shape = Experiment, group = Experiment))

p1 <- p1 + geom_point(position = position_dodge(width = .75), size = 3)

p1 <- p1 + theme(panel.border = element_blank(), 
               panel.background = element_blank(),
               panel.grid.minor = element_blank(),
               panel.grid.major.x = element_blank(),
               panel.grid.major.y = element_line(colour = "lightgrey", linetype = "dotted"),
               text = element_text(family = "Courier New", size = 12),
               axis.line = element_blank(),
               axis.ticks = element_blank(),
               axis.title.x = element_text(face="bold", margin = margin(t = 10)),
               axis.title.y = element_text(face="bold", margin = margin(r = 10)),
               strip.background = element_blank(),
               plot.title = element_text(face = "bold", hjust = 0.5),
               legend.key = element_blank(),
               legend.title = element_text(face = "bold", size = 10),
               legend.text = element_text(size = 10),
              # legend.position = "bottom",
               legend.box = "vertical",
               legend.margin = margin())

p1 <- p1 + scale_colour_manual(values = cols) +
  scale_shape_manual(values = shapes) +
  scale_y_continuous(breaks = seq(10, 90, 10), labels = paste0(seq(10, 90, 10), "%"), limits = c(40, 90))

p1 <- p1 + labs(x = "",
              y = "REDUCTION (%)",
              colour = expression(paste(bold("BASELINE ANNUAL "), bolditalic("Pf"), bold("PR"["2-10"]))),
              shape = "NUMBER OF ROUNDS OF SMC",
              title = "SHORT SEASONAL PROFILE")


# Long season

df_plot <- out[out$Seasonality == "LONG SEASON" & (out$Access == "HIGH ACCESS" & out$Agegroup == "CHILDREN 3M TO 5Y"), ]

p2 <- ggplot(df_plot, aes(x = Outcome, y = SPAQ, colour = annual_prev, shape = Experiment, group = Experiment))

p2 <- p2 + geom_point(position = position_dodge(width = .75), size = 3)

p2 <- p2 + theme(panel.border = element_blank(), 
                 panel.background = element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.grid.major.x = element_blank(),
                 panel.grid.major.y = element_line(colour = "lightgrey", linetype = "dotted"),
                 text = element_text(family = "Courier New", size = 12),
                 axis.line = element_blank(),
                 axis.ticks = element_blank(),
                 axis.title.x = element_text(face="bold", margin = margin(t = 10)),
                 axis.title.y = element_text(face="bold", margin = margin(r = 10)),
                 strip.background = element_blank(),
                 plot.title = element_text(face = "bold", hjust = 0.5),
                 legend.key = element_blank(),
                 legend.title = element_text(face = "bold", size = 10),
                 legend.text = element_text(size = 10),
                # legend.position = "bottom",
                 legend.box = "vertical",
                 legend.margin = margin())

p2 <- p2 + scale_colour_manual(values = cols) +
  scale_shape_manual(values = shapes) +
  scale_y_continuous(breaks = seq(10, 90, 10), labels = paste0(seq(10, 90, 10), "%"), limits = c(40, 90))

p2 <- p2 + labs(x = "",
                y = "REDUCTION (%)",
                colour = expression(paste(bold("BASELINE ANNUAL "), bolditalic("Pf"), bold("PR"["2-10"]))),
                title = "LONG SEASONAL PROFILE") +
  guides(shape = "none")

p1 + p2 + plot_layout(nrow = 2, ncol = 1)

ggsave(filename= paste0("./Experiments/", exp_list[2], "/Outputs/CLINICAL_TRANSLATION_OVERVIEW.jpg"),
       plot = last_plot(),
       width = 9,
       height = 6.5,
       dpi = 400)


