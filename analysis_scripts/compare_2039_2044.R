# Script to compare postprocessing results in 5 and 10 years after start of intervention

### SET UP ###

# Load packages
library(dplyr)
library(tidyr)
library(ggplot2)

# Set up
user <- strsplit(getwd(), "/", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]][5]
setwd(paste0("/scicore/home/penny/", user, "/M3TPP"))
GROUP_dr <- GROUP_dr <- "/scicore/home/penny/GROUP/M3TPP/"
exp <- "iTPP3_tradeoffs_4rounds"


### IMPORT DATA ###

# Import settings
setting <- Sys.glob(paste0(GROUP_dr, exp, "/postprocessing/agg*"))
(setting <- sub(".txt", "", basename(setting)))

# Import data
df <- data.frame()

for (i in setting) {
  temp <- read.table(paste0(GROUP_dr, exp, "/postprocessing_2039_yr/", i, ".txt"), header = TRUE, sep = "")
  temp$year <- 2039
  
  df<- rbind(df, temp)
  remove(temp)
  
  temp <- read.table(paste0(GROUP_dr, exp, "/postprocessing/", i, ".txt"), header = TRUE, sep = "")
  temp$year <- 2044
  
  df<- rbind(df, temp)
  remove(temp)
}

df$Scenario_Name <- as.character(df$Scenario_Name)
# df <- df[, !(names(df) %in% c("prev_red_int_May", "prev_red_int_Jun", "prev_red_int_Jul", "prev_red_int_Aug",
#                               "inc_red_int_May", "inc_red_int_Jun", "inc_red_int_Jul", "inc_red_int_Aug",
#                               "counterfactual", "intervention", "counterfactual.1", "intervention.1"))]
str(df)


### VISUALISE DIFFERENCE ###

# Incidence reduction
df_inc <- df %>%
  pivot_wider(id_cols = c(Seasonality, EIR, MaxAge, Access, Coverage1, Coverage2, Halflife, Efficacy),
              names_from = year,
              names_prefix = "year",
              values_from = inc_red_int)

plot_name <- "Incidence reduction - 2039 vs 2044"

p <- ggplot(df_inc, aes(x = year2039, y = year2044)) +
  geom_point(alpha = 0.25) +
  geom_abline(intercept = 0, slope = 1, colour = "red") +
  theme_bw() +
  facet_grid(MaxAge ~ EIR + Access + Seasonality) +
  labs(title = plot_name)

ggsave(filename = paste0("./Experiments/", exp, "/Outputs/", plot_name, ".jpg"),
       plot = p,
       width = 18,
       height = 8,
       dpi = 800)

# Severe disease reduction
df_sev <- df %>%
  pivot_wider(id_cols = c(Seasonality, EIR, MaxAge, Access, Coverage1, Coverage2, Halflife, Efficacy),
              names_from = year,
              names_prefix = "year",
              values_from = sev_red_int)

plot_name <- "Severe disease reduction - 2039 vs 2044"

p <- ggplot(df_sev, aes(x = year2039, y = year2044)) +
  geom_point(alpha = 0.25) +
  geom_abline(intercept = 0, slope = 1, colour = "red") +
  theme_bw() +
  facet_grid(MaxAge ~ EIR + Access + Seasonality) +
  labs(title = plot_name)

ggsave(filename = paste0("./Experiments/", exp, "/Outputs/", plot_name, ".jpg"),
       plot = p,
       width = 18,
       height = 8,
       dpi = 800)

# Mortality reduction
df_mor <- df %>%
  pivot_wider(id_cols = c(Seasonality, EIR, MaxAge, Access, Coverage1, Coverage2, Halflife, Efficacy),
              names_from = year,
              names_prefix = "year",
              values_from = mor_red_int)

plot_name <- "Mortality reduction - 2039 vs 2044"

p <- ggplot(df_mor, aes(x = year2039, y = year2044)) +
  geom_point(alpha = 0.25) +
  geom_abline(intercept = 0, slope = 1, colour = "red") +
  theme_bw() +
  facet_grid(MaxAge ~ EIR + Access + Seasonality)  +
  labs(title = plot_name)

ggsave(filename = paste0("./Experiments/", exp, "/Outputs/", plot_name, ".jpg"),
       plot = p,
       width = 18,
       height = 8,
       dpi = 800)
