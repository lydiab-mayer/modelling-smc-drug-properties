library(akima)
library(ggpubr)
seeds <- 5

smc_lai_df <- as.data.frame(read.table('/scicore/home/smith/GROUP/smc_lai/E4_SMCLAI/postprocessing_5/seeds_E4_SMCLAI_seasonal_Low_indoor.txt', 
                                       header = T))
smc_smc_df <- as.data.frame(read.table('/scicore/home/smith/GROUP/smc_lai/E3_SMCSMC/postprocessing_5/seeds_E3_SMCSMC_seasonal_Low_indoor.txt',
                                       header = T))


# remove simulations with no difference in cases between the two scencarios (these are similations
# where the health system alone is strong enough to reach elimination regardless of interventions)
index <- which(smc_lai_df$sum_clin - smc_smc_df$sum_clin == 0)
smc_lai_df <- smc_lai_df[- index, ] 
smc_smc_df <- smc_smc_df[- index, ]


# get mean and absulute difference between max and min outcome for each point (5 sims per point)

smc_lai_summary_df <- data.frame(matrix(0, nrow = nrow(smc_lai_df)/seeds, ncol = 0))
smc_lai_summary_df$EIR <- 0
smc_lai_summary_df$Coverage <- 0
smc_lai_summary_df$Halflife <- 0
smc_lai_summary_df$Efficacy <- 0
smc_lai_summary_df$mean_sum_clin_lai <- 0
smc_lai_summary_df$abs_max_min_sum_clin_lai <- 0
smc_lai_summary_df$mean_sum_clin_smc <- 0
smc_lai_summary_df$abs_max_min_sum_clin_smc <- 0

for (k in 1:(nrow(smc_lai_df)/seeds)){
  first_row <- 5*(k-1) + 1
  last_row <- 5*k
  
  smc_lai_summary_df[k, ]$EIR <- smc_lai_df[first_row, ]$EIR
  smc_lai_summary_df[k, ]$Coverage <- smc_lai_df[first_row, ]$Coverage
  smc_lai_summary_df[k, ]$Halflife <- smc_lai_df[first_row, ]$Halflife
  smc_lai_summary_df[k, ]$Efficacy <- smc_lai_df[first_row, ]$Efficacy
  
  smc_lai_summary_df[k, ]$mean_sum_clin_lai <- sum(smc_lai_df[first_row:last_row, ]$sum_clin)/seeds
  smc_lai_summary_df[k, ]$abs_max_min_sum_clin_lai <- abs(max(smc_lai_df[first_row:last_row, ]$sum_clin) - min(smc_lai_df[first_row:last_row, ]$sum_clin))
  smc_lai_summary_df[k, ]$mean_sum_clin_smc <- sum(smc_smc_df[first_row:last_row, ]$sum_clin)/seeds
  smc_lai_summary_df[k, ]$abs_max_min_sum_clin_smc <- abs(max(smc_smc_df[first_row:last_row, ]$sum_clin) - min(smc_smc_df[first_row:last_row, ]$sum_clin))
}





# max difference in mean clinical cases: 
upper_lim_col <- max(smc_lai_summary_df$mean_sum_clin_lai -smc_lai_summary_df$mean_sum_clin_smc )
lower_lim_col <- min(smc_lai_summary_df$mean_sum_clin_lai -smc_lai_summary_df$mean_sum_clin_smc )

# get interpolation and plot for EIR 0 - 12 & Coverage 0.5 - 1
index <- which(smc_lai_summary_df$EIR <= 12 & smc_lai_summary_df$Coverage > 0.5)

x <- smc_lai_summary_df[index, ]$Halflife
y <- smc_lai_summary_df[index, ]$Efficacy
z <- smc_lai_summary_df[index, ]$mean_sum_clin_lai - smc_lai_summary_df[index, ]$mean_sum_clin_smc

interpolation_mean <- interp(x, y, z, xo=seq(min(x), max(x), length = 20),
                             yo=seq(min(y), max(y), length = 20),
                             linear = TRUE, extrap=FALSE, duplicate = "error", dupfun = NULL,
                             nx = 10, ny = 10,
                             jitter = 10^-12, jitter.iter = 6, jitter.random = FALSE)




interpolation_mean_dataframe <- data.frame(expand.grid(x = interpolation_mean$x, y = interpolation_mean$y))
interpolation_mean_dataframe$z <- array(interpolation_mean$z)
colnames(interpolation_mean_dataframe) <- c("Halflife", "Efficacy", "Diff_Clin_Cases")

p1 <- ggplot(data = interpolation_mean_dataframe) +
  geom_raster(aes(x = Halflife, y = Efficacy, fill = Diff_Clin_Cases), alpha = 0.8)+
  geom_point(data = smc_lai_summary_df[index, ], aes(x = Halflife, y = Efficacy, col = mean_sum_clin_lai - mean_sum_clin_smc)) +
  scale_color_gradient(low="firebrick", high="springgreen4", limits = c(lower_lim_col, upper_lim_col)) + 
  scale_fill_gradient(low="firebrick", high="springgreen4", limits = c(lower_lim_col, upper_lim_col)) +
  guides( color = FALSE) +
  labs(title = "low EIR (0 - 12), high coverage (0.5 -1)", fill = "difference in \nclinical cases")



# get interpolation and plot for EIR 12 - 25 & Coverage 0.5 - 1
index <- which(smc_lai_summary_df$EIR > 12 & smc_lai_summary_df$Coverage > 0.5)

x <- smc_lai_summary_df[index, ]$Halflife
y <- smc_lai_summary_df[index, ]$Efficacy
z <- smc_lai_summary_df[index, ]$mean_sum_clin_lai - smc_lai_summary_df[index, ]$mean_sum_clin_smc

interpolation_mean <- interp(x, y, z, xo=seq(min(x), max(x), length = 20),
                             yo=seq(min(y), max(y), length = 20),
                             linear = TRUE, extrap=FALSE, duplicate = "error", dupfun = NULL,
                             nx = 10, ny = 10,
                             jitter = 10^-12, jitter.iter = 6, jitter.random = FALSE)




interpolation_mean_dataframe <- data.frame(expand.grid(x = interpolation_mean$x, y = interpolation_mean$y))
interpolation_mean_dataframe$z <- array(interpolation_mean$z)
colnames(interpolation_mean_dataframe) <- c("Halflife", "Efficacy", "Diff_Clin_Cases")

p2 <- ggplot(data = interpolation_mean_dataframe) +
  geom_raster(aes(x = Halflife, y = Efficacy, fill = Diff_Clin_Cases), alpha = 0.8)+
  geom_point(data = smc_lai_summary_df[index, ], aes(x = Halflife, y = Efficacy, col = mean_sum_clin_lai - mean_sum_clin_smc)) +
  scale_color_gradient(low="firebrick", high="springgreen4", limits = c(lower_lim_col, upper_lim_col)) + 
  scale_fill_gradient(low="firebrick", high="springgreen4", limits = c(lower_lim_col, upper_lim_col)) +
  guides( color = FALSE) +
  labs(title = "high EIR (12 - 25), high coverage (0.5 -1)", fill = "difference in \nclinical cases")



# get interpolation and plot for EIR 0 - 12 & Coverage 0 - 0.5 - 1
index <- which(smc_lai_summary_df$EIR <= 12 & smc_lai_summary_df$Coverage <= 0.5)

x <- smc_lai_summary_df[index, ]$Halflife
y <- smc_lai_summary_df[index, ]$Efficacy
z <- smc_lai_summary_df[index, ]$mean_sum_clin_lai - smc_lai_summary_df[index, ]$mean_sum_clin_smc

interpolation_mean <- interp(x, y, z, xo=seq(min(x), max(x), length = 20),
                             yo=seq(min(y), max(y), length = 20),
                             linear = TRUE, extrap=FALSE, duplicate = "error", dupfun = NULL,
                             nx = 10, ny = 10,
                             jitter = 10^-12, jitter.iter = 6, jitter.random = FALSE)




interpolation_mean_dataframe <- data.frame(expand.grid(x = interpolation_mean$x, y = interpolation_mean$y))
interpolation_mean_dataframe$z <- array(interpolation_mean$z)
colnames(interpolation_mean_dataframe) <- c("Halflife", "Efficacy", "Diff_Clin_Cases")

p3 <- ggplot(data = interpolation_mean_dataframe) +
  geom_raster(aes(x = Halflife, y = Efficacy, fill = Diff_Clin_Cases), alpha = 0.8)+
  geom_point(data = smc_lai_summary_df[index, ], aes(x = Halflife, y = Efficacy, col = mean_sum_clin_lai - mean_sum_clin_smc)) +
  scale_color_gradient(low="firebrick", high="springgreen4", limits = c(lower_lim_col, upper_lim_col)) + 
  scale_fill_gradient(low="firebrick", high="springgreen4", limits = c(lower_lim_col, upper_lim_col)) +
  guides( color = FALSE) +
  labs(title = "low EIR (0 - 12), low coverage (0 - 0.5)", fill = "difference in \nclinical cases")




# get interpolation and plot for EIR 12 - 20 & Coverage 0 - 0.5
index <- which(smc_lai_summary_df$EIR > 12 & smc_lai_summary_df$Coverage <= 0.5)

x <- smc_lai_summary_df[index, ]$Halflife
y <- smc_lai_summary_df[index, ]$Efficacy
z <- smc_lai_summary_df[index, ]$mean_sum_clin_lai - smc_lai_summary_df[index, ]$mean_sum_clin_smc

interpolation_mean <- interp(x, y, z, xo=seq(min(x), max(x), length = 20),
                             yo=seq(min(y), max(y), length = 20),
                             linear = TRUE, extrap=FALSE, duplicate = "error", dupfun = NULL,
                             nx = 10, ny = 10,
                             jitter = 10^-12, jitter.iter = 6, jitter.random = FALSE)




interpolation_mean_dataframe <- data.frame(expand.grid(x = interpolation_mean$x, y = interpolation_mean$y))
interpolation_mean_dataframe$z <- array(interpolation_mean$z)
colnames(interpolation_mean_dataframe) <- c("Halflife", "Efficacy", "Diff_Clin_Cases")

p4 <- ggplot(data = interpolation_mean_dataframe) +
  geom_raster(aes(x = Halflife, y = Efficacy, fill = Diff_Clin_Cases), alpha = 0.8)+
  geom_point(data = smc_lai_summary_df[index, ], aes(x = Halflife, y = Efficacy, col = mean_sum_clin_lai - mean_sum_clin_smc)) +
  scale_color_gradient(low="firebrick", high="springgreen4", limits = c(lower_lim_col, upper_lim_col)) + 
  scale_fill_gradient(low="firebrick", high="springgreen4", limits = c(lower_lim_col, upper_lim_col)) +
  guides( color = FALSE) +
  labs(title = "high EIR (12 - 25), low coverage (0 - 0.5)", fill = "difference in \nclinical cases")


p <- ggarrange(p1, p2, p3, p4, ncol = 2, nrow = 2)
annotate_figure(p, top = text_grob("Absoulute difference in clinical cases after 10 years (SCM/LAI - SMC only)", color = "black", face = "bold", size = 14))

