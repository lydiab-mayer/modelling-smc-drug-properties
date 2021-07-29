###############################################################
###
### Generating generic seasonality from normal distributions
### Used for MOCK example
### NN
### 15.06.2021
###
###############################################################


#########

rm(list = ls())
setwd("./M3TPP")
set.seed(42)

options(scipen=999)

#########

round_df <- function(x, digits) {
  numeric_columns <- sapply(x, mode) == 'numeric'
  x[numeric_columns] <-  round(x[numeric_columns], digits)
  x
}

library(ggplot2)

####### Generate seasonal distributions

plot(function(x) dnorm(x, mean = 6, sd = 1.5), 1, 12)
curve(log(dnorm(x)), add = TRUE, col = "red", lwd = 2)

dnorm(seq(1:12), mean = 6, sd = 0.2) # 1 month
dnorm(seq(1:12), mean = 6, sd = 0.5) # 2 months
dnorm(seq(1:12), mean = 6, sd = 0.75) # 3 months
dnorm(seq(1:12), mean = 6, sd = 1) # 4 months
dnorm(seq(1:12), mean = 6, sd = 1.5) # 5 months

####### Save data

Seasons = as.data.frame(rbind(seas1mo = dnorm(seq(1:12), mean = 6, sd = 0.45),
                        seas2mo = dnorm(seq(1:12), mean = 6, sd = 0.6),
                        seas3mo = dnorm(seq(1:12), mean = 6, sd = 0.75),
                        seas4mo = dnorm(seq(1:12), mean = 6, sd = 1),
                        seas5mo = dnorm(seq(1:12), mean = 6, sd = 1.5)))

Seasons$Seasonality = rownames(Seasons)

colnames(Seasons) = c("Jan","Feb","Mar","Apr","May","Jun",
                      "Jul","Aug","Sep","Oct","Nov","Dec","Seasonality")
Seasons = Seasons[,c("Seasonality","Jan","Feb","Mar","Apr","May","Jun",
                     "Jul","Aug","Sep","Oct","Nov","Dec")]

# Round numbers
Seasons = round_df(Seasons, 8)

# Plot data
Seasonality_plot_data=Seasons[,2:ncol(Seasons)]
Seasonality_plot_data=t(Seasonality_plot_data)
colnames(Seasonality_plot_data)=Seasons[,1]
Seasonality_plot_data=as.data.frame(Seasonality_plot_data)
Seasonality_plot_data$month = 1:12
Seasonality_plot_data_clean=as.data.frame(rbind(
  cbind(month=Seasonality_plot_data[,6], type=rep(colnames(Seasonality_plot_data)[1],12), value=Seasonality_plot_data[,1]),
  cbind(month=Seasonality_plot_data[,6], type=rep(colnames(Seasonality_plot_data)[2],12), value=Seasonality_plot_data[,2]),
  cbind(month=Seasonality_plot_data[,6], type=rep(colnames(Seasonality_plot_data)[3],12), value=Seasonality_plot_data[,3]),
  cbind(month=Seasonality_plot_data[,6], type=rep(colnames(Seasonality_plot_data)[4],12), value=Seasonality_plot_data[,4]),
  cbind(month=Seasonality_plot_data[,6], type=rep(colnames(Seasonality_plot_data)[5],12), value=Seasonality_plot_data[,5])))
Seasonality_plot_data_clean$month = as.numeric(as.character(Seasonality_plot_data_clean$month))
# Seasonality_plot_data_clean$month = factor(Seasonality_plot_data_clean$month, levels = as.character(c(seq(6,12,1),seq(1,5,1))))
Seasonality_plot_data_clean$value = as.numeric(as.character(Seasonality_plot_data_clean$value))

# Plot
ggplot(Seasonality_plot_data_clean) +
  geom_line(aes(x=month, y=value, group = type, color = type)) +
  scale_x_continuous(breaks = seq(1,12,1))

# Save data
Seasons$Seasonality = factor(Seasons$Seasonality)
# write.table(Seasons, file = "/scicore/home/penny/nekkab0000/M3TPP/Experiments/exp_5/seasonality.txt", 
#             sep="\t", quote = F, row.names = F)




#########################################
###
### Mock example seasonalities
### EXP_2MAB_5
### NN 26/04/2021
###
#########################################

# Step 1: run sims for flat, 9 mo, 6 mo, 3 mo, 2 mo, 1 mo, seasonal profiles
# Replot

# Load seasonality profiles and plot
Seasonality =read.table(paste0("./Experiments/seasonality.txt"), sep="\t", header = TRUE)

# Reshape
Seasonality_plot_data=Seasonality[,2:ncol(Seasonality)]
Seasonality_plot_data=t(Seasonality_plot_data)
colnames(Seasonality_plot_data)=Seasonality[,1]
Seasonality_plot_data=as.data.frame(Seasonality_plot_data)
Seasonality_plot_data$month = 1:12
Seasonality_plot_data = Seasonality_plot_data %>% 
  group_by(month) %>% 
  gather(key = season, value = value, -month) %>% 
  ungroup()

# Plot
library(ggplot2)
ggplot(Seasonality_plot_data) +
  geom_line(aes(x=month, y=value, group = season, color = season)) +
  scale_x_continuous(breaks = seq(1,12,1))


############## Choosing length
# https://www.mdpi.com/2306-5729/5/2/31/htm

# Choose 3, 4, 5, 6, 9
dnorm(seq(1:12), mean = 6, sd = 0.5) # 3 months
dnorm(seq(1:12), mean = 6, sd = 1) # 4 months
dnorm(seq(1:12), mean = 6, sd = 1.5) # 5 months
dnorm(seq(1:12), mean = 6, sd = 1.7) # 6 months
dnorm(seq(1:12), mean = 6, sd = 3) # 9 months

# Data
Seasons = as.data.frame(rbind(seas3mo = dnorm(seq(1:12), mean = 6, sd = 0.25),
                              seas4mo = dnorm(seq(1:12), mean = 6, sd = 1.1),
                              seas5mo = dnorm(seq(1:12), mean = 6, sd = 1.5),
                              seas6mo = dnorm(seq(1:12), mean = 6, sd = 1.7),
                              seas9mo = dnorm(seq(1:12), mean = 6, sd = 4)))

Seasons$Seasonality = rownames(Seasons)

colnames(Seasons) = c("Jan","Feb","Mar","Apr","May","Jun",
                      "Jul","Aug","Sep","Oct","Nov","Dec","Seasonality")
Seasons = Seasons[,c("Seasonality","Jan","Feb","Mar","Apr","May","Jun",
                     "Jul","Aug","Sep","Oct","Nov","Dec")]

# Round numbers
Seasons = round_df(Seasons, 8)

# Plot data
Seasonality_plot_data=Seasons[,2:ncol(Seasons)]
Seasonality_plot_data=t(Seasonality_plot_data)
colnames(Seasonality_plot_data)=Seasons[,1]
Seasonality_plot_data=as.data.frame(Seasonality_plot_data)
Seasonality_plot_data$month = 1:12
Seasonality_plot_data_clean=as.data.frame(rbind(
  cbind(month=Seasonality_plot_data[,6], type=rep(colnames(Seasonality_plot_data)[1],12), value=Seasonality_plot_data[,1]),
  cbind(month=Seasonality_plot_data[,6], type=rep(colnames(Seasonality_plot_data)[2],12), value=Seasonality_plot_data[,2]),
  cbind(month=Seasonality_plot_data[,6], type=rep(colnames(Seasonality_plot_data)[3],12), value=Seasonality_plot_data[,3]),
  cbind(month=Seasonality_plot_data[,6], type=rep(colnames(Seasonality_plot_data)[4],12), value=Seasonality_plot_data[,4]),
  cbind(month=Seasonality_plot_data[,6], type=rep(colnames(Seasonality_plot_data)[5],12), value=Seasonality_plot_data[,5])))
Seasonality_plot_data_clean$month = as.numeric(as.character(Seasonality_plot_data_clean$month))
# Seasonality_plot_data_clean$month = factor(Seasonality_plot_data_clean$month, levels = as.character(c(seq(6,12,1),seq(1,5,1))))
Seasonality_plot_data_clean$value = as.numeric(as.character(Seasonality_plot_data_clean$value))

# Plot
ggplot(Seasonality_plot_data_clean) +
  geom_line(aes(x=month, y=value, group = type, color = type)) +
  scale_x_continuous(breaks = seq(1,12,1))

# Save data
Seasons$Seasonality = factor(Seasons$Seasonality)
# write.table(Seasons, file = "/scicore/home/penny/nekkab0000/M3TPP/Seasonality_analysis/seasonality_34569_months_updatedMay302021.txt", 
#             sep="\t", quote = F, row.names = F)


################## remove 3 months

Seasons = Seasons[which(Seasons$Seasonality != "seas3mo"),]
Seasons$Seasonality = factor(Seasons$Seasonality)

# Save data
# write.table(Seasons, file = "/scicore/home/penny/nekkab0000/M3TPP/Seasonality_analysis/seasonality_4569_months_updatedMay302021.txt", 
#             sep="\t", quote = F, row.names = F)


################## only 4 months

Seasons = Seasons[which(Seasons$Seasonality == "seas4mo"),]
Seasons$Seasonality = factor(Seasons$Seasonality)

# Save data
# write.table(Seasons, file = "/scicore/home/penny/nekkab0000/M3TPP/Seasonality_analysis/seasonality_4_months_updatedMay302021.txt", 
#             sep="\t", quote = F, row.names = F)


