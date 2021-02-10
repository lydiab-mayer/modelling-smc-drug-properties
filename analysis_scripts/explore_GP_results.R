library(hetGP)
library(ggplot2)
library(ggpubr)
# load data imported from cluster
load('/scicore/home/smith/GROUP/smc_lai/E2_LAI/gp_5/trained/seeds_E2_LAI_seasonal_Low_indoor_cv.RData')


# get trained model
train_data <- cv_result$train_data
test_data <- cv_result$test_data
GP_model <- cv_result$GP_model

prediction <- predict(x = as.matrix(test_data[,c(2,3,4,5,6,7,8)]), object = GP_model)
plot(prediction$mean, test_data$sum_clin)

GP_model$theta
GP_model$g


#get plot of prediction with EIR on x-axis
x <- matrix(0, nrow = 241, ncol = 7)
# EIR
x[ ,1] <- seq(1,25,0.1)

# Coverage
x[ ,2] <- rep(0.5, nrow(x))

# Access
x[ ,3] <- rep(0.1, nrow(x))

# Halflife
x[ ,4] <- rep(50, nrow(x))

# Efficacy 
x[ ,5] <- rep(0.5, nrow(x))

# Drug_halflife
x[ ,6] <- rep(5, nrow(x))

# Drug_efficacy
x[ ,7] <- rep(0.95, nrow(x))

predicted_EIR <- predict(x = x, object = GP_model)

plot_prediction <- data.frame(predicted_EIR$mean)
colnames(plot_prediction) <- c("mean")
plot_prediction$EIR <- x[,1]
plot_prediction$sd <- predicted_EIR$sd2 + predicted_EIR$nugs


p1 <- ggplot() +
  geom_ribbon(data = plot_prediction, aes(x = EIR, ymin = mean - sqrt(sd), ymax = mean + sqrt(sd), fill = "grey" ), alpha = 0.5) +
  geom_point(data = train_data[train_data$Coverage > 0.4 & train_data$Coverage < 0.6,], aes(x = EIR, y = sum_clin, fill = "red"), alpha = 0.1, colour = "red") +
  geom_line(data = plot_prediction, aes(x = EIR, y = mean, fill = "black")) + 
  scale_fill_manual(values = c("black","grey","red"), labels = c("mean from GP", "mean +/- sd from GP","simulation results"), name = "") + 
  ylab("cumulative clinical cases (10 years)") +
  labs(title = "Relationship between clinical cases and EIR",
       subtitle = "coverage = 0.5, halflife = 50, efficacy = 0.5")+
  theme(plot.title = element_text(size = 12))


#get plot of prediction with Coverage on x-axis
x <- matrix(0, nrow = 101, ncol = 7)

# EIR
x[ ,1] <- rep(12, nrow(x))

# Coverage
x[ ,2] <- seq(0,1,0.01)

# Access
x[ ,3] <- rep(0.1, nrow(x))

# Halflife
x[ ,4] <- rep(50, nrow(x))

# Efficacy 
x[ ,5] <- rep(0.5, nrow(x))

# Drug_halflife
x[ ,6] <- rep(5, nrow(x))

# Drug_efficacy
x[ ,7] <- rep(0.95, nrow(x))

predicted_Coverage <- predict(x = x, object = GP_model)

plot_prediction <- data.frame(predicted_Coverage$mean)
colnames(plot_prediction) <- c("mean")
plot_prediction$Coverage <- x[,2]
plot_prediction$sd <- predicted_Coverage$sd2 + predicted_Coverage$nugs

index <- which(train_data$EIR < 14 & train_data$EIR > 10 )

p2 <- ggplot() +
  geom_ribbon(data = plot_prediction, aes(x = Coverage, ymin = mean - sqrt(sd), ymax = mean + sqrt(sd), fill = "grey" ), alpha = 0.5) +
  geom_point(data = train_data[index, ], aes(x = Coverage, y = sum_clin, fill = "red"), alpha = 0.1, colour = "red") +
  geom_line(data = plot_prediction, aes(x = Coverage, y = mean, fill = "black")) + 
  scale_fill_manual(values = c("black","grey","red"), labels = c("mean from GP", "mean +/- sd from GP","simulation results"), name = "") + 
  ylab("cumulative clinical cases (10 years)") +
  labs(title = "Relationship between clinical cases and coverage",
       subtitle = "EIR = 12, halflife = 50, efficacy = 0.5")+
  theme(plot.title = element_text(size = 12))






#get plot of prediction with Halflife on x-axis
x <- matrix(0, nrow = 161, ncol = 7)

# EIR
x[ ,1] <- rep(12, nrow(x))

# Coverage
x[ ,2] <- rep(0.5, nrow(x))

# Access
x[ ,3] <- rep(0.1, nrow(x))

# Halflife
x[ ,4] <- seq(20, 100, 0.5)

# Efficacy 
x[ ,5] <- rep(0.5, nrow(x))

# Drug_halflife
x[ ,6] <- rep(5, nrow(x))

# Drug_efficacy
x[ ,7] <- rep(0.95, nrow(x))

predicted_Halflife <- predict(x = x, object = GP_model)

plot_prediction <- data.frame(predicted_Halflife$mean)
colnames(plot_prediction) <- c("mean")
plot_prediction$Halflife <- x[,4]
plot_prediction$sd <- predicted_Halflife$sd2 + predicted_Halflife$nugs


p3 <- ggplot() +
  geom_ribbon(data = plot_prediction, aes(x = Halflife, ymin = mean - sqrt(sd), ymax = mean + sqrt(sd), fill = "grey" ), alpha = 0.5) +
  geom_point(data = train_data[index, ], aes(x = Halflife, y = sum_clin, fill = "red"), alpha = 0.1, colour = "red") +
  geom_line(data = plot_prediction, aes(x = Halflife, y = mean, fill = "black")) + 
  scale_fill_manual(values = c("black","grey","red"), labels = c("mean from GP", "mean +/- sd from GP","simulation results"), name = "") + 
  ylab("cumulative clinical cases (10 years)") +
  labs(title = "Relationship between clinical cases and halflife",
       subtitle = "EIR = 12, coverage = 0.5, efficacy = 0.5")+
  theme(plot.title = element_text(size = 12))





#get plot of prediction with Efficacy on x-axis
x <- matrix(0, nrow = 101, ncol = 7)

# EIR
x[ ,1] <- rep(12, nrow(x))

# Coverage
x[ ,2] <- rep(0.5, nrow(x))

# Access
x[ ,3] <- rep(0.1, nrow(x))

# Halflife
x[ ,4] <- rep(50, nrow(x))

# Efficacy 
x[ ,5] <- seq(0, 1, 0.01)

# Drug_halflife
x[ ,6] <- rep(5, nrow(x))

# Drug_efficacy
x[ ,7] <- rep(0.95, nrow(x))

predicted_Efficacy <- predict(x = x, object = GP_model)

plot_prediction <- data.frame(predicted_Efficacy$mean)
colnames(plot_prediction) <- c("mean")
plot_prediction$Efficacy <- x[,5]
plot_prediction$sd <- predicted_Efficacy$sd2 + predicted_Efficacy$nugs


p4 <- ggplot() +
  geom_ribbon(data = plot_prediction, aes(x = Efficacy, ymin = mean - sqrt(sd), ymax = mean + sqrt(sd), fill = "grey" ), alpha = 0.5) +
  geom_point(data = train_data[index,], aes(x = Efficacy, y = sum_clin, fill = "red"), alpha = 0.1, colour = "red") +
  geom_line(data = plot_prediction, aes(x = Efficacy, y = mean, fill = "black")) + 
  scale_fill_manual(values = c("black","grey","red"), labels = c("mean from GP", "mean +/- sd from GP","simulation results"), name = "") + 
  ylab("cumulative clinical cases (10 years)") +
  labs(title = "Relationship between clinical cases and efficacy",
       subtitle = "EIR = 12, halflife = 50, coverage = 0.5")+
  theme(plot.title = element_text(size = 12))




p <- ggarrange(p1, p2, p3, p4, ncol = 2, nrow = 2)
p <- annotate_figure(p, top = text_grob("Traning a GP (without adaptive sampling) using the LAI simulation outputs \nand projecting on one variable while fixing the other three", color = "black", face = "bold", size = 14))

p

