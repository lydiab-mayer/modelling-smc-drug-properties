setwd("~/smc_lai/analysis_workflow/analysis_scripts")

library(hetGP)
library(ggplot2)
library(ggpubr)
library(reshape2)
library(tgp)



plot_list = list()
EIR <- c(5,10,25)
SMC_halflife <- c(10,20,31.9)
plot_index <- 1
for (EIR_index in 1:length(EIR)){
  
  for (SMC_halflife_index in 1:length(SMC_halflife)){
    print(paste("generate plot", plot_index, "of ", length(EIR)*length(SMC_halflife)))
    #import dataset
    seeds <- read.table(paste('/scicore/home/smith/GROUP/smc_lai/E5_Comparison/postprocessing/seeds_E5_comparison_seasonal_Low_indoor_',
                              EIR[EIR_index],'_',SMC_halflife[SMC_halflife_index],'.txt', sep = ""), 
                        header = T, as.is = TRUE, stringsAsFactors = FALSE)
    
    non_inferiority_threshold <- mean(seeds$UL)
    
    #split into test and train
    n_seeds <- length(unique(seeds$seed))
    n_points = round(nrow(seeds)/n_seeds*0.75)
    
    seeds$id <- rep(seq(1,nrow(seeds)/n_seeds),each=n_seeds)
    index_train = sample(nrow(seeds)/n_seeds, n_points)
    index_test = setdiff(1:nrow(seeds)/n_seeds, index_train)
    
    train_data <-seeds[seeds$id %in% index_train,]
    test_data <- seeds[seeds$id %in% index_test,]
    
    # train GP on upper limit of the confidence interval
    param_col <- c("Halflife", "Efficacy")
    response_col <- "CI_high_HR"
    
    train_model_matern <- train_GP_matern(train_data[ ,c(rownames(param_ranges), predicted)])
    test_model_matern  <- test_GP(GP_model=train_model_matern, train_data=train_data[,c(rownames(param_ranges), predicted)],
                                  test_data=test_data2[,c(rownames(param_ranges), predicted)])
    test_GP_plot(GP_model=train_model_matern, test=test_data2[,c(rownames(param_ranges), predicted)])
    
    R2 <- cor(test_model_matern$test_data$CI_high_HR,test_model_matern$test_data$predicted_prev_red)
  
        Halflife =c(30,150)
    Efficacy= c(0,1)
    param_ranges = rbind(Halflife,Efficacy)
    
    Xcand = data.frame(lhs(3000, param_ranges))
    
         prediction_for_plot <- predict(x = as.matrix(Xcand), 
                                       object = GP_model)
         upper_limit_CI   <- data.frame(cbind(Xcand,pred=prediction_for_plot$mean ))
      
    
    
    colnames(upper_limit_CI) <- c("halflife", "efficacy", "upper_bound")
    df <- upper_limit_CI
    df$non_inferiority <-   ifelse(df$upper_bound>non_inferiority_threshold,0, 1)
    
    
    p <- ggplot(data = df, aes(x= efficacy, y= halflife, colour=as.factor(non_inferiority))) + 
      geom_point() +
      scale_color_manual(values=c("red","green"), labels = c("inferior", "non inferior")) +
      labs(title = paste("EIR = ", EIR[EIR_index], ", SMC halflife = ", SMC_halflife[SMC_halflife_index]), colour = "" ) +
      theme(title = element_text(size=14),
            axis.title = element_text(size = 14),
            axis.text = element_text(size = 12),
            legend.text = element_text(size = 12))

    
     plot_list[[plot_index]] = p
    plot_index <- plot_index + 1
  }
}


p <- ggarrange(plotlist = plot_list)
ggsave('outfiles/plot_non_inferiority_FroschimMixer.pdf', p, width = 15, height = 15)


