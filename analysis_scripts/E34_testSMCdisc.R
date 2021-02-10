
  
  EIRs = c(3)
  
  

  
    
    EIR <- 3
    
    #import dataset
    
    save_id <- paste0("non_inferiority_",setting[, "seasonality"],"_",   EIR  ,"_",setting[, "SMC_HL"],"_",setting[, "LAI_dec"])
    seeds <- read.table(paste("/scicore/home/penny/GROUP/smc_lai/E3_E4_testSMCdis/postprocessing_5/seeds_E3_E4_testSMCdis_Mali_4.9167_exp_0.1.txt"), 
                        header = T, as.is = TRUE, stringsAsFactors = FALSE)
    mean(seeds$prev_beg)
    
    load("/scicore/home/smith/burlyd00/smc_lai/param_ranges.RData")
    prev <- mean(seeds$prev_beg)
    
    non_inferiority_threshold <- mean(seeds$UL)
    
    #split into test and train
    n_seeds <- length(unique(seeds$seed))
    n_points = round(nrow(seeds)/n_seeds*0.85)
    
    seeds$id <- rep(seq(1,nrow(seeds)/n_seeds),each=n_seeds)
    index_train = sample(nrow(seeds)/n_seeds, n_points)
    index_test = setdiff(1:nrow(seeds)/n_seeds, index_train)
    
    train_data <-seeds[seeds$id %in% index_train,]
    test_data <- seeds[seeds$id %in% index_test,]
    test_data2 =   test_data %>% group_by(Scenario_Name) %>% summarise_at(c(names(test_data)[which(names(test_data)=="Coverage"):length(names(test_data) ) ]),mean,na.rm=TRUE)
    
    # train GP on upper limit of the confidence interval
    predicted= "CI_high_HR"
    
    train_model_matern <- train_GP_matern(train_data[ ,c(rownames(param_ranges), predicted)])
    test_model_matern  <- test_GP(GP_model=train_model_matern, train_data=train_data[,c(rownames(param_ranges), predicted)],
                                  test_data=test_data2[,c(rownames(param_ranges), predicted)])
    
    
    
    test_GP_plot(GP_model=train_model_matern, test=test_data2[,c(rownames(param_ranges), predicted)],save=paste0("/scicore/home/penny/GROUP/smc_lai/E3_E4_testSMCdis/",predicted,".jpeg"))
    
    R2 <- cor(test_model_matern$test_data$CI_high_HR,test_model_matern$test_data$predicted_prev_red)
  
    
    predicted= "UL"
    
    train_model_matern <- train_GP_matern(train_data[ ,c(rownames(param_ranges), predicted)])
    test_model_matern  <- test_GP(GP_model=train_model_matern, train_data=train_data[,c(rownames(param_ranges), predicted)],
                                  test_data=test_data2[,c(rownames(param_ranges), predicted)])
    
    
    
    test_GP_plot(GP_model=train_model_matern, test=test_data2[,c(rownames(param_ranges), predicted)],save=paste0("/scicore/home/penny/GROUP/smc_lai/E3_E4_testSMCdis/",predicted,".jpeg"))
    
    R2 <- cor(test_model_matern$test_data$CI_high_HR,test_model_matern$test_data$predicted_prev_red)
    
      
    #store output
    cv_result <- list(GP_model = train_model_matern)
    cv_result$prev <- prev
    cv_result$predicted <- predicted
    cv_result$R2 <- R2
    cv_result$train_data <- train_data
    cv_result$test_data <- test_data
    cv_result$UL =  mean(seeds$UL)
    
    cv_result$seasonality =  setting[i, "seasonality"]
    cv_result$SMC_HL =  setting[i, "SMC_HL"]
    cv_result$LAI_dec =  setting[i, "LAI_dec"]
    
    save(cv_result, file = paste0("/scicore/home/smith/GROUP/smc_lai/E5_Comparison/gp/",save_id,"_",predicted,".RData"))
    
  }
