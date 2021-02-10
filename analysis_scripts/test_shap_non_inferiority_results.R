setwd("~/smc_lai/analysis_workflow/analysis_scripts")
source('~/smc_lai/analysis_workflow/3_GP_train/GP_toolbox.R')


library(hetGP)
library(ggplot2)
library(ggpubr)
library(reshape2)
library(tgp)
library(xgboost)
library(pROC)
library(caret)
library(SHAPforxgboost)


plot_list = list()
EIR <- c(5,10,25)
SMC_halflife <- c(10,20,31.9)
plot_index <- 1
#for (EIR_index in 1:length(EIR)){
  
 # for (SMC_halflife_index in 1:length(SMC_halflife)){
   
    EIR_index=1
    SMC_halflife_index=1
    
     print(paste("generate plot", plot_index, "of ", length(EIR)*length(SMC_halflife)))
    #import dataset
  #  seeds <- read.table(paste('/scicore/home/smith/GROUP/smc_lai/E5_Comparison/postprocessing/seeds_E5_comparison_seasonal_Low_indoor_',
   #                           EIR[EIR_index],'_',SMC_halflife[SMC_halflife_index],'.txt', sep = ""), 
    #                    header = T, as.is = TRUE, stringsAsFactors = FALSE)
    
     split_file = "/scicore/home/smith/GROUP/smc_lai/Test_SMC/postprocessing_5/seeds_E3_E4_comparison_seasonal_Low_indoor_exp.txt"
     
     seeds <- read.table(split_file,   header = T, as.is = TRUE, stringsAsFactors = FALSE)
                               
     ranges_file = "/scicore/home/smith/GROUP/smc_lai/Test_SMC/param_ranges.RData"
     
    load(ranges_file)


    #split into test and train
    n_seeds <- length(unique(seeds$seed))
    n_points = round(nrow(seeds)/n_seeds*0.80)
    
    seeds$id <- rep(seq(1,nrow(seeds)/n_seeds),each=n_seeds)
    index_train = sample(nrow(seeds)/n_seeds, n_points)
    index_test = setdiff(1:nrow(seeds)/n_seeds, index_train)
    
    train <-seeds[seeds$id %in% index_train,]
    test <- seeds[seeds$id %in% index_test,]
    test_data2 =   test_data %>% group_by(Scenario_Name) %>% summarise_at(c(names(test_data)[which(names(test_data)=="Halflife"):length(names(test_data) ) ]),mean,na.rm=TRUE)
    
    # train GP on upper limit of the confidence interval
    predicted= "non_inferiority"
    
    
    xgb.train.data = xgb.DMatrix(data.matrix(train[,rownames(param_ranges)]), label = train[,predicted], missing = NA)
    
    param <- list(objective = "binary:logistic", base_score = 0.5)
    
    cv <- createFolds(train[,"id"],k=10)
    
    xgboost.cv = xgb.cv(param=param, data = xgb.train.data, folds =cv , nrounds = 4000, early_stopping_rounds = 100, metrics='auc')
   
    
     best_iteration = xgboost.cv$best_iteration
     xgb.model <- xgboost(param =param,  data = xgb.train.data, nrounds=best_iteration)
     xgb.test.data = xgb.DMatrix(data.matrix(test[,rownames(param_ranges)]), missing = NA)
    xgb.preds = predict(xgb.model, xgb.test.data)
    xgb.roc_obj <- roc(test[,predicted], xgb.preds)
    
    plot(xgb.roc_obj )
    
    
    cat("XGB AUC ", auc(xgb.roc_obj))
    
    
    
    dataX  <- test[,rownames(param_ranges)]
    shap_values <- shap.values(xgb_model = xgb.model, X_train = dataX)
    # The ranked features by mean |SHAP|
    shap_values$mean_shap_score
      
      
    shap_long <- shap.prep(xgb_model = xgb.model, X_train = dataX)
    
    shap.plot.summary(shap_long)
    
    shap.plot.summary(shap_long, x_bound  = 1.2, dilute = 10)
    
    xgboost::xgb.plot.shap(data = data.matrix(dataX), model = xgb.model, top_n = 4, n_col = 2)
   
    
    
    
    
    
    
    
    
    
    
    
    
    
    
     
    calc_shapley_effects = function(GP_model,param_names, param_spec, num_points){
      S_mat = T_mat = NULL
      print(param_spec)
      # define wrapper for the GP_model prediction to be called by shapley function
      GP_f = function(X){

        xgb.test.data = xgb.DMatrix(X, missing = NA)
        pred = predict(GP_model, xgb.test.data)
       out <- ifelse(pred < 0.5,0,1)
        
        return(out)
      }
      
      d <- nrow(param_spec)
      
      # Xall: a function to generate a n-sample of a d-dimensional input vector 
      # (following the required joint distribution). Here: joint distribution = uniform?
      Xall <- function(n) {
        a <- runif(n, param_spec[1,1], param_spec[1,2])
        b <- runif(n, param_spec[2,1], param_spec[2,2])
        
        mat <- matrix(cbind(a,b), nrow = n)
        
        colnames(mat) <- c("Halflife","Efficacy")
        
        return(mat)
      }
      
      # Xset: Xset(n, Sj, Sjc, xjc) is a function to generate a n-sample of a d-dimensional 
      # input vector corresponding to the indices in Sj conditional on the input values xjc with 
      # the index set Sjc (following the required joint distribution). Here: joint distribution = uniform?
      
      Xset <- function(n, Sj, Sjc, xjc){
        a <- runif(n, param_spec[1,1], param_spec[1,2])
        b <- runif(n, param_spec[2,1], param_spec[2,2])
        m <- matrix(cbind(a,b), nrow = n)

        colnames(m) <- c("Halflife","Efficacy")
        
        return(m[,Sj])
      }
      
      shapley <- shapleyPermEx(model = GP_f, Xall, Xset,d=d, Nv=1e4, No = 1e3, Ni = 2)
      
      return(shapley)
    }
    
    
    
    # shapley indices for non-inferiority
    param_ranges <- matrix(nrow = 2, ncol = 2)
    colnames(param_ranges) <- c("Halflife","Efficacy")
    param_ranges[,"Halflife" ] <- c(50,130)
    param_ranges[, "Efficacy"] <- c(0,1)
    
    
    a <- calc_shapley_effects(xgb.model,rownames(param_ranges), param_ranges)
    
    
    
    
    
    
    
    
      }
}


p <- ggarrange(plotlist = plot_list)
ggsave('./Outputs/E5_plot_non_inferiority.pdf', p, width = 15, height = 15)


