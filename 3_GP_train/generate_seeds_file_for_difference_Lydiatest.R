source('~/smc_lai/analysis_workflow/3_GP_train/GP_toolbox.R')

load("~/smc_lai/param_ranges.RData")
# generate a table to run the GP on the difference between two columns in seeds tables from different experiments

seeds_smc_smc <- read.table("/scicore/home/smith/GROUP/smc_lai/E3_SMCSMC/postprocessing_5/seeds_E3_SMCSMC_seasonal_Low_indoor.txt", header = T, as.is = TRUE, stringsAsFactors = FALSE)
seeds_smc_lai <- read.table("/scicore/home/smith/GROUP/smc_lai/E4_SMCLAI/postprocessing_5/seeds_E4_SMCLAI_seasonal_Low_indoor.txt", header = T, as.is = TRUE, stringsAsFactors = FALSE)


seeds <- length(unique(seeds_smc_lai$seed))
  
#remove runs where prevalence in year 5 is 0 (due to eliminating before switch)
scendellai <-seeds_smc_lai[which( seeds_smc_lai$iprev_y5_1 == 0),"Scenario_Name"]
scendelscm <- seeds_smc_smc[which( seeds_smc_smc$iprev_y5_1 == 0),"Scenario_Name"]
index <- c(scendellai,scendelscm)

seeds_comparison <- data.frame(seeds_smc_lai[- which(seeds_smc_lai$Scenario_Name %in% index), ])
seeds_comparison$abs_diff_pppy_y10_all <- seeds_smc_lai[- which(seeds_smc_lai$Scenario_Name %in% index),]$pppy_y10_all - seeds_smc_smc[- which(seeds_smc_smc$Scenario_Name %in% index),]$pppy_y10_all
seeds_comparison$rel_diff_pppy_y10_all <- (seeds_smc_lai[- which(seeds_smc_lai$Scenario_Name %in% index),]$pppy_y10_all - seeds_smc_smc[- which(seeds_smc_smc$Scenario_Name %in% index),]$pppy_y10_all)/
  seeds_smc_smc[- which(seeds_smc_smc$Scenario_Name %in% index),]$pppy_y10_all

# ml: first attempt to assign simulations to 0/1 based on non-inferiority measure
# for now use 1.64 as a cutoff (from Zongo et. al.). come back to this later.
seeds_comparison$odds_ratio <- (seeds_smc_lai[index,]$pppy_y10_all/(1-seeds_smc_lai[index,]$pppy_y10_all))/(seeds_smc_smc[index,]$pppy_y10_all/(1-seeds_smc_smc[index,]$pppy_y10_all))

seeds_comparison$non_inferiority <- 0
#seeds_comparison[seeds_comparison$odds_ratio > 1.64,]$non_inferiority <- 1

predicted= "rel_diff_pppy_y10_all"

#write.table(seeds_comparison, file = "/scicore/home/smith/GROUP/smc_lai/E3_E4_comparison/postprocessing_5/seeds_E3_E4_comparison_seasonal_Low_indoor_testLy.txt",  sep = "\t", col.names = T, row.names = F)


seeds_comparison$id <- rep(seq(1,nrow(seeds_comparison)/seeds),each=seeds)


num_points = round(nrow(seeds_comparison)/seeds*0.75)

index_train = sample(nrow(seeds_comparison)/seeds, num_points)
index_test = setdiff(1:nrow(seeds_comparison)/seeds, index_train)


train_data <-seeds_comparison[seeds_comparison$id %in% index_train,]
test_data <- seeds_comparison[seeds_comparison$id %in% index_test,]
test_data2 =   test_data %>% group_by(Scenario_Name) %>% summarise_at(c(names(test_data)[which(names(processed_OM_sim)=="seed"):length(names(test_data) ) ]),mean,na.rm=TRUE)


# Parameter columns


trained_model = train_GP(train_data[ , c(rownames(param_ranges), predicted)] ) 
test_model<- test_GP(GP_model=trained_model, train_data=train_data[,c(rownames(param_ranges), predicted)], test_data=test_data2[,c(rownames(param_ranges), predicted)])
test_GP_plot(GP_model=trained_model, test=test_data2[,c(rownames(param_ranges), predicted)]) 


train_model_matern <- train_GP_matern(train_data[ ,c(rownames(param_ranges), predicted)]) 
test_model_matern  <- test_GP(GP_model=train_model_matern, train_data=train_data[,c(rownames(param_ranges), predicted)],
                            test_data=test_data2[,c(rownames(param_ranges), predicted)])
test_GP_plot(GP_model=train_model_matern, test=test_data2[,c(rownames(param_ranges), predicted)]) 



