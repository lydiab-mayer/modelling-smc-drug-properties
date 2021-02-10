calcdifferencetrial  <- function(scenarios){

  for(i in 1:nrow(scenarios)){
    print(i)
seeds_smc <- read.table(paste0("/scicore/home/penny/GROUP/smc_lai/E5_1_CT_SMC/postprocessing/seeds_E5_1_CT_SMC_",scenarios[i,"seasonality"]
                               ,"_",4.9167,"_",scenarios[i,"SMC_HL"],"_",scenarios[i,"LAI_dec"],"_",scenarios[i,"EIR"],".txt"), header = T, as.is = TRUE, stringsAsFactors = FALSE)
seeds_lai <- read.table(paste0("/scicore/home/penny/GROUP/smc_lai/E5_2_CT_LAI/postprocessing/seeds_E5_2_CT_LAI_",scenarios[i,"seasonality"]
                               ,"_",4.9167,"__",scenarios[i,"LAI_dec"],"_",scenarios[i,"EIR"],".txt"), header = T, as.is = TRUE, stringsAsFactors = FALSE)
  
  
seeds <- length(unique(seeds_smc$seed))

seeds_smc <- seeds_smc[order(seeds_smc$seed),]
seeds_lai <- seeds_lai[order(seeds_lai$Scenario_Name,seeds_lai$seed),]

#remove runs where prevalence in year 5 is 0 (due to eliminating before switch)
scendellai <-seeds_lai[which( seeds_lai$prev_beg == 0),"Scenario_Name"]
scendelscm <- seeds_smc[which( seeds_smc$prev_beg == 0),"Scenario_Name"]
index <- c(scendellai,scendelscm)

if( length(index)>0){
  seeds_smc <- seeds_smc[- which(seeds_smc$Scenario_Name %in% index),]
  seeds_lai <- seeds_lai[- which(seeds_lai$Scenario_Name %in% index),]
}




seeds_smc <- do.call("rbind", replicate(nrow(seeds_lai)/seeds, seeds_smc, simplify = FALSE))

seeds_comparison <- seeds_lai
seeds_comparison$delta <- log(-log(seeds_lai[,"KM"])) -
  log(-log(seeds_smc[,"KM"]))



seeds_comparison$var <- (1/log(seeds_smc[,"KM"]))^2 * 1/(seeds_smc[,"KM"]^2)* seeds_smc[,"KM_var"] + 
  (1/log(seeds_lai[,"KM"]))^2 * 1/(seeds_lai[,"KM"]^2)* seeds_lai[,"KM_var"]
seeds_comparison$CI_low_delta <- seeds_comparison$delta-1.96* sqrt(seeds_comparison$var)
seeds_comparison$CI_high_delta <-seeds_comparison$delta+1.96* sqrt(seeds_comparison$var)

seeds_comparison$HR <- exp(seeds_comparison$delta)
seeds_comparison$CI_low_HR <- exp(seeds_comparison$CI_low_delta)

seeds_comparison$CI_high_HR <-exp(seeds_comparison$CI_high_delta)

seeds_comparison$Eff_SMC <- seeds_smc$KM
seeds_comparison$margin <- 0.05 
seeds_comparison$UL <- log(seeds_comparison$Eff_SMC-seeds_comparison$margin)/log(seeds_comparison$Eff_SMC)
seeds_comparison$non_inferiority <-   ifelse(seeds_comparison$CI_high_HR>seeds_comparison$UL,0,1) 
seeds_comparison$incred_SMC <- mean(seeds_smc$incred)
write.table(seeds_comparison, 
            file = paste0("/scicore/home/penny/GROUP/smc_lai/E5_Comparison/postprocessing/seeds_E5_comparison_" , scenarios[i,"seasonality"],"_",scenarios[i,"EIR"],"_",scenarios[i,"SMC_HL"],"_",scenarios[i,"LAI_dec"],".txt"),  sep = "\t", col.names = T, row.names = F)
 
seeds_comparison =   seeds_comparison %>% group_by(Scenario_Name) %>% 
  mutate(prop = sum(non_inferiority) / 10)

aggregated_OM =   seeds_comparison %>% group_by(Scenario_Name) %>% 
summarise_all(funs(if(is.numeric(.)) median(., na.rm = TRUE) else unique(.)))


write.table(aggregated_OM, 
            file = paste0("/scicore/home/penny/GROUP/smc_lai/E5_Comparison/postprocessing/agg_E5_comparison_", scenarios[i,"seasonality"],"_",scenarios[i,"EIR"],"_",scenarios[i,"SMC_HL"],"_",scenarios[i,"LAI_dec"],".txt"),  sep = "\t", col.names = T, row.names = F)


}
}
