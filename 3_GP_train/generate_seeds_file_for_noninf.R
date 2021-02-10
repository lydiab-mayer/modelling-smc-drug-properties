# generate a table to run the GP on the difference between two columns in seeds tables from different experiments
as_folder = "/scicore/home/smith/GROUP/smc_lai/Test_SMC/gp_5/as/"
generate_noninf_file(as_folder,processing_results1,processing_results2 ,cv_result$seasonality, cv_result$biting_pattern,cv_result$Decay_Scen,cv_result$maxGroup,cv_result$Access)

decay = "hill"
generate_noninf_file <- function(seasonality,biting_pattern,Decay_Scen,maxGroup,Access){
seeds_smc_smc <- processing_results2
seeds_smc_smc <- seeds_smc_smc[order(seeds_smc_smc$Scenario_Name, seeds_smc_smc$seed),]
seeds_smc_lai <- processing_results1
seeds_smc_lai <- seeds_smc_lai[order(seeds_smc_lai$Scenario_Name, seeds_smc_lai$seed),]

n_seeds <- length(unique(seeds_smc_lai$seed))

#ml: somehow the seeds got shuffeled in the postprocessing. Solve this here.
#get back to this later.
seeds_smc_smc$seed <- rep(1:n_seeds, nrow(seeds_smc_smc)/n_seeds)
seeds_smc_lai$seed <- rep(1:n_seeds, nrow(seeds_smc_lai)/n_seeds)

seeds_smc_smc$scen_id <- rep(1:(nrow(seeds_smc_smc)/n_seeds),each = n_seeds)
seeds_smc_lai$scen_id <- rep(1:(nrow(seeds_smc_lai)/n_seeds),each = n_seeds)

#remove runs where prevalence in year 5 is 0 (due to eliminating before switch)
scen_id_del_smc <- seeds_smc_smc[which( seeds_smc_smc$iprev_y5_1 < 0.01),"scen_id"]
scen_id_del_lai <- seeds_smc_lai[which( seeds_smc_lai$iprev_y5_1 < 0.01),"scen_id"]

index <- c(scen_id_del_smc, scen_id_del_lai)

seeds_smc_smc_reduced <- seeds_smc_smc#[- which(seeds_smc_smc$scen_id %in% index), ]
seeds_smc_lai_reduced <- seeds_smc_lai#[- which(seeds_smc_lai$scen_id %in% index), ]

seeds_comparison <- seeds_smc_lai_reduced
seeds_comparison$abs_diff_pppy_y10_all <- seeds_smc_lai_reduced$pppy_y10_all - seeds_smc_smc_reduced$pppy_y10_all
seeds_comparison$rel_diff_pppy_y10_all <- (seeds_smc_lai_reduced$pppy_y10_all - seeds_smc_smc_reduced$pppy_y10_all)/
  seeds_smc_smc_reduced$pppy_y10_all

# calculate non-inferiority criteria

seeds_comparison$delta <- log(-log(seeds_smc_lai_reduced[,"KM"])) -
  log(-log(seeds_smc_smc_reduced[,"KM"]))


seeds_comparison$var <- (1/log(seeds_smc_smc_reduced[,"KM"]))^2 * 1/(seeds_smc_smc_reduced[,"KM"]^2)* seeds_smc_smc_reduced[,"KM_var"] + 
  (1/log(seeds_smc_lai_reduced[,"KM"]))^2 * 1/(seeds_smc_lai_reduced[,"KM"]^2)* seeds_smc_lai_reduced[,"KM_var"]
seeds_comparison$CI_low_delta <- seeds_comparison$delta-1.96* sqrt(seeds_comparison$var)
seeds_comparison$CI_high_delta <-seeds_comparison$delta+1.96* sqrt(seeds_comparison$var)

seeds_comparison$HR <- exp(seeds_comparison$delta)
seeds_comparison$CI_low_HR <- exp(seeds_comparison$CI_low_delta)

seeds_comparison$CI_high_HR <-exp(seeds_comparison$CI_high_delta)

seeds_comparison$Eff_SMC <- seeds_smc_smc_reduced$KM

seeds_comparison$margin <-  0.05 

seeds_comparison$UL <- log(seeds_comparison$Eff_SMC-seeds_comparison$margin)/log(seeds_comparison$Eff_SMC)
seeds_comparison$non_inferiority <-   ifelse(seeds_comparison$CI_high_HR>seeds_comparison$UL,0,1)


agg_comparison_file <- paste0("/scicore/home/smith/GROUP/smc_lai/E3_E4_comparison/postprocessing_5/agg_E3_E4_comparison_seasonal_Low_indoor_",decay,".txt")
seeds_comparison =   seeds_comparison %>% group_by(Scenario_Name) %>% 
  mutate(prop = sum(non_inferiority) / 20)

aggregated_OM =   seeds_comparison %>% group_by(Scenario_Name) %>% 
  summarise_all(funs(if(is.numeric(.)) median(., na.rm = TRUE) else unique(.)))
file.remove(agg_comparison_file)


#write.table(aggregated_OM, file = agg_comparison_file, sep = "\t", col.names = T, row.names = F)


return(aggregated_OM)

}