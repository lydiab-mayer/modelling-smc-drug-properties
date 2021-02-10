# generate a table to run the GP on the difference between two columns in seeds tables from different experiments

Access = c(0.1)
Seasonality = c("Mali")
LAI_dec = c("exp")
Access = c(0.1)
EIR=3
#Seasonality <- c("Sen","Mali")
#LAI_dec <- c("hill","exp")
IntAge <- c(4.9167)
scenarios <- expand.grid(Access,Seasonality,LAI_dec,IntAge,EIR)
names(scenarios) <- c("Access","Seasonality","LAI_dec","IntAge","EIR")


ranges_file = "/scicore/home/penny/GROUP/smc_lai/E3_E4_comparison/param_ranges.RData"

param_ranges = rbind(c(1, 250), c(0.4, 1),c(0.4, 1), c(30, 150), c(0.7, 1))
row.names(param_ranges) = c("EIR", "Coverage","Coverage_SMC_Res" ,"Halflife", "Efficacy")
save(param_ranges, file = ranges_file)




for (i in 1:nrow(scenarios)){ 
scenario <- scenarios[i, ]
seeds_smc_smc <- data.frame(read.table(paste("/scicore/home/penny/GROUP/smc_lai/E4_SMCLAIdisc/postprocessing_5/seeds_E4_SMCLAIdisc_",scenario[,"Seasonality"],"_",scenario[,"IntAge"],"_",scenario[,"LAI_dec"],"_",scenario[,"Access"] ,"_3.txt", sep = ""), header = T, as.is = TRUE, stringsAsFactors = FALSE))

seeds_smc_smc <- seeds_smc_smc[order(seeds_smc_smc$Scenario_Name, seeds_smc_smc$seed),]
seeds_smc_lai <- data.frame(read.table(paste("/scicore/home/penny/GROUP/smc_lai/E4_testSMCdis/postprocessing_5/seeds_E4_testSMCdis_",scenario[,"Seasonality"],"_",scenario[,"IntAge"],"_",scenario[,"LAI_dec"],"_",scenario[,"Access"] ,"_3.txt", sep = ""), header = T, as.is = TRUE, stringsAsFactors = FALSE))
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

seeds_smc_smc_reduced <- seeds_smc_smc[- which(seeds_smc_smc$scen_id %in% index), ]
seeds_smc_lai_reduced <- seeds_smc_lai[- which(seeds_smc_lai$scen_id %in% index), ]

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


seeds_comparison_file <- paste0("/scicore/home/penny/GROUP/smc_lai/E3_E4_testSMCdis/postprocessing_5/seeds_E3_E4_testSMCdis_",scenario[,"Seasonality"],"_",scenario[,"IntAge"],"_",scenario[,"LAI_dec"],"_",scenario[,"Access"] ,".txt")

file.remove(seeds_comparison_file)
write.table(seeds_comparison, file = seeds_comparison_file,
			    sep = "\t", col.names = T, row.names = F)

agg_comparison_file <- paste0("/scicore/home/penny/burlyd00/smc_lai/E3_E4_testSMCdis/postprocessing_5/agg_E3_E4_testSMCDis_",scenario[,"Seasonality"],"_",scenario[,"IntAge"],"_",scenario[,"LAI_dec"],"_",scenario[,"Access"] ,".txt")
seeds_comparison$seeds_group <- ifelse(seeds_comparison$seed %in% seq(1,5),1,0)
seeds_comparison =   seeds_comparison %>% group_by(Scenario_Name, seeds_group) %>% 
  mutate(prop = sum(non_inferiority) / 5)

aggregated_OM =   seeds_comparison %>% group_by(Scenario_Name,seeds_group) %>% 
  summarise_all(funs(if(is.numeric(.)) median(., na.rm = TRUE) else unique(.)))
file.remove(agg_comparison_file)


write.table(aggregated_OM, file = agg_comparison_file,
            sep = "\t", col.names = T, row.names = F)

}
