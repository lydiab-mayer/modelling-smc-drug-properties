rm(list = ls())

seasonalities <- c("Mali","Sen")
decays <- c("exp","hill","wei")
ages <- c(4.9167)
acesses <- c(0.1,0.5)
EIRs <- c(3,4,5,8,9,20,28,47,150)

EIRs_low_acces <- c(3,4,8,28)
EIRs_high_acces <- c(5,9,20,47)


scenarios <- expand.grid(seasonalities,decays,ages,acesses,EIRs)
colnames(scenarios) <- c("Seasonality","Decay","Age","Acess","EIR")

index <- which(scenarios$Acess == 0.1 & scenarios$EIR %in% EIRs_high_acces)
scenarios <- scenarios[-index,]

index <- which(scenarios$Acess == 0.5 & scenarios$EIR %in% EIRs_low_acces)
scenarios <- scenarios[-index,]


for (k in 1:nrow(scenarios)){
  
  seasonality <- scenarios[k,1]
  decay <- scenarios[k,2]
  age <- scenarios[k,3]
  acess <- scenarios[k,4]
  EIR <- scenarios[k,5]
  
  scen_name <- paste(seasonality,"_",age,"_",decay,"_",acess,"_",EIR, sep = "")
  
  
  infile_SMC <- paste('/scicore/home/penny/GROUP/smc_lai/E3_SMCSMCdisc/postprocessing_5/seeds_E3SMCSMCdisc_', scen_name, '.txt', sep = "")
  infile_LAI <- paste('/scicore/home/penny/GROUP/smc_lai/E4_SMCLAIdisc/postprocessing_5/seeds_E4SMCLAIdisc_', scen_name, '.txt', sep = "")
  
  outfile <- paste("/scicore/home/penny/GROUP/smc_lai/E3E4disc_comparison/postprocessing_5/seeds_E3E4disc_comparison_", scen_name,".txt", sep = "")
  
  #print(c(scen_name,file.exists(infile_SMC),  file.exists(infile_LAI)))

  if (file.exists(infile_SMC) & file.exists(infile_LAI)){

  seeds_SMC <- read.table(infile_SMC, header = T)
  seeds_LAI <- read.table(infile_LAI, header = T)
  
  print(c(scen_name, nrow(seeds_SMC), nrow(seeds_LAI)))
  
  if (nrow(seeds_SMC) == 7500 & nrow(seeds_LAI) == 7500){
  
  ##### generate difference file
  n_seeds <- length(unique(seeds_SMC$seed))
  
  seeds_SMC$scen_id <- rep(1:(nrow(seeds_SMC)/n_seeds),each = n_seeds)
  seeds_LAI$scen_id <- rep(1:(nrow(seeds_LAI)/n_seeds),each = n_seeds)
  
  #remove runs where prevalence in year 5 is 0 (due to eliminating before switch)
  scen_id_del_smc <- seeds_SMC[which(seeds_SMC$iprev_y5_2 < 0.01),"scen_id"]
  scen_id_del_lai <- seeds_LAI[which(seeds_LAI$iprev_y5_2 < 0.01),"scen_id"]
  
  index <- c(scen_id_del_smc, scen_id_del_lai)
  seeds_SMC_reduced <- seeds_SMC
  seeds_LAI_reduced <- seeds_LAI
  
  if (length(index) > 0){
    seeds_SMC_reduced <- seeds_SMC_reduced[- which(seeds_SMC_reduced$scen_id %in% index), ]
    seeds_LAI_reduced <- seeds_LAI_reduced[- which(seeds_LAI_reduced$scen_id %in% index), ]
  }
  
  seeds_comparison <- seeds_LAI_reduced
  seeds_comparison$abs_diff_pppy_y10_all <- seeds_LAI_reduced$pppy_y10_all - seeds_SMC_reduced$pppy_y10_all
  seeds_comparison$rel_diff_pppy_y10_all <- (seeds_LAI_reduced$pppy_y10_all - seeds_SMC_reduced$pppy_y10_all)/
  seeds_SMC_reduced$pppy_y10_all
  
  # calculate non-inferiority criteria
  seeds_comparison$delta <- log(-log(seeds_LAI_reduced[,"KM"])) - log(-log(seeds_SMC_reduced[,"KM"]))
  
  seeds_comparison$var <- (1/log(seeds_SMC_reduced[,"KM"]))^2 * 1/(seeds_SMC_reduced[,"KM"]^2)* seeds_SMC_reduced[,"KM_var"] + 
  (1/log(seeds_LAI_reduced[,"KM"]))^2 * 1/(seeds_LAI_reduced[,"KM"]^2)* seeds_LAI_reduced[,"KM_var"]
  
  seeds_comparison$CI_low_delta <- seeds_comparison$delta-1.96* sqrt(seeds_comparison$var)
  seeds_comparison$CI_high_delta <-seeds_comparison$delta+1.96* sqrt(seeds_comparison$var)

  
  seeds_comparison$HR <- exp(seeds_comparison$delta)
  seeds_comparison$CI_low_HR <- exp(seeds_comparison$CI_low_delta)
  
  seeds_comparison$CI_high_HR <-exp(seeds_comparison$CI_high_delta)
  
  seeds_comparison$Eff_SMC <- seeds_SMC_reduced$KM
  seeds_comparison$margin <-  0.05 
  
  seeds_comparison$UL <- log(seeds_comparison$Eff_SMC-seeds_comparison$margin)/log(seeds_comparison$Eff_SMC)
  
  #in some (very few) cases, Eff_SMC < margin, leading to log of negative numbers (NaN). 
  #in these cases, non inferiority would not be achieved
  #new: drop these simulations as they lead to bad GP performance
  index <- which(is.nan(seeds_comparison$UL))
  if (length(index) > 0){
  #seeds_comparison[index ,]$UL <- 0
  seeds_comparison <- seeds_comparison[-index,]
  }
  
  seeds_comparison$non_inferiority <-   ifelse(seeds_comparison$CI_high_HR>seeds_comparison$UL,0,1)
  
  write.table(seeds_comparison, file = outfile, sep = "\t", col.names = T, row.names = F)
  
  }
  }
}




