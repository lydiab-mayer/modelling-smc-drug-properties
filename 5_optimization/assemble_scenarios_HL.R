
rm(list = ls())

seasonalities <- c("Mali")
decays <- c("exp","hill")
ages <- c(4.9167)
acesses <- c(0.1)
EIRs <- c(3,4,8,28,150)

EIRs_low_acces <- c(3,4,8,28)
EIRs_high_acces <- c(5,9,20,47)


scenarios <- expand.grid(seasonalities,decays,ages,acesses,EIRs)
colnames(scenarios) <- c("Seasonality","Decay","Age","Acess","EIR")


settings <- scenarios

rows <- seq(1,364,25)

list_df <- list()
for (k in 1:nrow(settings)){

  seasonality <- settings[k,1]
  decay <- settings[k,2]
  age <- settings[k,3]
  acess <- settings[k,4]
  EIR <- settings[k,5]
  
  scen_name <- paste(seasonality,"_",age,"_",decay,"_",acess,"_",EIR, sep = "")
  file_scenarios_one <- paste('/scicore/home/penny/GROUP/smc_lai/E3E4_lb/outfiles_optimisation_HL/scenarios_', scen_name,".txt", sep = "")
  
  outfile_scenario_assembled <- paste('/scicore/home/penny/GROUP/smc_lai/E3E4_lb/outfiles_optimisation_summary_HL/scenarios_all_',scen_name,".txt",sep = "")
  
  
  filename <- paste('/scicore/home/penny/GROUP/smc_lai/E3E4_lb/outfiles_optimisation_HL/scenarios_', scen_name,'.txt', sep = "")
  
  if (file.exists(filename)){
    scenarios <- read.table(filename, header = T)
  }
  scenarios <- data.frame(cbind(scenarios,settings[k,]))
  list_df[[k]] <- scenarios
 
}
