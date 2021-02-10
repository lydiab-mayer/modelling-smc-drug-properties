
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
settings <- scenarios[-index,]

rows <- seq(1,364,25)


for (k in 1:nrow(settings)){
  files_exist <- 1
  
  seasonality <- settings[k,1]
  decay <- settings[k,2]
  age <- settings[k,3]
  acess <- settings[k,4]
  EIR <- settings[k,5]
  
  scen_name <- paste(seasonality,"_",age,"_",decay,"_",acess,"_",EIR, sep = "")
  file_scenarios_one <- paste('/scicore/home/penny/GROUP/smc_lai/E3E4disc_comparison/outfiles_optimisation_all/scenarios_', scen_name, '_', rows[1],'.txt', sep = "")
  
  if (file.exists(file_scenarios_one)){
  print(scen_name)
  outfile_scenario_assembled <- paste('/scicore/home/penny/GROUP/smc_lai/E3E4disc_comparison/outfiles_optimisation_summary/scenarios_all_',scen_name,".txt",sep = "")
  
  
  row <- rows[1]
  filename <- paste('/scicore/home/penny/GROUP/smc_lai/E3E4disc_comparison/outfiles_optimisation_all/scenarios_', scen_name, '_', row,'.txt', sep = "")
  
  if (file.exists(filename)){
    scenarios <- read.table(filename, header = T)
    
    
    for (j in 2:(length(rows)-1)){
      row <- rows[j]
      filename <- paste('/scicore/home/penny/GROUP/smc_lai/E3E4disc_comparison/outfiles_optimisation_all/scenarios_',scen_name,"_",row,".txt",sep = "")
      
      if (file.exists(filename)){
      a <- read.table(filename, header = T)
      scenarios[row:(row+25),] <- a[row:(row+25),]
      } else {
        files_exist <- 0
      }
      
      
    }
    
    row <- rows[length(rows)]
    filename <- paste('/scicore/home/penny/GROUP/smc_lai/E3E4disc_comparison/outfiles_optimisation_all/scenarios_',scen_name,"_",row,".txt",sep = "")
    
    if (file.exists(filename)){
    a <- read.table(filename, header = T)
    scenarios[row:nrow(scenarios),] <- a[row:nrow(scenarios),]
    } else{
      file_exists <- 0
    }
    
  } else{
    files_exist <- 0
  }
    

  } else{
    files_exist <- 0
  }
  
  if (files_exist == 1){
  write.table(scenarios, file = outfile_scenario_assembled)
  }
}
