import_procsim <- function(Exp,follow_up,setting,agg,scen_id){
  
  if(agg==TRUE){
    res <- read.table(file=paste0("/scicore/home/penny/GROUP/smc_lai/",Exp,"/postprocessing_",follow_up,"/agg_",Exp,"_",
                                  setting,".txt"),header=TRUE)
  } else {
    res <- read.table(file=paste0("/scicore/home/penny/GROUP/smc_lai/",Exp,"/postprocessing_",follow_up,"/seeds_",Exp,"_",
                                  setting,".txt"),header=TRUE)
    
  }
  
  if (scen_id[1]!="all" ) {res <- res[scen_id, ]}
  
  # colums with parameters
  
colparams  <- which(names(res)== "sum_clin")-1
 
 # extract clinical cases
df <- res
colInt <- grep(pattern="clin_y",colnames(df))
  names(df)[colInt] <- seq(1,10)
  
  ClinCases <-   gather(df, year, cases, colInt, factor_key=TRUE)
  ClinCases <-  ClinCases[, c(names(ClinCases)[1:colparams], "year","cases") ]
  
  rm(df)
  
  # extract prevalence reduction
  df <- res
 colInt <- grep(pattern="pppy",colnames(res))
  names(df)[colInt] <- seq(1,3)
  
  meanpppy <- data.frame(cbind(mean=rowMeans(df[,colInt]),
                               df[,c(names(df)[1:colparams])]))
  
  CasesPppy <-   gather(df, ageGroup, pppy, colInt, factor_key=TRUE)
  CasesPppy <-  CasesPppy[, c(names(CasesPppy)[1:colparams], "ageGroup","pppy") ]
  
  
  
  rm(df)
  
  
out <- list(ClinCases=ClinCases,
            meanpppy=meanpppy,
            CasesPppy=CasesPppy
           )
  return(out)
  
}



import_sim_trial <- function(Exp,EIR,SMC_HL,LAI_dec,scen_id){
  
  source('~/smc_lai/analysis_workflow/2_postprocessing/postprocessing_resources.R')
  
i=scen_id
  if (Exp== "E5_1_CT_SMC") { param_table_file =  paste0("/scicore/home/penny/GROUP/smc_lai/",Exp,"/postprocessing/split/",
                                      Exp,"_seasonal_Low_indoor_",EIR,"_",SMC_HL,".txt")} else {
                                        param_table_file =  paste0("/scicore/home/penny/GROUP/smc_lai/",Exp,"/postprocessing/split/",
                                                                   Exp,"_seasonal_Low_indoor_",EIR,"__",LAI_dec,".txt")   
                                      }

  
  om_results_folder = paste0("/scicore/home/penny/GROUP/smc_lai/",Exp,"/om/")
  
  param_table = read.table(param_table_file, sep= "\t", as.is = TRUE, header = TRUE, stringsAsFactors = FALSE)
  
 
#  param_table <- param_table[order(param_table$Scenario_Name), ]
  OM_result_file = paste(om_results_folder, param_table[i,]$Scenario_Name, "_",
                         +                                param_table[i,]$SEED, "_out.txt", sep="")
  OM_result = read.table(OM_result_file, sep="\t")
  
  #OM_result = read.table("./scicore/home/penny/GROUP/smc_lai/E5_2_CT_LAI/om/E5_2_CT_LAI.txt", sep="\t")
  
  om_result <- OM_result
  
  scenario_params <- param_table[i,]
  
  start= 1
    end= 73
  
 out <-  calculate_trial_outputs_plots(om_result, scenario_params,start,end)
  out$params <- scenario_params
 return(out)
  
}

import_EIRs <- function(setting){
  


     param_table_file =  paste("/scicore/home/penny/GROUP/smc_lai/E0_Cont/postprocessing_5/split/E0_Cont_",setting[, "seasonality"],"_4.9167_",
                               setting[, "LAI_dec"],"_",setting[, "Access"],'.txt', sep = "")
  
  om_results_folder = paste0("/scicore/home/penny/GROUP/smc_lai/E0_Cont/om/")
  
  param_table = read.table(param_table_file, sep= "\t", as.is = TRUE, header = TRUE, stringsAsFactors = FALSE)
  res_df <- list()
  for (i in 1:nrow(param_table)) {
  
  #  param_table <- param_table[order(param_table$Scenario_Name), ]
  OM_result_file = paste(om_results_folder, param_table[i,]$Scenario_Name, "_",
                         +                                param_table[i,]$SEED, "_cts.txt", sep="")
  OM_result = read.table(OM_result_file, sep="\t",header=TRUE)
  
  res <- subset(OM_result, timestep %in% seq(1,73))
  inpmax <- which(res$input.EIR==max(res$input.EIR))
  simmax <- which(res$simulated.EIR==max(res$simulated.EIR))
  diff <- simmax-inpmax
  
  
  res_df[[i]]<- data.frame(cbind(param_table[i,],simmax,inpmax,diff))
  
  }
  df <- do.call(rbind, res_df)
  return(df)
  
}

import_prev <- function(Exp,scen_id,param_table_file){
 follow_up =5

  
  

  
  
  om_results_folder = paste0("/scicore/home/penny/GROUP/smc_lai/",Exp,"/om/")
  
  param_table = read.table(param_table_file, sep= "\t", as.is = TRUE, header = TRUE, stringsAsFactors = FALSE)
  
  #  param_table <- param_table[order(param_table$Scenario_Name), ]
  OM_result_file = paste(om_results_folder, param_table[i,]$Scenario_Name, "_",
                         +                                param_table[i,]$SEED, "_out.txt", sep="")
  OM_result = read.table(OM_result_file, sep="\t",header=TRUE)
  
  #OM_result = read.table("./scicore/home/penny/GROUP/smc_lai/E5_2_CT_LAI/om/E5_2_CT_LAI.txt", sep="\t")
  
  om_result <- OM_result
  
  scenario_params <- param_table[i,]
  colnames(om_result) = c("time", "age_group", "measure", "value")
  year_to_5day = 73
  month_to_5day = 6
  years_before_interv = 5
  
  # Remove first measurement as it includes all cases until then
  to_remove = which(om_result$time == 1)
  om_result = om_result[-to_remove,]
  
  # population per time step per age group
  total_pop = as.data.frame(om_result[om_result$measure == 0, ])
  
  # add years tot time-steps, 73 time-steps are one year 
  total_pop$year <- rep(seq(1,max(total_pop$time)/year_to_5day),each=year_to_5day*max(om_result$age_group) )
  
  # calculate total population per age group per year
  total_pop_age = data.frame(total_pop %>% group_by(age_group, year) %>% summarise(n = mean(value) ))
  
  # sum up total population per year 
  total_pop_all = total_pop %>% group_by( time,year) %>% summarise(sum = sum(value) )%>% group_by( year)  %>% summarise(n = mean(sum) )
  
  # sum up population intervention age-groups over the years 
  pop_int <- total_pop_age[total_pop_age$age_group %in% seq(2,as.numeric(scenario_params["maxGroup"])),]%>% group_by( year) %>% summarise(n = sum(n) )
  # mean intervention population size 
  meanpopint <- mean(pop_int$n)
  
  
  # Summarize all measures results by summing up over age groups
  agg_om_result_total = om_result[,-which(names(om_result)=="age_group")]%>% group_by(time, measure) %>% summarise(val = sum(value))
  
  
  # Extract the clinical case numbers 
  trialpop <- sum( total_pop_age[total_pop_age$age_group %in% seq(2,as.numeric(scenario_params["maxGroup"])) & total_pop_age$year == 1,"n"])
  
  ######################################
  #calculate output prevalence reduction
  ######################################   
  
  # Calculate the prevalence for all the monitored years 
  # Prevalence = total number of infected people (in age group)/ total population (in age group)
  
  # number of infected per age-group per time-step
  n_infected = as.data.frame(om_result[om_result$measure == 1,])
  n_infected$year <- rep(seq(1,max(n_infected$time)/year_to_5day),each=year_to_5day*max(om_result$age_group) )
  
  # total number of infected per time-step
  n_infected_total = as.data.frame(agg_om_result_total[agg_om_result_total$measure == 1,])
  n_infected_total$year = rep(seq(1,max(n_infected_total$time)/year_to_5day),each=year_to_5day )
  
  # calculate prevalence 
  
  # add the population size in respective age-groups
  n_infected <- inner_join(total_pop_age,n_infected, by=c("age_group","year"))
  n_infected_total <- inner_join(total_pop_all,n_infected_total ,by=c("year"))
  
  # divide number of cases by respective age-group
  
  prev = n_infected %>% mutate(prev = value/n ) 
  
}

import_cpppy <- function(Exp,scen_ids,param_table_file){
  follow_up =5
  

  om_results_folder = paste0("/scicore/home/penny/GROUP/smc_lai/",Exp,"/om/")
  param_table_files <-paste0("/scicore/home/penny/GROUP/smc_lai/",Exp,"/postprocessing_5/",param_table_file)
  param_table = read.table(param_table_files, sep= "\t", as.is = TRUE, header = TRUE, stringsAsFactors = FALSE)
 
  survival_list <- list()
  df_list <- list()
  
  
  
  for (j in 1:length(scen_ids)) { 
   
  i <- scen_ids[j]
  
  
  #  param_table <- param_table[order(param_table$Scenario_Name), ]
  OM_result_file = paste(om_results_folder, param_table[i,]$Scenario_Name, "_",
                         +                                param_table[i,]$seed, "_out.txt", sep="")
  OM_result = read.table(OM_result_file, sep="\t",header=TRUE)
  
  #OM_result = read.table("./scicore/home/penny/GROUP/smc_lai/E5_2_CT_LAI/om/E5_2_CT_LAI.txt", sep="\t")
  
  om_result <- OM_result
  
  scenario_params <- param_table[i,]
  colnames(om_result) = c("time", "age_group", "measure", "value")
  year_to_5day = 73
  month_to_5day = 6
  years_before_interv = 5
  
  # Remove first measurement as it includes all cases until then
  to_remove = which(om_result$time == 1)
  om_result = om_result[-to_remove,]
  
  # population per time step per age group
  total_pop = as.data.frame(om_result[om_result$measure == 0, ])
  
  # add years tot time-steps, 73 time-steps are one year 
  total_pop$year <- rep(seq(1,max(total_pop$time)/year_to_5day),each=year_to_5day*max(om_result$age_group) )
  
  # calculate total population per age group per year
  total_pop_age = data.frame(total_pop %>% group_by(age_group, year) %>% summarise(n = mean(value) ))
  
  # sum up total population per year 
  total_pop_all = total_pop %>% group_by( time,year) %>% summarise(sum = sum(value) )%>% group_by( year)  %>% summarise(n = mean(sum) )
  
  # sum up population intervention age-groups over the years 
  pop_int <- total_pop_age[total_pop_age$age_group %in% seq(2,as.numeric(scenario_params["maxGroup"])),]%>% group_by( year) %>% summarise(n = sum(n) )
  # mean intervention population size 
  meanpopint <- mean(pop_int$n)
  
  
  # Summarize all measures results by summing up over age groups
  agg_om_result_total = om_result[,-which(names(om_result)=="age_group")]%>% group_by(time, measure) %>% summarise(val = sum(value))
  
  
  # Extract the clinical case numbers 
  trialpop <- sum( total_pop_age[total_pop_age$age_group %in% seq(2,as.numeric(scenario_params["maxGroup"])) & total_pop_age$year == 1,"n"])
  
  ######################################
  #calculate output prevalence reduction
  ######################################   
  nclin = as.data.frame(om_result[om_result$measure == 14,])
  nclin$year <- c(rep(seq(1,max(nclin$time)/year_to_5day),each=year_to_5day*max(om_result$age_group) ))
  
  
  # yearly average clinical cases per age-group
  nclin_yearly <- nclin %>% group_by(age_group ,year) %>% summarise(sum = sum(value) )
  
  
  # summed yearly clinical cases in intervention age-groups
  final_nclin_intGroup <- nclin_yearly[nclin_yearly$age_group  %in% c(2,as.numeric(scenario_params["maxGroup"])),]
  
  # sum of intervntion age-groups
  final_nclin_yearly <- final_nclin_intGroup[final_nclin_intGroup$year %in% seq(1,years_before_interv+follow_up),] %>% 
    group_by( year) %>% 
    summarise(sumclin = sum(sum) )
  
  
  
  final_nclin_pppy <- final_nclin_yearly[final_nclin_yearly$year ==10,"sumclin"] /pop_int[pop_int$year ==10,"n"]$n 
  beg_nclin_pppy <- final_nclin_yearly[final_nclin_yearly$year ==4,"sumclin"] /pop_int[pop_int$year ==5,"n"]$n 
  
  
  # calculate cpppy in 2-10 years old for reference 
  
  
  final_nclin_210 <- nclin_yearly[nclin_yearly$age_group %in% c(3,4),]
  
  # sum of intervntion age-groups
  final_nclin_yearly210 <- final_nclin_210[final_nclin_210$year %in% seq(1,years_before_interv+follow_up),] %>% 
    group_by( year) %>% 
    summarise(sumclin = sum(sum) )
  

  # cases pppy at last survey in the interventaion age-groups
  final_nclin_age_group <- final_nclin_intGroup[final_nclin_intGroup$year == years_before_interv +follow_up,"sum"]/
    total_pop_age[total_pop_age$age_group %in% seq(2,as.numeric(scenario_params["maxGroup"])) & total_pop_age$year == years_before_interv +follow_up,"n"]
  
  # cumulative clinical cases over all age-groups over all years
  final_nclin = c(sum(final_nclin_yearly$sumclin))
  
  
  # extract cases in intervention year 
  nclintrial <-  nclin[nclin$age_group %in% c(2,as.numeric(scenario_params["maxGroup"]))  & nclin$year== years_before_interv+follow_up, ]
  nclintrial <- nclintrial %>% group_by(time) %>% summarise(sum = sum(value) )
  
  nclintrial$timeyeartrial <- nclintrial$time-min(nclintrial$time)
  nclintrial$cpp <- nclintrial$sum/trialpop
  
  
  nclinint <-  nclintrial
  nclinint$trialtime <-  nclinint$timeyeartrial-min(nclinint$timeyeartrial)
  
  nclinint <-  nclinint  %>% mutate(cppcum = cumsum(cpp))
  
  
  
  dfsurv <- nclinint[, c("trialtime", "sum")]
  
  # dfcases <- as.data.frame(lapply(dfsurv, rep, dfsurv$sum))
  # dfcases$surv <- 1
  # dfcases$id <- seq(1:nrow(dfcases))
  # 
  # dfcens <- data.frame(cbind(trialtime=year_to_5day, sum=0,surv=0, id= seq(max(dfcases$id)+1:trialpop) ))
  # 
  # dfsurv2 <- data.frame(rbind(dfcases, dfcens))
  # 
  # 
  # km_fit <- survfit(Surv(trialtime, surv) ~ 1, data=dfsurv2)
  # surv <- summary(km_fit, times = max(dfsurv$trialtime))
  # survival <-  surv$surv
  # survival_var <- (surv$std.err)^2
  
  calc_KM <- function(df,trialpop) {
    
    df$natrisk <- trialpop-df$sum
    
    processed <- df %>%
      arrange(trialtime) %>%
      mutate(KM = cumprod((natrisk - sum) / natrisk),
             KM_se = KM*sqrt(cumsum(sum/(natrisk * (natrisk-sum)) ) ) )
    
    
    return(processed)
  }
  km_fit <- calc_KM(dfsurv, trialpop)
  
  survival_list[[j]] <- km_fit
  
  df_list[[j]] <-  nclinint 
  
  }
  
  dfs <- rbindlist(df_list)[,lapply(.SD,mean), list(timeyeartrial, trialtime)]  
  
  dfs2 <- rbindlist(survival_list)[,lapply(.SD,mean), list( trialtime)]  
  
  return(list(dfs,dfs2 ))
} 
