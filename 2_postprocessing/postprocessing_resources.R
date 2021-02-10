##########################
# Auxiliary functions for postprocessing OpenMalaria results
#
# created 02.10.2019
# monica.golumbeanu@unibas.ch
# modified by lydia.burgert@unibas.ch 
##########################
library(rapportools)
library(survival)
library(cmprsk)

# Function which calculates prevalence reduction given an OpenMalaria
# simulation result.
calculate_outputs = function(om_result, scenario_params, follow_up) {
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
 
  pop_210 <- total_pop_age[total_pop_age$age_group %in% c(3,4),]%>% group_by( year) %>% summarise(n = sum(n) )
  
   # mean intervention population size 
  meanpopint <- mean(pop_int$n)
  
  
  # Summarize all measures results by summing up over age groups
  agg_om_result_total = om_result[,-which(names(om_result)=="age_group")]%>% group_by(time, measure) %>% summarise(val = sum(value))
  
  
  # Extract the clinical case numbers 
  trialpop <- sum( total_pop_age[total_pop_age$age_group %in% seq(2,as.numeric(scenario_params["maxGroup"])) & total_pop_age$year == years_before_interv +follow_up,"n"])
  
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
  prev_total = n_infected_total %>% mutate(prev = val/n )
  
  cases_int_y10 <- data.frame(n_infected[n_infected$year==10 & n_infected$age_group %in% seq(2,as.numeric(scenario_params["maxGroup"])),]%>% 
    group_by(time)%>% mutate(totcas = sum(value) ))
  prev_int_y10 <- mean(cases_int_y10[cases_int_y10$age_group==2 , "totcas"] /  pop_int[pop_int$year ==10,"n"]$n )
  
  cases_int_y5 <- data.frame(n_infected[n_infected$year==5 & n_infected$age_group %in% seq(2,as.numeric(scenario_params["maxGroup"])),]%>% 
                                group_by(time)%>% mutate(totcas = sum(value) ))
  
  prev_int_y5 <- mean(cases_int_y5[cases_int_y5$age_group==2 , "totcas"] /  pop_int[pop_int$year ==10,"n"]$n )
  
 
  
  cases_210 <- data.frame(n_infected[n_infected$year==5 & n_infected$age_group %in% c(3,4),]%>% 
                                group_by(time)%>% mutate(totcas = sum(value) ))
  
  prev_210_y5 <- mean(cases_210[cases_210$age_group==4 , "totcas"] /  pop_210[pop_210$year ==5,"n"]$n )
  
  
  # yearly average prevalence
  yearly_avg_prev <- prev %>% group_by(age_group, year) %>% summarise(avg = mean(prev) )
  yearly_avg_prev_total = prev_total %>% group_by( year) %>% summarise(avg = mean(prev) )
  
  # initial prevalence for prevalence reduction
  initial_prev = yearly_avg_prev[yearly_avg_prev$year==years_before_interv,"avg"]
    #prevalence reduction for follow up point
  
  #Select the last time point according to the type of follow-up
  final_prev = yearly_avg_prev[yearly_avg_prev$year==years_before_interv + follow_up,"avg"]
  

  # prevalence reduction per age group
  prev_red = (initial_prev - final_prev)/initial_prev * 100
  prev_red = prev_red*(prev_red >= 0) 
  
  #Select the last time point according to the type of follow-up
  final_prevall = yearly_avg_prev[yearly_avg_prev$year==years_before_interv + follow_up,"avg"]
  
  # prevalence reduction per age group
  prev_red_all = (initial_prev - final_prev)/initial_prev * 100
  prev_red_all = prev_red*(prev_red >= 0) 
  
  # prevalence reduction in intervention age-group
  prev_red_int = (prev_int_y5 - prev_int_y10)/prev_int_y5 * 100
  prev_red_int = prev_red_int*(prev_red_int >= 0) 
  
  
  ######################################
  #calculate output clinical cases
  ######################################   
  
  # Extract the clinical case numbers 
  
  # number of clinical cases
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
  
  beg_nclin_pppy210 <- final_nclin_yearly210[final_nclin_yearly210$year ==4,"sumclin"] /pop_210[pop_210$year ==5,"n"]$n 
  
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
  
   survival <- tail(km_fit,1)$KM
   
   survival_var <- (tail(km_fit,1)$KM_se)^2
  # Final row with outputs to return
  return_row = cbind.data.frame(scenario_params$Scenario_Name, 
                                scenario_params$SEED,
                                final_nclin,meanpopint,final_nclin_pppy,beg_nclin_pppy,beg_nclin_pppy210,
                                t(final_nclin_yearly$sumclin),
                                t(final_nclin_age_group),
                                t(initial_prev[ 2:as.numeric(scenario_params["maxGroup"]),]) ,prev_210_y5,
                                t(prev_red[ 2:as.numeric(scenario_params["maxGroup"]),]), prev_int_y5,prev_red_int,survival,survival_var
  ) 
  colnames(return_row) = c("Scenario_Name",
                           "seed", 
                           "sum_clin","mean_popint","pppy_y10_all","pppy_y4_all","pppy_y4_210",
                           paste0("clin_y",final_nclin_yearly$year),
                           paste0("pppy_y", follow_up+years_before_interv,"_", seq(2,as.numeric(scenario_params["maxGroup"]))),
                           paste0("iprev_y",years_before_interv,"_", seq(2,as.numeric(scenario_params["maxGroup"]))),"prev_210_y5",
                           paste0("prevred_y", follow_up+years_before_interv,"_", seq(2,as.numeric(scenario_params["maxGroup"]))),
                           "iprev_int_y5","prevred_int_y10","KM","KM_var"
  )
  
  return(return_row)
}


# function which calculates outputs for a clinical trial run in OpenMalaria
calculate_trial_outputs = function(om_result, scenario_params,start,end) {
  colnames(om_result) = c("time", "age_group", "measure", "value")
  year_to_5day = 73
  month_to_5day = 6
  years_before_interv = 5
  follow_up=1
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
  
  # year of intervention : 2009 ---> year 10 
  
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
  prev_total = n_infected[n_infected$age_group %in% c(2,as.numeric(scenario_params["maxGroup"])), ]  %>% mutate(prev = value/n )
  
  
  # yearly average prevalence
  
  yearly_avg_prev <- prev %>% group_by(age_group, year) %>% summarise(avg = mean(prev) )
  
  prev_total = prev_total%>% group_by(time,year)  %>%summarise( prev = weighted.mean(prev,n)) 
  
  yearly_avg_prev_total = prev_total %>% group_by(year) %>% summarise(avg = mean(prev) )
  
  
  
  # prevalence during intervention
  prevtrial <- prev_total[prev_total$year==years_before_interv+follow_up, ]
  prevtrial$trialtime <- prevtrial$time-min(prevtrial$time)
  prevint <-  prevtrial[prevtrial$trialtime>=(start-1), ]
  
  ######################################
  #calculate output clinical cases
  ######################################   
  
  # Extract the clinical case numbers 
  trialpop <- sum( total_pop_age[total_pop_age$age_group %in% seq(2,as.numeric(scenario_params["maxGroup"])) & total_pop_age$year == years_before_interv +follow_up,"n"])
  # number of clinical cases
  nclin = as.data.frame(om_result[om_result$measure == 14,])
  nclin$year <- c(rep(seq(1,max(nclin$time)/year_to_5day),each=year_to_5day*max(om_result$age_group) ))
  
  
  # yearly average clinical cases per age-group
  nclin_yearly <- nclin %>% group_by(age_group ,year) %>% summarise(sum = sum(value) )
  # summed yearly clinical cases in intervention age-groups
  final_nclin_intGroup <- nclin_yearly[nclin_yearly$age_group <=3,]
  # sum of intervntion age-groups
  final_nclin_yearly <- final_nclin_intGroup[final_nclin_intGroup$year %in% seq(1,years_before_interv+follow_up),] %>% 
    group_by( year) %>% 
    summarise(sumclin = sum(sum) )
  
  

  
  
  
  
  # extract cases in intervention year 
  nclintrial <-  nclin[nclin$age_group %in% c(2,as.numeric(scenario_params["maxGroup"]))  & nclin$year== years_before_interv+follow_up, ]
  nclintrial <- nclintrial %>% group_by(time) %>% summarise(sum = sum(value) )
  
  nclintrial$timeyeartrial <- nclintrial$time-min(nclintrial$time)
  nclintrial$cpp <- nclintrial$sum/trialpop
  
  
  nclinint <-  nclintrial[nclintrial$timeyeartrial %in% seq(start,end), ]
  nclinint$trialtime <-  nclinint$timeyeartrial-min(nclinint$timeyeartrial)
  
  nclinint <-  nclinint  %>% mutate(cppcum = cumsum(cpp))

  
  prev_beg <- prevint[1,"prev"]
prev_end <- prevint[end-start,"prev"]
  CPP <- nclinint[end-start,"cppcum"]
  
dfsurv <- nclinint[seq(1:(end-start)),c("trialtime", "sum")]

calc_KM <- function(df,trialpop) {
  
  df$natrisk <- trialpop-df$sum
  
  processed <- df %>%
    arrange(trialtime) %>%
    mutate(KM = cumprod((natrisk - sum) / natrisk),
           KM_se = KM*sqrt(cumsum(sum/(natrisk * (natrisk-sum)) ) ) )
  
  
  return(processed)
}
km_fit <- calc_KM(dfsurv, trialpop)

survival <- tail(km_fit,1)$KM

survival_var <- (tail(km_fit,1)$KM_se)^2


#extract number of cases in year before clinical trial for trial period

nclinNI <-  nclin[nclin$age_group %in% c(2,as.numeric(scenario_params["maxGroup"])) & nclin$year== years_before_interv, ]
nclinNI <- nclinNI %>% group_by(time) %>% summarise(sum = sum(value) )

nclinNI$timeyeartrial <- nclinNI$time-min(nclinNI$time)
nclinNI$cpp <- nclinNI$sum/trialpop


nclinintNI <-  nclinNI[nclinNI$timeyeartrial %in% seq(start,end), ]
nclinintNI$trialtime <-  nclinintNI$timeyeartrial-min(nclinintNI$timeyeartrial)

nclinintNI <-  nclinintNI  %>% mutate(cppcum = cumsum(cpp))



CPPNI <- nclinintNI[end-start,"cppcum"]
incred <- (CPPNI-CPP)/CPPNI
 if(incred<0) {incred <- 0} 




#   dfcases <- as.data.frame(lapply(dfsurv, rep, dfsurv$sum))
#   dfcases$surv <- 1
#   dfcases$id <- seq(1:nrow(dfcases))
#   
#   dfcens <- data.frame(cbind(trialtime=end, sum=0,surv=0, id= seq(max(dfcases$id)+1:trialpop) ))
#   
#   dfsurv2 <- data.frame(rbind(dfcases, dfcens))
#   
# 
#   km_fit <- survfit(Surv(trialtime, surv) ~ 1, data=dfsurv2)
#   surv <- summary(km_fit, times = max(dfsurv$trialtime))
#  survival <-  surv$surv
#   survival_var <- (surv$std.err)^2
#   
#   
#   xx_cif<-cuminc(dfsurv2$trialtime, dfsurv2$surv,cencode=0)
# surv_cif <-   timepoints(xx_cif, max(dfsurv$trialtime))
#   
# cif <- as.numeric(surv_cif$est)
# cif_var <- as.numeric(surv_cif$var)
  ####################

  return_row = cbind.data.frame(scenario_params$Scenario_Name, 
                                scenario_params$SEED,
                                prev_beg, prev_end, CPP,CPPNI,incred,
                                survival, survival_var
  ) 
  colnames(return_row) = c("Scenario_Name",
                           "seed", 
                           "prev_beg","prev_end","cppcum","cppcumNI","incred",
                           "KM","KM_var"
  )
  return(return_row)
  

}

calculate_trial_outputs_plots = function(om_result, scenario_params,start,end) {
  colnames(om_result) = c("time", "age_group", "measure", "value")
  year_to_5day = 73
  month_to_5day = 6
  years_before_interv = 5
  follow_up=1
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
  
  # year of intervention : 2009 ---> year 10 
  
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
  prev_total = n_infected[n_infected$age_group %in% c(2,as.numeric(scenario_params["maxGroup"])), ]  %>% mutate(prev = value/n )
  
  
  # yearly average prevalence
  
  yearly_avg_prev <- prev %>% group_by(age_group, year) %>% summarise(avg = mean(prev) )
  
  prev_total = prev_total%>% group_by(time,year)  %>%summarise( prev = weighted.mean(prev,n)) 
  
  yearly_avg_prev_total = prev_total %>% group_by(year) %>% summarise(avg = mean(prev) )
  
  
  
  # prevalence during intervention
  prevtrial <- prev_total[prev_total$year==years_before_interv+follow_up, ]
  prevtrial$trialtime <- prevtrial$time-min(prevtrial$time)
  prevint <-  prevtrial[prevtrial$trialtime>=(start-1), ]
  
  ######################################
  #calculate output clinical cases
  ######################################   
  
  # Extract the clinical case numbers 
  trialpop <- sum( total_pop_age[total_pop_age$age_group %in% seq(2,as.numeric(scenario_params["maxGroup"])) & total_pop_age$year == years_before_interv +follow_up,"n"])
  # number of clinical cases
  nclin = as.data.frame(om_result[om_result$measure == 14,])
  nclin$year <- c(rep(seq(1,max(nclin$time)/year_to_5day),each=year_to_5day*max(om_result$age_group) ))
  
  
  # yearly average clinical cases per age-group
  nclin_yearly <- nclin %>% group_by(age_group ,year) %>% summarise(sum = sum(value) )
  # summed yearly clinical cases in intervention age-groups
  final_nclin_intGroup <- nclin_yearly[nclin_yearly$age_group <=2,]
  # sum of intervntion age-groups
  final_nclin_yearly <- final_nclin_intGroup[final_nclin_intGroup$year %in% seq(1,years_before_interv+follow_up),] %>% 
    group_by( year) %>% 
    summarise(sumclin = sum(sum) )
  
  # extract cases in intervention year 
  nclintrial <-  nclin[nclin$age_group %in% c(2,as.numeric(scenario_params["maxGroup"]))  & nclin$year== years_before_interv+follow_up, ]
  nclintrial <- nclintrial %>% group_by(time) %>% summarise(sum = sum(value) )
  
  nclintrial$timeyeartrial <- nclintrial$time-min(nclintrial$time)
  nclintrial$cpp <- nclintrial$sum/trialpop
  
  
  nclinint <-  nclintrial[nclintrial$timeyeartrial %in% seq(start,end), ]
  nclinint$trialtime <-  nclinint$timeyeartrial-min(nclinint$timeyeartrial)
  
  nclinint <-  nclinint  %>% mutate(cppcum = cumsum(cpp))
  
  
  prev_beg <- prevint[1,"prev"]
  prev_end <- prevint[end-start,"prev"]
  CPP <- nclinint[end-start,"cppcum"]
  
  dfsurv <- nclinint[seq(1:(end-start)),c("trialtime", "sum")]
  
  calc_KM <- function(df,trialpop) {
    
    df$natrisk <- trialpop-df$sum
    
    processed <- df %>%
      arrange(trialtime) %>%
      mutate(KM = cumprod((natrisk - sum) / natrisk),
             KM_se = KM*sqrt(cumsum(sum/(natrisk * (natrisk-sum)) ) ) )
    
    
    return(processed)
  }
  km_fit <- calc_KM(dfsurv, trialpop)
  
  survival <- tail(km_fit,1)$KM
  
  survival_var <- (tail(km_fit,1)$KM_se)^2

  return_row = list(Scenario_Name=scenario_params$Scenario_Name, 
                       seed=scenario_params$SEED,
                       prev_beg= prev_beg, 
                       prev_end=prev_end, 
                       cppcum=CPP,
                               KM= survival, 
                       KM_var=survival_var,
                       df_prev=prevint,
                       df_cas=nclinint
  ) 
    return(return_row)
  
  
}

# Wrapper for looping across all simulation results and gathering postprocessing results in a table
postprocess_OM = function(results_folder, param_table_file, final_table_dest, 
                          final_seed_table_dest, follow_up) {
  
  param_table = read.table(param_table_file, sep= "\t", as.is = TRUE, header = TRUE, stringsAsFactors = FALSE)
  processed_OM_sim = NULL
  for( i in 1:nrow(param_table)) {
    
    print(i)
    # Read the OM simulation result
    OM_result_file = paste(results_folder, param_table[i,]$Scenario_Name, "_", 
                           param_table[i,]$SEED, "_out.txt", sep="")
    # Calculate the necessary outputs
    if(file.exists(OM_result_file) & file.info(OM_result_file)$size > 0) {
      # if (str_detect(OM_result_file, "_1220_")) {
      #     write("detected")
      # }
      OM_result = read.table(OM_result_file, sep="\t")
      # om_result, scenario_params, total_pop, survey_start, survey_end, int_start, int_end, pulsed_int_start
      scenario_row = calculate_outputs(OM_result, param_table[i,], follow_up)
      processed_OM_sim = data.frame(rbind(processed_OM_sim, scenario_row),stringsAsFactors = FALSE)
    }
  }
  
  
  # Summarize results over seeds and create final results tables
  
  
  
  aggregated_OM =   processed_OM_sim %>% group_by(Scenario_Name) %>% summarise_at(c(names(processed_OM_sim)[which(names(processed_OM_sim)=="sum_clin"):length(names(processed_OM_sim) ) ]),median,na.rm=TRUE)
  
  
  
  no_seed_table = param_table[,-c(which(colnames(param_table)=="SEED"))]
  no_seed_table = unique(no_seed_table)
  final_seed_table = merge(no_seed_table, processed_OM_sim, by = c("Scenario_Name"))
  final_table = merge(no_seed_table, aggregated_OM, by = c("Scenario_Name"))
  
  # Write result tables (summarized and with seeds) to files
  write.table(final_table, final_table_dest, sep="\t",
              col.names = TRUE, row.names = FALSE, quote=FALSE)
  write.table(final_seed_table, final_seed_table_dest, sep="\t",
              col.names = TRUE, row.names = FALSE, quote=FALSE)
}



# Wrapper for looping across all simulation results and gathering postprocessing results in a table
postprocess_OM_trial = function(results_folder, param_table_file, final_table_dest, 
                          final_seed_table_dest) {
  
  param_table = read.table(param_table_file, sep= "\t", as.is = TRUE, header = TRUE, stringsAsFactors = FALSE)
  processed_OM_sim = NULL
  for( i in 1:nrow(param_table)) {
    # Read the OM simulation result
    OM_result_file = paste(results_folder, param_table[i,]$Scenario_Name, "_", 
                           param_table[i,]$SEED, "_out.txt", sep="")
    # Calculate the necessary outputs
    if(file.exists(OM_result_file) & file.info(OM_result_file)$size > 0) {
      # if (str_detect(OM_result_file, "_1220_")) {
      #     write("detected")
      # }
      OM_result = read.table(OM_result_file, sep="\t")
      # om_result, scenario_params, total_pop, survey_start, survey_end, int_start, int_end, pulsed_int_start
      if(param_table[i,"Seasonality"]=="Mali"){
      timing_df = data.frame(EIR=c(5,9,20,47,150),
                             timingInt=c(5,5,5,4,3),
                             start= c(47,47,47,46,45),
                             end=c(71,71,71,70,69))
      } else {
        timing_df = data.frame(EIR=c(5,9,20,47,150),
                               timingInt=c(5,5,5,4,3),
                               start= c(53,53,53,52,51),
                               end=c(71,71,71,70,69))
        
     }
      
      start <- timing_df[which(param_table[i,"EIR"]==timing_df[,"EIR"]),"start"]
      end <- timing_df[which(param_table[i,"EIR"]==timing_df[,"EIR"]),"end"]
      
      
        scenario_row = calculate_trial_outputs(OM_result, param_table[i,], start,end)
      processed_OM_sim = data.frame(rbind(processed_OM_sim, scenario_row),stringsAsFactors = FALSE)
    }
  }
  
  
  # Summarize results over seeds and create final results tables
  no_seed_table = param_table[,-c(which(colnames(param_table)=="SEED"))]
  no_seed_table = unique(no_seed_table)
  final_seed_table = merge(no_seed_table, processed_OM_sim, by = c("Scenario_Name"))
  
  # Write result tables (summarized and with seeds) to files
  write.table(final_seed_table, final_seed_table_dest, sep="\t",
              col.names = TRUE, row.names = FALSE, quote=FALSE)
  
  if (results_folder =="/scicore/home/penny/GROUP/smc_lai/E5_2_CT_LAI/postprocessing/"){
  aggregated_OM =   processed_OM_sim %>% group_by(Scenario_Name) %>% summarise_at(c(names(processed_OM_sim)[which(names(processed_OM_sim)=="prev_beg"):length(names(processed_OM_sim) ) ]),median,na.rm=TRUE)
  write.table(final_table, final_table_dest, sep="\t",
              col.names = TRUE, row.names = FALSE, quote=FALSE)
  
  final_table = merge(no_seed_table, aggregated_OM, by = c("Scenario_Name"))
  }

}

# Wrapper for looping across all simulation results and gathering postprocessing results in a table
# Version to use during adaptive sampling, returns the postprocessing results instead 
# of writing them to a file
postprocess_OM_as = function(results_folder, param_table, follow_up) {
  processed_OM_sim = NULL
  for( i in 1:nrow(param_table)) {
    # Read the OM simulation result
    OM_result_file = paste(results_folder, param_table[i,]$Scenario_Name, "_", 
                           param_table[i,]$SEED, "_out.txt", sep="")
    # Calculate the necessary outputs
    if(file.exists(OM_result_file) & file.info(OM_result_file)$size > 0) {
      OM_result = read.table(OM_result_file, sep="\t")
      # om_result, scenario_params, total_pop, survey_start, survey_end, int_start, int_end, pulsed_int_start
      scenario_row = calculate_outputs(OM_result, param_table[i,], follow_up)
      processed_OM_sim = rbind(processed_OM_sim, scenario_row)
    }
  }
  
  # # Summarize results over seeds and create final results tables
  aggregated_OM =   processed_OM_sim %>% group_by(Scenario_Name) %>% summarise_at(c(names(processed_OM_sim)[which(names(processed_OM_sim)=="sum_clin"):length(names(processed_OM_sim) ) ]),median,na.rm=TRUE)
  no_seed_table = param_table[,-c(which(colnames(param_table)=="SEED"))]
  no_seed_table = unique(no_seed_table)
  final_seed_table = merge(no_seed_table, processed_OM_sim, by = c("Scenario_Name"))
  # final_table = merge(no_seed_table, aggregated_OM, by = c("Scenario_Name"))
  
  return(final_seed_table)
}

# # # # Only for testing:
#om_results_folder = "/scicore/home/penny/GROUP/smc_lai/E4_SMCLAIdisc/om/"
#split_file = "/scicore/home/penny/GROUP/smc_lai/E4_SMCLAIdisc/postprocessing_5/seeds_E4SMCLAIdisc_Mali_4.9167_hill_0.1_3.txt"
#dest_dir = "/scicore/home/penny/GROUP/smc_lai/E4_SMCLAIdisc/postprocessing_5/"
#follow_up = 5

# om_results_folder = "/scicore/home/penny/GROUP/smc_lai/E3_SMCSMCdisc/om/"
# split_file = "/scicore/home/penny/GROUP/smc_lai/E3_SMCSMCdisc/postprocessing_5/seeds_E3SMCSMCdisc_Mali_4.9167_hill_0.1_150.txt"
# dest_dir = "/scicore/home/penny/GROUP/smc_lai/E3_SMCSMCdisc/postprocessing_5/"
# follow_up = 5
# 
# 
# # Create output file names
# split_name = basename(split_file)
# dest_table_agg = paste(dest_dir, "agg_", split_name, sep="")
# dest_table_seeds = paste(dest_dir, "seeds_", split_name, sep="")
# 
# 
# results_folder = om_results_folder
# param_table_file = split_file
# final_table_dest = dest_table_agg
# final_seed_table_dest = dest_table_seeds
# param_table = read.table(param_table_file, sep= "\t", as.is = TRUE, header = TRUE, stringsAsFactors = FALSE)
# 
# #
# i=4501
#  OM_result_file = paste(results_folder, param_table[i,]$Scenario_Name, "_",
#                         +                                param_table[i,]$seed, "_out.txt", sep="")
#  OM_result = read.table(OM_result_file, sep="\t")
# 
#  om_result <- OM_result
#  scenario_params <- param_table[i,]
# # # # #
# # scenario_params <- param_table[i,]
# # # # # # # Postprocess the OpenMalaria simulations
#  res <- postprocess_OM(results_folder = om_results_folder, param_table_file = split_file,
#                        final_table_dest = dest_table_agg, final_seed_table_dest = dest_table_seeds,
#                        follow_up)
# 
# 
# # Only for trial testing:
# om_results_folder = "/scicore/home/penny/GROUP/smc_lai/E5_2_CT_LAI/om/"
# split_file = "/scicore/home/penny/GROUP/smc_lai/E5_2_CT_LAI/postprocessing/split/E5_2_CT_LAI_Mali_4.9167__exp_47.txt"
# dest_dir = "/scicore/home/penny/GROUP/smc_lai/E5_2_CT_LAI/postprocessing/"
# follow_up = 1

#  om_results_folder = "/scicore/home/penny/GROUP/smc_lai/E5_1_CT_SMC/om/"
#  split_file = "/scicore/home/penny/GROUP/smc_lai/E5_1_CT_SMC/postprocessing/split/E5_1_CT_SMC_Sen_4.9167_31.295_exp_47.txt"
#  dest_dir = "/scicore/home/penny/GROUP/smc_lai/E5_1_CT_SMC/postprocessing/"
#  follow_up = 1
# 
# # Create output file names
#  split_name = basename(split_file)
#  dest_table_agg = paste(dest_dir, "agg_", split_name, sep="")
#  dest_table_seeds = paste(dest_dir, "seeds_", split_name, sep="")
# 
# #
#  results_folder = om_results_folder
#  param_table_file = split_file
#  final_table_dest = dest_table_agg
#  final_seed_table_dest = dest_table_seeds
#  param_table = read.table(param_table_file, sep= "\t", as.is = TRUE, header = TRUE, stringsAsFactors = FALSE)
# 
# #
#  i=3
#  OM_result_file = paste(results_folder, param_table[i,]$Scenario_Name, "_",
#                       +                                param_table[i,]$SEED, "_out.txt", sep="")
#  OM_result = read.table(OM_result_file, sep="\t")
# #
#  om_result <- OM_result
# #
#  scenario_params <- param_table[i,]
#  
#  if(param_table[i,"Seasonality"]=="Mali"){
#    timing_df = data.frame(EIR=c(5,9,20,47,150),
#                           timingInt=c(5,5,5,4,3),
#                           start= c(47,47,47,46,45),
#                           end=c(71,71,71,70,69))
#  } else {
#    timing_df = data.frame(EIR=c(5,9,20,47,150),
#                           timingInt=c(5,5,5,4,3),
#                           start= c(53,53,53,52,51),
#                           end=c(71,71,71,70,69))
#    
#  }
#  
#  start <- timing_df[which(param_table[i,"EIR"]==timing_df[,"EIR"]),"start"]
#  end <- timing_df[which(param_table[i,"EIR"]==timing_df[,"EIR"]),"end"]
 
# # # # # Postprocess the OpenMalaria simulations
# # # res <- postprocess_OM_trial(results_folder = om_results_folder, param_table_file = split_file,
# # #                       final_table_dest = dest_table_agg, final_seed_table_dest = dest_table_seeds,
# # #                       start,end)