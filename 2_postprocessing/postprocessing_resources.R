##########################
# Auxiliary functions for postprocessing OpenMalaria results
#
# created 02.10.2019
# monica.golumbeanu@unibas.ch
# modified by lydia.burgert@unibas.ch 
##########################

# library(rapportools)
# library(survival)
# library(cmprsk)
# require(plyr)
require(dplyr)

# Function which calculates prevalence reduction given an OpenMalaria
# simulation result.
calculate_outputs = function(om_result, scenario_params, follow_up, years_before_interv, cont) {
  colnames(om_result) = c("time", "age_group", "measure", "value")
  year_to_5day = 73
  month_to_5day = 6
  years_before_interv = years_before_interv
  
  # define age groups 
  age_groups <- c(0,0.25,2,5,10,15,20,100)
  minIntAge=0.25
  age210 = seq(which(age_groups==2),which(age_groups==10)-1)
  ageint = seq(which(age_groups==minIntAge),as.numeric(scenario_params["maxGroup"]))
  age05 = seq(which(age_groups==0),which(age_groups==5)-1)
  
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
  pop_int <- total_pop_age[total_pop_age$age_group %in% ageint,]%>% group_by( year) %>% summarise(n = sum(n) )
 
  pop_210 <- total_pop_age[total_pop_age$age_group %in% age210,]%>% group_by( year) %>% summarise(n = sum(n) )
  
  pop_05 <- total_pop_age[total_pop_age$age_group %in% age05,]%>% group_by( year) %>% summarise(n = sum(n) )
  
   # mean intervention population size 
  meanpopint <- mean(pop_int$n)
  
  # Summarize all measures results by summing up over age groups
  agg_om_result_total = om_result[,-which(names(om_result)=="age_group")]%>% group_by(time, measure) %>% summarise(value = sum(value))
  
  # Extract the clinical case numbers 
  trialpop <- sum( total_pop_age[total_pop_age$age_group %in% ageint & total_pop_age$year == years_before_interv +follow_up,"n"])
  
  ######################################
  #calculate output prevalence reduction
  ######################################   
  
  # Calculate the prevalence for all the monitored years 
  # Prevalence = total number of infected people (in age group)/ total population (in age group)
  
  # number of infected per age-group per time-step
  # n_infected = as.data.frame(om_result[om_result$measure == 1,]) # The number of human hosts with an infection (patent or not) on the reporting timestep
  n_infected = as.data.frame(om_result[om_result$measure == 3,]) # The number of human hosts whose total (blood-stage) parasite density is above the detection threshold
  n_infected$year <- rep(seq(1,max(n_infected$time)/year_to_5day),each=year_to_5day*max(om_result$age_group) )
  
  # total number of infected per time-step
  # n_infected_total = as.data.frame(agg_om_result_total[agg_om_result_total$measure == 1,]) # The number of human hosts with an infection (patent or not) on the reporting timestep
  n_infected_total = as.data.frame(agg_om_result_total[agg_om_result_total$measure == 3,]) # The number of human hosts whose total (blood-stage) parasite density is above the detection threshold
  n_infected_total$year = rep(seq(1,max(n_infected_total$time)/year_to_5day),each=year_to_5day )
  
  # calculate prevalence 
  
  # add the population size in respective age-groups
  n_infected <- inner_join(total_pop_age,n_infected, by=c("age_group","year"))
  n_infected_total <- inner_join(total_pop_all,n_infected_total ,by=c("year"))
  
  # divide number of infected  by respective age-group
  prev_agegroups = n_infected %>% mutate(prev = value/n )
  prev_allages = n_infected_total %>% mutate(prev = value/n )
  
  prev_int = n_infected[ n_infected$age_group %in% ageint, ] %>% group_by(time) %>%
    summarise(intinf = sum(value), intn=sum(n),year=mean(year))%>% mutate(prev = intinf/intn )
  
  prev_210 = n_infected[ n_infected$age_group %in% age210,] %>% group_by(time) %>%
    summarise(intinf = sum(value), intn=sum(n),year=mean(year))%>% mutate(prev = intinf/intn ) 
  
  # yearly average prevalence
  prev_agegroups_yearly <- prev_agegroups %>% group_by(age_group, year) %>% summarise(avg = mean(prev) )
  prev_allages_yearly = prev_allages %>% group_by( year) %>% summarise(avg = mean(prev) ) %>% ungroup()

  prev_all_before=as.numeric(prev_allages_yearly[which(prev_allages_yearly$year==years_before_interv),"avg"])
  prev_all_followup=as.numeric(prev_allages_yearly[which(prev_allages_yearly$year==(years_before_interv+follow_up)),"avg"])
  prev_red_all = ( prev_all_before - prev_all_followup )/ prev_all_before* 100
  prev_red_all = prev_red_all*(prev_red_all >= 0)
  
  # calculate the number of infected at follow up  in intervention age groups 
  inf_int_followup <- data.frame(n_infected[n_infected$year==years_before_interv+follow_up & n_infected$age_group %in% ageint,] %>% 
                                     group_by(time) %>%
                                     summarise(totcas = sum(value)) )
                                  
  # calculate yearly average prevalence at follow up for intervention age groups 
  prev_int_followup <- mean(inf_int_followup[ , "totcas"] /  pop_int[pop_int$year ==years_before_interv+follow_up,"n"]$n )
  
 
  # calculate the number of infected  before and yearly average prevalence interventions 
  
  inf_int_beg <- data.frame(n_infected[n_infected$year==years_before_interv & n_infected$age_group %in% ageint,] %>% 
                              group_by(time) %>%
                              summarise(totcas = sum(value)) )
  
  prev_int_beg <- mean(inf_int_beg[  , "totcas"] /  pop_int[pop_int$year ==years_before_interv+follow_up,"n"]$n )
  
  # prevalence reduction per age group
  prev_red_int = (prev_int_beg - prev_int_followup)/prev_int_beg * 100
  prev_red_int = prev_red_int*(prev_red_int >= 0) 
  
  
  # calculate the number of infected in age_group 2-10 and yearly average prevalence before interventions 
  
  inf_210_beg <- data.frame(n_infected[n_infected$year==years_before_interv & n_infected$age_group %in% age210,]%>% 
                                group_by(time)%>%
                            summarise(totcas = sum(value)) )
  
  prev_210_beg <- mean(inf_210_beg[  , "totcas"] /  pop_210[pop_210$year ==years_before_interv,"n"]$n )
  
  # calculate the number of infected  and yearly average prevalence in age_group 2-10 at followup  
  
  inf_210_followup <- data.frame(n_infected[n_infected$year==years_before_interv+follow_up & n_infected$age_group %in% age210,]%>% 
                              group_by(time)%>%
                              summarise(totcas = sum(value)) )
  
  prev_210_followup <- mean(inf_210_followup[  , "totcas"] /  pop_210[pop_210$year ==years_before_interv+follow_up,"n"]$n )
  
  # prevalence reduction per age group
  prev_red_210 = (prev_210_beg - prev_210_followup)/prev_210_beg * 100
  prev_red_210 = prev_red_210*(prev_red_210 >= 0) 
  
  ######################################
  #calculate output clinical cases
  ######################################   
  
  # Extract the clinical case numbers 
  
  # number of clinical cases
  nclin = as.data.frame(om_result[om_result$measure == 14,])
  nclin$year <- c(rep(seq(1,max(nclin$time)/year_to_5day),each=year_to_5day*max(om_result$age_group) ))
  
  # yearly average clinical cases per age-group
  inc_agegroups_yearly <- nclin %>% group_by(age_group ,year) %>% summarise(sum = sum(value) )
  
  # summed yearly clinical cases in intervention age-groups
  inc_all_yearly <- inc_agegroups_yearly%>% group_by( year) %>% summarise(inc = sum(sum) )
  inc_all_yearly$cpp <- inc_all_yearly$inc / tail(total_pop_all$n,1)
  
  # incidence reduction in intervention age group
  inc_red_all = (as.numeric(inc_all_yearly[inc_all_yearly$year==years_before_interv,"cpp"]) - as.numeric(inc_all_yearly[inc_all_yearly$year==(follow_up+years_before_interv),"cpp"]))/as.numeric(inc_all_yearly[inc_all_yearly$year==years_before_interv,"cpp"]) * 100
  inc_red_all = inc_red_all*(inc_red_all >= 0)
  
  # summed yearly clinical cases in intervention age-groups
  inc_int_yearly <- inc_agegroups_yearly[inc_agegroups_yearly$age_group  %in% ageint,]%>% group_by( year) %>% summarise(inc = sum(sum) )
  inc_int_yearly$cpp <- inc_int_yearly$inc / pop_int$n
  
  # incidence reduction in intervention age group
  inc_red_int = (as.numeric(inc_int_yearly[inc_int_yearly$year==years_before_interv,"cpp"]) - as.numeric(inc_int_yearly[inc_int_yearly$year==(follow_up+years_before_interv),"cpp"]))/as.numeric(inc_int_yearly[inc_int_yearly$year==years_before_interv,"cpp"]) * 100
  inc_red_int = inc_red_int*(inc_red_int >= 0) 
  
  # summed yearly clinical cases in  age-group 0-5
  inc_05_yearly <- inc_agegroups_yearly[inc_agegroups_yearly$age_group  %in% age05,]%>% group_by( year) %>% summarise(inc = sum(sum) )
  inc_05_yearly$cpp <- inc_05_yearly$inc / pop_05$n

  # incidence reduction in intervention age group
  inc_red_05 = (as.numeric(inc_05_yearly[inc_05_yearly$year==years_before_interv,"cpp"]) - as.numeric(inc_05_yearly[inc_05_yearly$year==(follow_up+years_before_interv),"cpp"]))/as.numeric(inc_05_yearly[inc_05_yearly$year==years_before_interv,"cpp"]) * 100
  inc_red_05 = inc_red_05*(inc_red_05 >= 0)
  
  # incidence over time 
  nclin <- inner_join(total_pop_age,nclin, by=c("age_group","year"))
  inc_agegroups = nclin %>% mutate(cpp = value/n )
  
  inc_allages = nclin %>% group_by(time) %>%
    summarise(intinc = sum(value), intn=sum(n),year=mean(year))%>% mutate(cpp = intinc/intn )%>% group_by(year)  %>% mutate(cppcum = cumsum(cpp))
  
  inc_int = nclin[ nclin$age_group %in% ageint, ] %>% group_by(time) %>%
    summarise(intinc = sum(value), intn=sum(n),year=mean(year))%>% mutate(cpp = intinc/intn )%>% group_by(year)  %>% mutate(cppcum = cumsum(cpp))
  
  inc_05 = nclin[ nclin$age_group %in% age05, ] %>% group_by(time) %>%
    summarise(intinc = sum(value), intn=sum(n),year=mean(year))%>% mutate(cpp = intinc/intn )%>% group_by(year)  %>% mutate(cppcum = cumsum(cpp))
  
  ######################################
  #calculate 5 months post intervention
  ######################################   
  
  ### this here needs to be implemented in the calculate_outputs function in the postprocesing_resources 
  inc_int_before <- subset(inc_int, year==years_before_interv)
  inc_int_before$timeyear <- inc_int_before$time-min(inc_int_before$time)
  
  inc_int_before_5mo = inc_int_before[ inc_int_before$timeyear %in% seq(18,48), ] 
  inc_int_before_5mo$cppcumnew <-  inc_int_before_5mo$cppcum - min(inc_int_before_5mo$cppcum) 
  
  inc_int_before_5mo_cont <-  max(inc_int_before_5mo$cppcumnew)
  
  inc_int_after <- subset(inc_int, year==years_before_interv+follow_up)
  inc_int_after$timeyear <- inc_int_after$time-min(inc_int_after$time)
  
  inc_int_after_5mo = inc_int_after[ inc_int_after$timeyear %in% seq(18,48), ] 
  inc_int_after_5mo$cppcumnew <-  inc_int_after_5mo$cppcum - min(inc_int_after_5mo$cppcum) 
  
  inc_int_after_5mo_int <-  max(inc_int_after_5mo$cppcumnew)
  
  # implement the inc_red_int_5mo in the return_row so it is also visible in the post-processing 
  inc_red_int_5mo <- ((inc_int_before_5mo_cont-inc_int_after_5mo_int)/inc_int_before_5mo_cont) *100
  
  ######################################
  # return results
  ######################################
  
   if(cont==FALSE) {
    # Final row with outputs to return
  return_row = cbind.data.frame(scenario_params$Scenario_Name,
                                scenario_params$SEED,prev_red_all,prev_red_210,prev_red_int,
                                inc_red_05,inc_red_int, inc_red_all,inc_red_int_5mo

  )
  colnames(return_row) = c("Scenario_Name",
                           "seed","prev_red_all","prev_red_210","prev_red_int", "inc_red_05","inc_red_int", "inc_red_all","inc_red_int_5mo"

  )

  return(return_row)}else{
    out_df= list("prevalence_210"=prev_210,
                 "prevalence_int"=prev_int,
                 "prevalence_allages"=prev_allages,
                 "prevalence_agegroups"=prev_agegroups,
                 "incidence_05"=inc_05,
                 "incidence_int"=inc_int,
                 "incidence_allages"=inc_allages,
                 "incidence_agegroups"=inc_agegroups,
                 "incidence_int_5mo"=inc_int
                 )
    return(out_df)
  }
  if(cont==FALSE) { 
    # Final row with outputs to return
    return_row = cbind.data.frame(scenario_params$Scenario_Name, 
                                  scenario_params$SEED,
                                  prev_red_210,
                                  prev_red_int,
                                  inc_red_int,
                                  inc_red_int_5mo
                                  
    ) 
    colnames(return_row) = c("Scenario_Name",
                             "seed",
                             "prev_red_210",
                             "prev_red_int", 
                             "inc_red_int", 
                             "inc_red_int_5mo"
                             
    )
    
    return(return_row)}else{
      out_df= list("prevalence_210"=prev_210,
                   "prevalence_int"=prev_int,
                   "incidence_int"=inc_int,
                   "incidence_int_5mo"=inc_int
      )
      return(out_df)
    }
}

# function which calculates outputs for a clinical trial run in OpenMalaria
calculate_trial_outputs = function(om_result, scenario_params,start,end,follow_up, years_before_interv) {
  colnames(om_result) = c("time", "age_group", "measure", "value")
  year_to_5day = 73
  month_to_5day = 6

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
  prev_total = n_infected[n_infected$age_group %in% pop, ]  %>% mutate(prev = value/n )
  
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
  nclin = as.data.frame(om_result[om_result$measure == 14,]) # number of episodes (uncomplicated) An episode of uncomplicated malaria is a period during which an individual has symptoms caused by malaria parasites present at the time of illness, where the symptoms do not qualifying as severe malari
  nclin$year <- c(rep(seq(1,max(nclin$time)/year_to_5day),each=year_to_5day*max(om_result$age_group) ))
  
  # yearly average clinical cases per age-group
  nclin_yearly <- nclin %>% group_by(age_group ,year) %>% summarise(sum = sum(value) )
  
  # summed yearly clinical cases in intervention age-groups
  final_nclin_intGroup <- nclin_yearly[nclin_yearly$age_group <=3,]
  
  # sum of intervention age-groups
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

  # Prevalence
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
                          final_seed_table_dest, follow_up,years_before_interv) {
  
  param_table = read.table(param_table_file, sep= "\t", as.is = TRUE, header = TRUE, stringsAsFactors = FALSE)
  processed_OM_sim = NULL
  for( i in 1:nrow(param_table)) {
    
    print(i)
    # Read the OM simulation result
    OM_result_file = paste(results_folder, param_table[i,]$Scenario_Name, "_", 
                           param_table[i,]$SEED, "_out.txt", sep="")
    # Calculate the necessary outputs
    if(file.exists(OM_result_file) & file.info(OM_result_file)$size > 0) {
      
      # Read in file
      OM_result = read.table(OM_result_file, sep="\t")
      
      # Identify error to skip
      mtry <- try(calculate_outputs(OM_result, param_table[i,], follow_up,years_before_interv,cont=FALSE),
                  silent = TRUE)
      
      # Skip error or calculate outputs
      if (class(mtry) != "try-error") {
        scenario_row = calculate_outputs(OM_result, param_table[i,], follow_up,years_before_interv,cont=FALSE)
        processed_OM_sim = data.frame(rbind(processed_OM_sim, scenario_row),stringsAsFactors = FALSE)
      } else {
        message("Error")
      }
    }
  }
  
  
  # Summarize results over seeds and create final results tables
  
  
  
  aggregated_OM =   processed_OM_sim %>% group_by(Scenario_Name) %>% summarise_at(c(names(processed_OM_sim)[which(names(processed_OM_sim)=="prev_red_210"):length(names(processed_OM_sim) ) ]),median,na.rm=TRUE)
  
  
  
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


############################################################ BY MONTH

# Function which calculates prevalence reduction given an OpenMalaria
# simulation result.
calculate_outputs_month = function(om_result, scenario_params, follow_up,years_before_interv,cont) {
  colnames(om_result) = c("time", "age_group", "measure", "value")
  year_to_5day = 73
  month_to_5day = 6
  years_before_interv = 5
  
  # define age groups 
  
  age_groups <- c(0,0.25,2,5,10,15,20,100)
  minIntAge=0.25
  
  age210 = seq(which(age_groups==2),which(age_groups==10)-1)
  ageint = seq(which(age_groups==minIntAge),as.numeric(scenario_params["maxGroup"]))
  age05 = seq(which(age_groups==0),which(age_groups==5)-1)
  # Remove first measurement as it includes all cases until then
  to_remove = which(om_result$time == 1)
  om_result = om_result[-to_remove,]
  
  # population per time step per age group
  total_pop = as.data.frame(om_result[om_result$measure == 0, ])
  
  # add years tot time-steps, 73 time-steps are one year 
  total_pop$year <- rep(seq(1,max(total_pop$time)/year_to_5day),each=year_to_5day*max(om_result$age_group) )
  # add months tot time-steps, 6 time-steps are one month 
  total_pop$month <- rep(seq(1,max(total_pop$time)/month_to_5day),each=month_to_5day*max(om_result$age_group) )
  
  # calculate total population per age group per month
  total_pop_age = data.frame(total_pop %>% group_by(age_group, month) %>% summarise(n = mean(value) ))
  
  # sum up total population per month 
  total_pop_all = total_pop %>% group_by( time,month) %>% summarise(sum = sum(value) )%>% group_by( month)  %>% summarise(n = mean(sum) )
  
  # sum up population intervention age-groups over the months 
  pop_int <- total_pop_age[total_pop_age$age_group %in% ageint,]%>% group_by( month) %>% summarise(n = sum(n) )
  
  pop_210 <- total_pop_age[total_pop_age$age_group %in% age210,]%>% group_by( month) %>% summarise(n = sum(n) )
  
  pop_05 <- total_pop_age[total_pop_age$age_group %in% age05,]%>% group_by( month) %>% summarise(n = sum(n) )
  
  # mean intervention population size 
  meanpopint <- mean(pop_int$n)
  
  
  # Summarize all measures results by summing up over age groups
  agg_om_result_total = om_result[,-which(names(om_result)=="age_group")]%>% group_by(time, measure) %>% summarise(value = sum(value))
  
  
  # Extract the clinical case numbers 
  trialpop <- sum( total_pop_age[total_pop_age$age_group %in% ageint & total_pop_age$month == years_before_interv*12 +follow_up,"n"])
  
  ######################################
  #calculate output prevalence reduction
  ######################################   
  
  # Calculate the prevalence for all the monitored months 
  # Prevalence = total number of infected people (in age group)/ total population (in age group)
  
  # number of infected per age-group per time-step
  n_infected = as.data.frame(om_result[om_result$measure == 1,])
  n_infected$year <- rep(seq(1,max(n_infected$time)/year_to_5day),each=year_to_5day*max(om_result$age_group) )
  n_infected$month <- rep(seq(1,max(n_infected$time)/month_to_5day),each=month_to_5day*max(om_result$age_group) )
  
  # total number of infected per time-step
  n_infected_total = as.data.frame(agg_om_result_total[agg_om_result_total$measure == 1,])
  n_infected_total$year = rep(seq(1,max(n_infected_total$time)/year_to_5day),each=year_to_5day )
  n_infected_total$month = rep(seq(1,max(n_infected_total$time)/month_to_5day),each=month_to_5day )
  
  # calculate prevalence 
  
  # add the population size in respective age-groups
  n_infected <- inner_join(total_pop_age,n_infected, by=c("age_group","month"))
  n_infected_total <- inner_join(total_pop_all,n_infected_total ,by=c("month"))
  
  
  # divide number of infected  by respective age-group
  prev_agegroups = n_infected %>% mutate(prev = value/n )
  prev_allages = n_infected_total %>% mutate(prev = value/n )
  
  # Prevalence
  prev_int = n_infected[ n_infected$age_group %in% ageint, ] %>% group_by(time) %>%
    summarise(intinf = sum(value), intn=sum(n),month=mean(month))%>% mutate(prev = intinf/intn )
  
  prev_210 = n_infected[ n_infected$age_group %in% age210,] %>% group_by(time) %>%
    summarise(intinf = sum(value), intn=sum(n),month=mean(month))%>% mutate(prev = intinf/intn ) 
  
  
  # monthly average prevalence
  prev_agegroups_monthly <- prev_agegroups %>% group_by(age_group, month) %>% summarise(avg = mean(prev) )
  prev_allages_monthly = prev_allages %>% group_by( month) %>% summarise(avg = mean(prev) ) %>% ungroup()
  
  prev_all_before=as.numeric(prev_allages_monthly[which(prev_allages_monthly$month==years_before_interv*12),"avg"])
  prev_all_followup=as.numeric(prev_allages_monthly[which(prev_allages_monthly$month==(years_before_interv*12+follow_up)),"avg"])
  prev_red_all = ( prev_all_before - prev_all_followup )/ prev_all_before* 100
  prev_red_all = prev_red_all*(prev_red_all >= 0) 
  
  
  # calculate the number of infected at follow up  in intervention age groups 
  inf_int_followup <- data.frame(n_infected[n_infected$month==years_before_interv*12+follow_up & n_infected$age_group %in% ageint,] %>% 
                                   group_by(time) %>%
                                   summarise(totcas = sum(value)) )
  
  # calculate monthly average prevalence at follow up for intervention age groups 
  prev_int_followup <- mean(inf_int_followup[ , "totcas"] /  pop_int[pop_int$month ==years_before_interv*12+follow_up,"n"]$n )
  
  
  # calculate the number of infected  before and monthly average prevalence interventions 
  
  inf_int_beg <- data.frame(n_infected[n_infected$month==years_before_interv*12 & n_infected$age_group %in% ageint,] %>% 
                              group_by(time) %>%
                              summarise(totcas = sum(value)) )
  
  prev_int_beg <- mean(inf_int_beg[  , "totcas"] /  pop_int[pop_int$month ==years_before_interv*12+follow_up,"n"]$n )
  
  # prevalence reduction per age group
  prev_red_int = (prev_int_beg - prev_int_followup)/prev_int_beg * 100
  prev_red_int = prev_red_int*(prev_red_int >= 0) 
  
  
  # calculate the number of infected in age_group 2-10 and monthly average prevalence before interventions 
  
  inf_210_beg <- data.frame(n_infected[n_infected$month==years_before_interv*12 & n_infected$age_group %in% age210,]%>% 
                              group_by(time)%>%
                              summarise(totcas = sum(value)) )
  
  prev_210_beg <- mean(inf_210_beg[  , "totcas"] /  pop_210[pop_210$month ==years_before_interv*12,"n"]$n )
  
  # calculate the number of infected  and monthly average prevalence in age_group 2-10 at followup  
  
  inf_210_followup <- data.frame(n_infected[n_infected$month==years_before_interv*12+follow_up & n_infected$age_group %in% age210,]%>% 
                                   group_by(time)%>%
                                   summarise(totcas = sum(value)) )
  
  prev_210_followup <- mean(inf_210_followup[  , "totcas"] /  pop_210[pop_210$month ==years_before_interv*12+follow_up,"n"]$n )
  
  # prevalence reduction per age group
  prev_red_210 = (prev_210_beg - prev_210_followup)/prev_210_beg * 100
  prev_red_210 = prev_red_210*(prev_red_210 >= 0) 
  
  
  
  ######################################
  #calculate output clinical cases
  ######################################   
  
  # Extract the clinical case numbers 
  
  # number of clinical cases
  nclin = as.data.frame(om_result[om_result$measure == 14,])
  nclin$year <- c(rep(seq(1,max(nclin$time)/year_to_5day),each=year_to_5day*max(om_result$age_group) ))
  nclin$month <- c(rep(seq(1,max(nclin$time)/month_to_5day),each=month_to_5day*max(om_result$age_group) ))
  
  
  # monthly average clinical cases per age-group
  inc_agegroups_monthly <- nclin %>% group_by(age_group ,month) %>% summarise(sum = sum(value) )
  
  # summed monthly clinical cases in intervention age-groups
  inc_all_monthly <- inc_agegroups_monthly%>% group_by( month) %>% summarise(inc = sum(sum) )
  inc_all_monthly$cpp <- inc_all_monthly$inc / tail(total_pop_all$n,1)
  
  # incidence reduction in intervention age group
  inc_red_all = (as.numeric(inc_all_monthly[inc_all_monthly$month==years_before_interv*12,"cpp"]) - as.numeric(inc_all_monthly[inc_all_monthly$month==(follow_up+years_before_interv*12),"cpp"]))/as.numeric(inc_all_monthly[inc_all_monthly$month==years_before_interv*12,"cpp"]) * 100
  inc_red_all = inc_red_all*(inc_red_all >= 0) 
  
  
  # summed monthly clinical cases in intervention age-groups
  inc_int_monthly <- inc_agegroups_monthly[inc_agegroups_monthly$age_group  %in% ageint,]%>% group_by( month) %>% summarise(inc = sum(sum) )
  inc_int_monthly$cpp <- inc_int_monthly$inc / pop_int$n
  
  # incidence reduction in intervention age group
  inc_red_int = (as.numeric(inc_int_monthly[inc_int_monthly$month==years_before_interv*12,"cpp"]) - as.numeric(inc_int_monthly[inc_int_monthly$month==(follow_up+years_before_interv*12),"cpp"]))/as.numeric(inc_int_monthly[inc_int_monthly$month==years_before_interv*12,"cpp"]) * 100
  inc_red_int = inc_red_int*(inc_red_int >= 0) 
  
  
  # summed monthly clinical cases in  age-group 0-5 
  inc_05_monthly <- inc_agegroups_monthly[inc_agegroups_monthly$age_group  %in% age05,]%>% group_by( month) %>% summarise(inc = sum(sum) )
  inc_05_monthly$cpp <- inc_05_monthly$inc / pop_05$n
  
  # incidence reduction in intervention age group
  inc_red_05 = (as.numeric(inc_05_monthly[inc_05_monthly$month==years_before_interv*12,"cpp"]) - as.numeric(inc_05_monthly[inc_05_monthly$month==(follow_up+years_before_interv*12),"cpp"]))/as.numeric(inc_05_monthly[inc_05_monthly$month==years_before_interv*12,"cpp"]) * 100
  inc_red_05 = inc_red_05*(inc_red_05 >= 0) 
  
  
  # incidence over time 
  nclin <- inner_join(total_pop_age,nclin, by=c("age_group","month"))
  inc_agegroups = nclin %>% mutate(cpp = value/n )
  
  inc_allages = nclin %>% group_by(time) %>%
    summarise(intinc = sum(value), intn=sum(n),month=mean(month))%>% mutate(cpp = intinc/intn )%>% group_by(month)  %>% mutate(cppcum = cumsum(cpp))
  
  inc_int = nclin[ nclin$age_group %in% ageint, ] %>% group_by(time) %>%
    summarise(intinc = sum(value), intn=sum(n),month=mean(month))%>% mutate(cpp = intinc/intn )%>% group_by(month)  %>% mutate(cppcum = cumsum(cpp))
  
  inc_05 = nclin[ nclin$age_group %in% age05, ] %>% group_by(time) %>%
    summarise(intinc = sum(value), intn=sum(n),month=mean(month))%>% mutate(cpp = intinc/intn )%>% group_by(month)  %>% mutate(cppcum = cumsum(cpp))
  
  
  
  if(cont==FALSE) { 
    # Final row with outputs to return
    return_row = cbind.data.frame(scenario_params$Scenario_Name, 
                                  scenario_params$SEED,prev_red_all,prev_red_210,prev_red_int,
                                  inc_red_05,inc_red_int, inc_red_all
                                  
    ) 
    colnames(return_row) = c("Scenario_Name",
                             "seed","prev_red_all","prev_red_210","prev_red_int", "inc_red_05","inc_red_int", "inc_red_all"
                             
    )
    
    return(return_row)}else{
      out_df= list("prevalence_210"=prev_210,
                   "prevalence_int"=prev_int,
                   "prevalence_allages"=prev_allages,
                   "prevalence_agegroups"=prev_agegroups,
                   "incidence_05"=inc_05,
                   "incidence_int"=inc_int,
                   "incidence_allages"=inc_allages,
                   "incidence_agegroups"=inc_agegroups
      )
      return(out_df)
    }
}


# Wrapper for looping across all simulation results and gathering postprocessing results in a table
postprocess_OM_month = function(results_folder, param_table_file, final_table_dest, 
                          final_seed_table_dest, follow_up,years_before_interv) {
  
  param_table = read.table(param_table_file, sep= "\t", as.is = TRUE, header = TRUE, stringsAsFactors = FALSE)
  processed_OM_sim = NULL
  for( i in 1:nrow(param_table)) {
    
    print(i)
    # Read the OM simulation result
    OM_result_file = paste(results_folder, param_table[i,]$Scenario_Name, "_", 
                           param_table[i,]$SEED, "_out.txt", sep="")
    # Calculate the necessary outputs
    if(file.exists(OM_result_file) & file.info(OM_result_file)$size > 0) {
      
      # Read in file
      OM_result = read.table(OM_result_file, sep="\t")
      
      # Identify error to skip
      mtry <- try(calculate_outputs_month(OM_result, param_table[i,], follow_up,years_before_interv,cont=FALSE),
                  silent = TRUE)
      
      # Skip error or calculate outputs
      if (class(mtry) != "try-error") {
        scenario_row = calculate_outputs_month(OM_result, param_table[i,], follow_up,years_before_interv,cont=FALSE)
        processed_OM_sim = data.frame(rbind(processed_OM_sim, scenario_row),stringsAsFactors = FALSE)
      } else {
        message("Error")
      }
    }
  }
  
  
  # Summarize results over seeds and create final results tables
  
  
  
  aggregated_OM =   processed_OM_sim %>% group_by(Scenario_Name) %>% summarise_at(c(names(processed_OM_sim)[which(names(processed_OM_sim)=="prev_red_all"):length(names(processed_OM_sim) ) ]),median,na.rm=TRUE)
  
  
  
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