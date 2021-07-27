########################################
# script generate_experiment_param_table.R
#
# construct a table with parameter values for each scenario of an OM experiment
#
# INPUTS:
#   table_file - name of the file where the table should be saved
#   experiment - name of the experiment
#
# OUTPUTS:
#	- a csv table with all the parameters for OM experiment

# created 21.01.2018
# monica.golumbeanu@unibas.ch
########################################
rm(list = ls())

library(tgp)

set.seed(42)

ranges_file = "~/smc_lai/param_ranges.RData"


# Seasonality and biting patterns values
seasons =read.table("~/smc_lai/analysis_workflow/resource_files/seasonality.txt", sep="\t", header = TRUE)
coverages <- data.frame(Coverage_SMC_1=c(1,1,0),Coverage_SMC_2=c(0.96,0.96,1),
                        Coverage_SMC_3=c(0.93,0.93,0.96),Coverage_SMC_4=c(0.90,0.90,0.93),Coverage_LAI_1=c(1,1,0),Coverage_LAI_2=c(0,0,1))

seasons <- cbind(seasons,coverages)
biting_pattern <- data.frame(Biting_pattern=c("Mali"),indoor=c(0.6),outdoor=c(0.4))
# for now: only low indoor biting and seasonal transmission 

seasons <- seasons[c(1,3),]
biting_patterns <- biting_pattern

# max age intervention

IntAge = data.frame(IntAge=c(4.9167),maxGroup=c(3))

# intervention decay

LAIdecay <- data.frame(fundecay=c("weibull","weibull","hill"),kdecay=c(1,0.69,8 ),Decay_Scen=c("exp","wei","hill" ) )#"weibull", 2

points_low <-c(3,4,8,28,150)
points_high <- c(5,9,20,47,150)
Access_df = data.frame(Access=c(rep(0.1,5),rep(0.5,5) ), EIR=c(points_low,points_high),timingInt=rep(c(25,25,25,20,15),2))

# Name of the experiment and parameters

Coverage = c(0.4, 1)
Halflife =c(30,150)
Efficacy= c(0.7,1)
Coverage_SMC = c(0.4,1)
param_ranges = rbind( Coverage,Coverage_SMC,Halflife,Efficacy)

Xcand = as.data.frame(lhs(1500, param_ranges))


colnames(Xcand) = c("Coverage","Coverage_SMC","Halflife","Efficacy")

SEED = c(1:5)


# Data frame with parameter ranges for the experiment
 param_ranges = rbind( c(0.4, 1),c(0.4, 1), c(30, 150), c(0.7, 1))
 row.names(param_ranges) = c( "Coverage","Coverage_SMC" ,"Halflife", "Efficacy")
save(param_ranges, file = ranges_file)

param_ranges_all <- param_ranges
# Table with the parameter values
param_tab_all = Reduce(merge, list(seasons, biting_patterns,IntAge,LAIdecay,Access_df, 
                               as.data.frame(Xcand)))

experiment = "E3_SMCSMCdisc"
dir.create(file.path(paste0("../../GROUP/smc_lai/",experiment,"/")), showWarnings = FALSE)

param_tab = cbind(scenarios_names= paste(experiment, 1:nrow(param_tab_all), sep="_"),
                  param_tab_all)
colnames(param_tab)[1] = "Scenario_Name"

# Add seed column at the beginning
param_tab = merge(param_tab, as.data.frame(SEED))

# remove unnescesary parameters


to_del <- which(names(param_tab)%in% c("Halflife", "Coverage","Efficacy"  ,"fundecay","kdecay"))
param_tab <- param_tab[,-c(to_del)]

table_file = paste0("../../GROUP/smc_lai/",experiment,"/param_tab.txt")

# Write table to specified destination file

file.remove(table_file)

param_tab <- param_tab[order(param_tab$IntAge,
                             param_tab$Access,
                             param_tab$Decay_Scen,
                             param_tab$Seasonality,
                             param_tab$EIR,
                             param_tab$Scenario_Name,
                             param_tab$SEED),]

write.table(param_tab, table_file, sep = "\t", quote = FALSE, col.names = TRUE,
            row.names = FALSE)
unique(param_tab[,c("Seasonality",'EIR','Decay_Scen',"Access")])

chunk <- 37500
n <- nrow(param_tab)
r  <- rep(1:ceiling(n/chunk),each=chunk)[1:n]
split_param_tab <- split(param_tab,r)

num_tab <- max(r)
for(j in 0:(num_tab-1)){ 
  
  table_file = paste0("../../GROUP/smc_lai/",experiment,"/param_tab_",j,".txt")
  
  # Write table to specified destination file
  file.remove(table_file)
  write.table(split_param_tab[[j+1]], table_file, sep = "\t", quote = FALSE, col.names = TRUE,
              row.names = FALSE)
  
}

# copy scaffold from git tracked folder 

file.remove(paste0("../../GROUP/smc_lai/",experiment,"/scaffold.xml"))

file.copy(paste0("~/smc_lai/",experiment,"/scaffold.xml"), paste0("../../GROUP/smc_lai/",experiment,"/scaffold.xml"),overwrite=TRUE)

# write ranges file to specified destination file
file.remove(paste0("../../GROUP/smc_lai/",experiment,"/param_ranges.RData"))
param_ranges = matrix(param_ranges_all[which(rownames(param_ranges_all)%in% c(names(param_tab))),], 
                         ncol = 2, 
                         nrow = length(which(rownames(param_ranges_all)%in% c(names(param_tab)))) )
row.names(param_ranges) = rownames(param_ranges_all)[which(rownames(param_ranges_all)%in% c(names(param_tab)))]

save(param_ranges, file = paste0("../../GROUP/smc_lai/",experiment,"/param_ranges.RData"))


# save parameters for Switch experiment

experiment = "E4_SMCLAIdisc"
dir.create(file.path(paste0("../../GROUP/smc_lai/",experiment,"/")), showWarnings = FALSE)

param_tab = cbind(scenarios_names= paste(experiment, 1:nrow(param_tab_all), sep="_"),
                  param_tab_all)
colnames(param_tab)[1] = "Scenario_Name"

# Add seed column at the beginning
param_tab = merge(param_tab, as.data.frame(SEED))
param_tab <- param_tab[order(param_tab$IntAge,
                             param_tab$Access,
                             param_tab$Decay_Scen,
                             param_tab$Seasonality,
                             param_tab$EIR,
                             param_tab$Scenario_Name,
                             param_tab$SEED),]

# remove unnescesary parameters


table_file = paste0("../../GROUP/smc_lai/",experiment,"/param_tab.txt")

# Write table to specified destination file
file.remove(table_file)
write.table(param_tab, table_file, sep = "\t", quote = FALSE, col.names = TRUE,
            row.names = FALSE)

chunk <- 37500
n <- nrow(param_tab)
r  <- rep(1:ceiling(n/chunk),each=chunk)[1:n]
split_param_tab <- split(param_tab,r)

num_tab <- max(r)
for(j in 0:(num_tab-1)){ 
  
  table_file = paste0("../../GROUP/smc_lai/",experiment,"/param_tab_",j,".txt")
  
  # Write table to specified destination file
  file.remove(table_file)
  write.table(split_param_tab[[j+1]], table_file, sep = "\t", quote = FALSE, col.names = TRUE,
              row.names = FALSE)
  
}

# copy scaffold from git tracked folder 


file.remove(paste0("../../GROUP/smc_lai/",experiment,"/scaffold.xml"))

file.copy(paste0("~/smc_lai/",experiment,"/scaffold.xml"), paste0("../../GROUP/smc_lai/",experiment,"/scaffold.xml"),overwrite=TRUE)


# write ranges file to specified destination file
param_ranges <- param_ranges
file.remove(paste0("../../GROUP/smc_lai/",experiment,"/param_ranges.RData"))
param_ranges = as.matrix(param_ranges_all[which(rownames(param_ranges_all)%in% c(names(param_tab))),])
row.names(param_ranges) = rownames(param_ranges_all)[which(rownames(param_ranges_all)%in% c(names(param_tab)))]

save(param_ranges, file = paste0("../../GROUP/smc_lai/",experiment,"/param_ranges.RData"))

