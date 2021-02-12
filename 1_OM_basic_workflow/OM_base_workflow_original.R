##############################
# Main script for specifying parameter values and running OpenMalaria simulations on the cluster. 
# 
#
# created 12.02.2021
#lydia.burgert@unibas.ch 
#############################
rm(list = ls())
setwd("~/M3TPP")
set.seed(42)

library(tgp)

source("./analysisworkflow/1_OM_basic_workflow/genOMsimscripts.R")
source("./analysisworkflow/1_OM_basic_workflow/generate_param_table.R")
source("./analysisworkflow/1_OM_basic_workflow/create_folders.R")

# insert experiment name here
exp ="..."

create_folders(exp)

chunk_size = 90000

# !!! copy the scaffold.xml file for your experiment into the ./M3TPP/Experiments/"exp"/OM_jobs folder!!!

#############################
# Specify the desired parameter values 
#############################

# categorical variables
#############################

# Seasonality and biting patterns values
seasons =read.table(paste0("./Experiments/",exp,"/seasonality.txt"), sep="\t", header = TRUE)

biting_pattern <- data.frame(Biting_pattern=c("Mali"),indoor=c(0.6),outdoor=c(0.4))

EIR= data.frame(EIR=c(10))

# max age intervention

IntAge = data.frame(IntAge=c(4.9167),maxGroup=c(3))

# intervention decay

LAIdecay <- data.frame(fundecay=c("weibull"),kdecay=c(1 ),Decay_Scen=c("exp" ) )

Access_df = data.frame(Access=c(0.1))

param_cat = list(seasons,
                 biting_pattern,
                 EIR,
                 IntAge,
                 LAIdecay,
                 Access_df)


# continuous variables and their ranges
#############################

# Name of the experiment and parameters
param_ranges_cont = rbind( Coverage = c(0.4, 1),
                      Halflife =c(30,150),
                      Efficacy= c(0.7,1) )




# generate the parameter tables  
#############################

# no. of continuous parameters sampled via lhs
noSamples = 5

# no. of OM seeds per sample
noSeeds=  2

gen_paramtable(exp, param_ranges_cont,param_cat, noSamples, noSeeds,chunk_size)
  
# generate submission scripts and run the OM simulations
#############################
genOMsimscripts(exp, chunk_size)


  
