##############################
# Main script for running OpenMalaria simulations on the cluster. 
# 
#	SIM_FOLDER = simulation  folder with all necessary input files:
#			param_tab.txt = file where each line is a parameter configuration to simulate
#			scaffold.xml = xml file containing @parameter@ wildcards for varied values across simulations
# OUTPUT:
#	The script creates the following folders in the specified SIM_FOLDER:
#		base/ = folder containing replacement patterns and base xml files
#		scenarios/ = folder with scenario files
#		om/ = folder with OpenMalaria simulations
#
# SYNTHAX: 
#	bash OM_base_workflow.sh SIM_FOLDER
# 
#
# created 03.05.2020
#lydia.burgert@unibas.ch adapted from theresa.reiker@unibas.ch
#############################
rm(list = ls())
setwd("~/M3TPP")
set.seed(42)

library(tgp)

source("./analysisworkflow/1_OM_basic_workflow/genOMsimscripts.R")
source("./analysisworkflow/1_OM_basic_workflow/generate_param_table.R")

exp ="E0_test"

chunk_size = 90000

#############################
# Specify the desired parameter values 
#############################

# categorical variables
#############################

# Seasonality and biting patterns values
seasons =read.table(paste0("./Experiments/",exp,"/seasonality.txt"), sep="\t", header = TRUE)

biting_pattern <- data.frame(Biting_pattern=c("Mali"),indoor=c(0.6),outdoor=c(0.4))

# max age intervention

IntAge = data.frame(IntAge=c(4.9167),maxGroup=c(3))#,9.9167 and ,4

# intervention decay

LAIdecay <- data.frame(fundecay=c("weibull"),kdecay=c(1 ),Decay_Scen=c("exp" ) )#"weibull", 2

Access_df = data.frame(Access=c(0.1))

param_cat = list(seasons,
                 biting_pattern,
                 IntAge,
                 LAIdecay,
                 Access_df)


# continuous variables and their ranges
#############################

# Name of the experiment and parameters
param_ranges_cont = rbind( EIR = c(1, 25),
                      Coverage = c(0.4, 1),
                      Halflife =c(30,150),
                      Efficacy= c(0.7,1) )


noSamples = 5

noSeeds=  2


gen_paramtable(exp, param_ranges_cont,param_cat, noSamples, noSeeds,chunk_size)
  


genOMsimscripts(exp, chunk_size)


  
