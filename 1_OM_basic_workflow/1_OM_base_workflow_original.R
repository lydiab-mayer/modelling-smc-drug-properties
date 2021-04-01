##############################
# Main script for specifying parameter values and running OpenMalaria simulations on the cluster. 
# 
#
# created 12.02.2021
# lydia.burgert@unibas.ch 
#############################

# Setup
rm(list = ls())
library(tgp)
set.seed(42)

# User 
user = strsplit(getwd(), "/", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]][5]

# Working directory
setwd(paste0("/scicore/home/penny/",user,"/M3TPP"))

# Source function scripts
source(paste0("./analysisworkflow/1_OM_basic_workflow/genOMsimscripts.R"))
source(paste0("./analysisworkflow/1_OM_basic_workflow/generate_param_table.R"))
source(paste0("./analysisworkflow/1_OM_basic_workflow/create_folders.R"))

# Insert experiment name here
exp ="..."

# Create folder in working directory
create_folders(exp) # <----- run 1st time then comment

# Seasonality type: monthly EIR vs. Fourier transformation (choose 1)
# seasonality_type = "monthly"
seasonality_type = "Fourier"

# !!! Depending on which seasonality is used, choose the correct scaffold.xml to upload
# !!! If using Fourier, make sure number of coefficients correspond (example here I use a 6 coeff model with requires 5 a and 5 b values)

# !!! Copy the scaffold.xml file for your experiment into the ./M3TPP/Experiments/"exp"/OM_jobs folder!!!
# !!! Copy the same scaffold.xml file for your experiment into the .GROUP//M3TPP/Experiments/"exp"/OM_jobs folder!!!
# !!! Copy a seasonality.csv file for your experiment into the ./M3TPP/Experiments/"exp"/ folder!!!

######################################
# Specify the desired parameter values 
######################################

chunk_size = 90000

#############################
# Categorical variables

# Seasonality
if(seasonality_type == "monthly"){
  Seasonality =read.table(paste0("./Experiments/",exp,"/seasonality_monthly.txt"), sep="\t", header = TRUE)
}
if(seasonality_type == "Fourier"){
  Seasonality =read.table(paste0("./Experiments/",exp,"/seasonality_Fourier_6_coeff.txt"), sep="\t", header = TRUE)
}

# Biting patterns 
Biting_pattern <- data.frame(Biting_pattern=c("Mali"),indoor=c(0.6),outdoor=c(0.4))

# EIR
EIR= data.frame(EIR=c(10))

# Max age intervention
MaxAge = data.frame(MaxAge=c(4.9167),maxGroup=c(3))

# Intervention decay
LAIdecay <- data.frame(fundecay=c("weibull"),kdecay=c(1 ),LAIdecay=c("exp" ) )

# Coverage of healthcare system
Access = data.frame(Access=c(0.1))  

# Combine
param_cat = list(Seasonality=Seasonality,
                 Biting_pattern=Biting_pattern,
                 EIR=EIR,
                 MaxAge=MaxAge,
                 LAIdecay=LAIdecay,
                 Access=Access)

#############################
# Continuous variables and their ranges

# Name of the experiment and parameters
param_ranges_cont = rbind( Coverage = c(0.4, 1),
                           Halflife =c(30,150),
                           Efficacy= c(0.7,1) )



###############################
# Generate the parameter tables  
###############################

# no. of continuous parameters sampled via lhs
noSamples = 5

# no. of OM seeds per sample
noSeeds=  2

# Generate
gen_paramtable(exp, param_ranges_cont,param_cat, noSamples, noSeeds,chunk_size)

########################################################
# Generate submission scripts and run the OM simulations
########################################################\

# Run
genOMsimscripts(exp, chunk_size)



