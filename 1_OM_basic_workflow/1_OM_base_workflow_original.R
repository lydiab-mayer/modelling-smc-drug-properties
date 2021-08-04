#####################################
#####################################
###                               ###
### STEP 1: OM SETUP & SIMULATION ###
###                               ###
#####################################
#####################################

##########################################################################
###
### M3TPP PROJECT:
### Main script for specifying parameter values 
### and running OpenMalaria simulations on the cluster 
### 
### Original script:
### Created 12.02.2021
### lydia.burgert@unibas.ch 
###
### Adapted script:
### 30.07.2021
### narimane.nekkab@unibas.ch
###
### R version 3.6.0
###
##########################################################################


##############
### HEADER ###
##############

# Clear environment
rm(list = ls())

# Set seed for replication
set.seed(42)

# Library
library(tgp)

# User 
user = strsplit(getwd(), "/", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]][5]

# Working directory
setwd(paste0("/scicore/home/penny/",user,"/M3TPP"))

# Source function scripts
source(paste0("./analysisworkflow/1_OM_basic_workflow/genOMsimscripts.R"))
source(paste0("./analysisworkflow/1_OM_basic_workflow/generate_param_table.R"))
source(paste0("./analysisworkflow/1_OM_basic_workflow/create_folders.R"))
source(paste0("./Experiments/iTPP1_TestCase/convert_access.R"))


##################
### EXPERIMENT ###
##################

# Insert experiment name here
exp ="..."

# Create folder in working directory
create_folders(exp) # <----- run 1st time then comment

# Optimize file size
chunk_size = 90000


##################
### PARAMETERS ###
##################

# Specify the desired parameter values 
# Note some values are categorical and others continuous
# Values are sampled for continuous variables from defined range

##################
### DEMOGRAPHY ###

# Population size
pop_size = 10000

# Demography dataframe
demography=data.frame(pop_size=pop_size)

##################
### MONITORING ###

# Years of analysis
start_year = 2030
end_year   = 2046

# Number of years burn in
burn_in_years = 80

# Set the burn in relative to start year
burn_in = start_year - burn_in_years

# Reference survey age group (in years) (default for now)
# 0 <= 0.25 = group 1
# 0.25 <= 2 = group 2
# 2 <= 5    = group 3
# 5 <= 10   = group 4
# 10 <= 15  = group 5
# 15 <= 20  = group 6
# 20 <= 100 = group 7

# Monitoring dataframe
monitoring=data.frame(start_year=start_year,
                      end_year=end_year,
                      burn_in=burn_in)

####################
### INTERVENTION ###

# Efficacy range (in %)
Efficacy= c(0.7,1)

# Half-life range (in days)
Halflife =c(30,150)

# Coverage range (in %) (single range, can have multiple)
Coverage = c(0.4, 1)

# Maximum age (refer to reference survey age group)
MaxAge = data.frame(MaxAge=c(10), 
                    maxGroup=c(4))

# Decay: define function, k and, efficacyB
# ---- function: "constant" or "step" or "linear" or "exponential" or "weibull" or "hill" or "smooth-compact"
# ---- k:         shape parameter of distribution. If not specified, default value of 1 is used.
# ---- efficacyB: measure of variation
Decay = data.frame(fundecay=c("hill"),
                   kdecay=c(8),
                   effB=c(1000))

#####################
### HEALTH SYSTEM ###

# Coverage of healthcare system: probability of healthcare seeking of uncomplicated malaria cases (in %) 
initial_access = data.frame(access=c(0.1,0.5)) 

# Convert access to care to 5 day probabilities for use in XML files
# Sources code from MMC project
access = pmax(convert_access(initial_access * 100), 0.04)

##################
### ENTOMOLOGY ###

# Seasonality type: monthly EIR vs. Fourier transformation (by month or with coefficients)
# Check XML structure to match
seasonality_type = "monthly"
# seasonality_type = "Fourier"

# Load seasonality file (path should be inside experiment)
if(seasonality_type == "monthly"){
  # seasonality = read.table(paste0("./Experiments/",exp,"/seasonality_4_months_updatedMay302021.txt"), sep="\t", header = TRUE)
  seasonality = read.table(paste0("./Experiments/",exp,"/exampleseason.txt"), sep="\t", header = TRUE)
}
if(seasonality_type == "Fourier"){
  seasonality = read.table(paste0("./Experiments/",exp,"/seasonality_Fourier_6_coeff.txt"), sep="\t", header = TRUE)
}

# Biting patterns 
biting_pattern <- data.frame(indoor=c(0.5),outdoor=c(0.5))

# EIR
EIR= data.frame(EIR=c(10,100))

###################
### DIAGNOSTICS ###

#############
### MODEL ###


##############################
### GENERATE PARAMS TABLES ###
##############################

# Categorical variables
param_cat = list(demography=demography,
                 monitoring=monitoring,
                 MaxAge=MaxAge,
                 Decay=Decay,
                 access=access,
                 seasonality=seasonality,
                 biting_pattern=biting_pattern,
                 EIR=EIR)

# Continuous variables
param_ranges_cont = rbind(Coverage,Halflife,Efficacy)

# Number of continuous parameters to sample via lhs
noSamples = 1

# Number of OM seeds per sample
noSeeds=  1

# Generate
gen_paramtable(exp, param_ranges_cont, param_cat, noSamples, noSeeds, chunk_size)


###########################################################
### GENERATE SCENARIOS AND RUN OPEN MALARIA SIMULATIONS ###
###########################################################

# Number of outputs to get in OM folder
2*nrow(demography)*nrow(monitoring)*nrow(MaxAge)*nrow(LAIdecay)*nrow(access)*nrow(seasonality)*nrow(biting_pattern)*nrow(EIR)*noSamples*noSeeds

# Run
genOMsimscripts(exp, chunk_size)



