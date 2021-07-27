######################################
###
### Decay shapes: functions & plots
###
######################################

##### Melissa comment on mock mAb example

# <PEV> <!--assuming csp (from RTS.S-A01), anti-circumsporozoit -->
#   <decay L=“@Halflife@d” k=“@kdecay@” function=“@fundecay@” />
#     <efficacyB value=“10.0" />
#           <initialEfficacy value=“@Efficacy@” />
#         </PEV>
# @fundecay = “hill”
# k=“@kdecay@” choose 8 to be strong sigmoidal
# efficacyB value=“10.0" SHOULD be 1000 so no inter-individual variation
# and vary L from 60 to 180

##### Description of PEV parameters

# <decay
# function=("constant" or "step" or "linear" or "exponential" or "weibull" or "hill" or "smooth-compact")
# [ L=string ]
# [ k=double ] DEFAULT VALUE 1.0
# [ CV=double ] DEFAULT VALUE 0
# />

# DecayFunction type or https://github.com/SwissTPH/openmalaria/wiki/ModelDecayFunctions

# function
# function=("constant" or "step" or "linear" or "exponential" or "weibull" or "hill" or "smooth-compact")
# Units: None Min: 0 Max: 1
# Determines which decay function to use. Available decay functions, for age t in years: 
# constant: 1 step: 1 for t less than L, otherwise 0 linear: 1 - t/L for t less than L, 
# otherwise 0 exponential: exp( - t/L * log(2) ) weibull: exp( -(t/L)^k * log(2) ) 
# hill: 1 / (1 + (t/L)^k) smooth-compact: exp( k - k / (1 - (t/L)^2) ) for t less than L, otherwise 0

# L
# L=string
# Units: User-defined (defaults to years) Min: 0
# (Time) scale parameter of distribution: this is either the age of complete decay (smooth-compact, 
# step and linear functions) or the age at which the parameter has decayed to half its original value 
# (exponential, weibull and hill). Not used when function="constant" (i.e. no decay). This value can 
# be specified in years, days or steps (e.g. 2y, 180d or 100t). When the unit is not specified years 
# are assumed. The value is used without rounding except when sampling an age of decay, when the rounding 
# happens as late as possible.

# k
# k=double
# Units: none Min: 0
# Default value: 1.0
# Shape parameter of distribution. If not specified, default value of 1 is used. 
# Meaning depends on function; not used in some cases.

# Coefficient of Variation
# CV=double
# Min: 0
# Default value: 0
# If CV is non-zero, heterogeneity of decay is introduced via a random variable sampled 
# from the log-normal distribution. This distribution is parameterised with mean=1 and CV as given. 
# The effective age of decay is the real age multiplied by this variable (for decay functions 
# with a half-life, this is equivalent to dividing the half-life by the variable).

# <efficacyB
# value=double
# />
#   
# Units: Positive real Min: 0.001 Max: 1.00E+06
# Measure of variation in vaccine efficacy: efficacy is sampled from a beta distribution 
# with efficacyB its beta parameter and its alpha parameter fixed such that the mean is 
# that given by initialEfficacy.
# Input parameter value
# value=double
# A double-precision floating-point value.

# <initialEfficacy
# value=double
# />
#   
# Units: dimensionless Min: 0 Max: 1
# Mean efficacy values before decay (see efficacyB and decay parameter descriptions for 
# sampling and decay). The i-th value in this list is used for the efficacy of the vaccine 
# after the i-th dose. Where more doses are given than there are values in this list, the 
# last value is repeated.
# Input parameter value
# value=double

########### FUNCTIONS

# exponential:
# exp( - t/L * log(2) )

# weibull:
# exp( -(t/L)^k * log(2) )

# hill:
# 1 / (1 + (t/L)^k) 

# smooth-compact:
# exp( k - k / (1 - (t/L)^2) )

########### PLOT EXAMPLE DISTRIBUTIONS (without heterogeneity)

# test exponential
L= 100
t_seq= seq(1,100)
exp( - t_seq/L * log(2) )
eff_init = 100
eff_exponential = eff_init*exp( - t_seq/L * log(2) )
plot(t_seq,eff_exponential)

# test weibull
L= 100
t_seq= seq(1,100)
k = 2
exp( -(t_seq/L)^k * log(2) )
eff_init = 100
eff_weibull = eff_init* exp( -(t_seq/L)^k * log(2) )
plot(t_seq,eff_weibull)

# test smooth-compact
L= 100
t_seq= seq(1,100)
k = 1
exp( k - k / (1 - (t_seq/L)^2) )
eff_init = 100
eff_smoothcompact = eff_init*exp( k - k / (1 - (t_seq/L)^2) )
plot(t_seq,eff_smoothcompact)

# test Hill function 
L=50
t_seq= seq(1,100)
k_hill = 2
eff_hill =eff_init*(1 / (1 + (t_seq/L)^k_hill))
plot(t_seq,eff_hill,ylim=c(0,100))

# test Hill function (strong sigmoidal b/c k=8)
L=50
t_seq= seq(1,100)
k_hill_sig = 8
eff_hill_sig =eff_init*(1 / (1 + (t_seq/L)^k_hill_sig))
plot(t_seq,eff_hill_sig,ylim=c(0,100))




