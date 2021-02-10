###########################
# Optimization 
# 
# created 06.12.2019
# monica.golumbeanu@unibas.ch
###########################
library(nloptr)
library(Rsolnp)
library(metaheuristicOpt)

args = commandArgs(TRUE)
gp_file = args[1]
ranges_file = args[2]
results_folder = args[3]

gp_file = "~/MMC/TPP/simulations/MAB_once_3years_avg_prev/gp_4/as/seeds_MAB_once_3_years_perennial_High_indoor_cv_as.RData"
# gp_file = "~/MMC/TPP/simulations/MAB_twice_3years/gp_4/as/seeds_MAB_twice_3_years_seasonal_High_indoor_cv_as.RData"
ranges_file = "~/MMC/TPP/simulations/MAB_once_3years_cont_EIR/param_ranges.RData"
results_folder = "~/MMC/TPP/simulations/MAB_once_3years_cont_EIR/gp_4/optimization/"

eval_f0 <- function(x, EIR, minPrevReduc) { 
    return(x)
}

eval_g0 <- function(x, EIR, minPrevReduc){
    #fixed_effi0.85_hl0.1667"
    x_input = t(c(EIR, 0.1667, 0.85, x, 0.04))
    # print(x_input)
    prediction <- predict(x = x_input, object = gp_result$GP_model)
    diffBetweenMinReductionRequiredAndPred = minPrevReduc - prediction$mean
    return(diffBetweenMinReductionRequiredAndPred)
}

eval_g0_eq <- function(x, EIR, minPrevReduc){
    #fixed_effi0.85_hl0.1667"
    x_input = t(c(EIR, 0.1667, 0.85, x, 0.04))
    # print(x_input)
    prediction <- predict(x = x_input, object = gp_result$GP_model)$mean
    if(prediction>100) {
        prediction = 100
    }
    diffBetweenMinReductionRequiredAndPred = minPrevReduc - prediction
    return(diffBetweenMinReductionRequiredAndPred)
}

get_param = function(X, GP_model) {
    return(X)
}

get_p_red = function(x, GP_model) {
    x_input = t(c(1, 0.1667, 0.85, x, 0.04))
    prev_red = predict(x = x_input, GP_model)$mean
    if(prev_red>100) {
        prev_red = 100
    }
    return(100-prev_red)
}

get_p_red_eq = function(x, GP_model) {
    x_input = t(c(1, 0.1667, 0.85, x, 0.04))
    prev_red = predict(x = x_input, GP_model)$mean
    if(prev_red>100) {
        prev_red = 100
    }
    return(prev_red)
}

get_p_red2 = function(X) {
    prev_red = predict(x = X, GP_model)$mean
    if(prev_red>100) {
        prev_red = 100
    }
    return(100-prev_red)
}

# Load GP model
gp_result_name = load(gp_file)
gp_result = get(gp_result_name)
rm(gp_result_name)

# Load parameter ranges
load(ranges_file)

lower_bounds = param_ranges[which(row.names(param_ranges)=="Coverage"),1]
upper_bounds = param_ranges[which(row.names(param_ranges)=="Coverage"),2]
x0 = lower_bounds+0.01

EIRValuesRange = seq(1,25,by=1)
# outcomeValuesRange =seq(0.1,0.95,by=0.05) 
outcomeValuesRange = c(0.10, 0.25, 0.5, 0.75, 0.95)
resultsShort = list()
resultsShort$EIR = NULL
resultsShort$prev.red = NULL
resultsShort$mincoverage = NULL
resultsShort$efficacy = NULL
resultsShort$halflife = NULL

fname="fixed_effi0.85_hl0.1667"
opts <- list("algorithm"="NLOPT_GN_ISRES",maxeval=5000,"xtol_rel"=1.0e-6)

desired_red = 90

for (EIR_lvl in 1:length(EIRValuesRange)){
    param_ranges2 = param_ranges
    param_ranges2["EIR",] = c(EIR_lvl, EIR_lvl)
    param_ranges2["Halflife",] = c(0.1667, 0.1667)
    param_ranges2["Efficacy",] = c(0.85, 0.85)
    param_ranges2["Access",] = c(0.04, 0.04)
    
    
    LB = param_ranges["Coverage",1]
    UB = param_ranges["Coverage",2]
    # LB[1] = EIR_lvl
    # UB[1] = EIR_lvl + 0.01
    # pars = LB
    # pars["Halflife"] = 0.1667
    # pars["Efficacy"] = 0.85
    # pars["Access"] = 0.04
    # pars[1] = EIR_lvl
    # fixed = c(1, 2, 3, 5)#, rep(0, length(LB)-1))
    ans = gosolnp(pars  = NULL, fixed = NULL, fun = get_param, eqfun = get_p_red_eq, eqB = c(70), LB = LB, UB = UB,
                  distr = rep(1, length(LB)), distr.opt = list(), n.restarts = 2, control = list(maxit = 100),
                  n.sim = 500, rseed = 443, GP_model = gp_result$GP_model)
    GP_model = gp_result$GP_model
    result = metaOpt(FUN = get_p_red2, optimType = "MIN", algorithm = "DA", numVar = 5, rangeVar = t(param_ranges2))
    
    # res_temp1 =  nloptr(x0=x0, eval_f = eval_f0, lb = lower_bounds, ub = upper_bounds,
    #                     eval_g_eq = eval_g0_eq, opts = opts, EIR = EIRValuesRange[EIR_lvl],
    #                     minPrevReduc = desired_red)
    # 
    # res_temp2 =  nloptr(x0=x0, eval_f = eval_f0, lb = lower_bounds, ub = upper_bounds,
    #                     eval_g_ineq = eval_g0, opts = opts, EIR = EIRValuesRange[EIR_lvl],
    #                     minPrevReduc = desired_red)
    # x1 = as.vector(res_temp1$solution) 
    # 
    # # now local minimisation,there is risk to fall into a local optimium,therefore run twice
    # res_temp = nloptr(x0=x1, eval_f= eval_f0, lb = lower_bounds, ub = upper_bounds, 
    #                   eval_g_ineq = eval_g0, opts = opts, EIR = EIRValuesRange[k], 
    #                   minPrevReduc = outcomeValuesRange[j])
    
    # resultsShort$EIR = c(resultsShort$EIR, EIRValuesRange[k])
    # resultsShort$prev.red = c(resultsShort$prev.red, outcomeValuesRange[j])
    # res_temp$solution = as.vector(res_temp$solution)
    # resultsShort$minvalue = c(resultsShort$minvalue,res_temp$solution[1])
}
resultsShort=as.data.frame(resultsShort)


