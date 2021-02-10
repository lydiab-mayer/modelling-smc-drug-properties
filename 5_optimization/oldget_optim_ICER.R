####################
# Find optimal intervention coverages for lowest ICER
#
# created 01.01.2019
# monica.golumbeanu@unibas.ch
####################
library(stringr)
library(metaheuristicOpt)
library(Rsolnp)

source("~/MMC/elimination_modeling/scripts/processing/cost_functions.R")

rep.row = function(x,n){
    return(matrix(rep(x,each=n),nrow=n))
}
rep.col<-function(x,n){
    matrix(rep(x,each=n), ncol=n, byrow=TRUE)
}

get_ICER = function(X) {
        prev_red = predict(x = X, object = prev_red_GP)
        total_cost = predict(x = X, object = cost_GP)
        ICER = (total_cost$mean/10000)/(prev_red$mean*100*(prev_red$mean>0.7))
        if (is.infinite(ICER)) ICER = 10000000
        return(ICER)
}

calc_cost = function(MDA_cov, IRS_cov, Access, RCD_nu, RCD_iota) {
    par = as.data.frame(cbind(Access, MDA_cov, IRS_cov, RCD_nu, RCD_iota))
    #CM_cost = predict(GP_CM, as.matrix(par))$mean
    total_cost = attach_CI_Costs(5, 10000, RCD_iota, RCD_nu) + 
        attach_MDA_Costs(MDA_years, 10000, 2, MDA_cov) +
        attach_IRS_Costs(IRS_years, 2, 10000, IRS_cov)
    return(total_cost)
}

get_prev_red = function(X) {
    prev_red = predict(x = X, object = prev_red_GP)
    #total_cost = predict(x = X, object = cost_GP)
    total_cost = calc_cost(X[1], X[2], X[3], X[4], X[5])
    if (total_cost/(10000*5)>=cost_constraint) {
        est_prev_red = -1
        # print(est_prev_red)
    }else {
        # print(total_cost)
        est_prev_red = prev_red$mean
    }
    return(est_prev_red)
}

get_prev_red_cost = function(X) {
    prev_red = predict(x = X, object = prev_red_GP)
    total_cost = predict(x = X, object = cost_GP)
    # total_cost = calc_cost(X[1], X[2], X[3], X[4], X[5])
    if (total_cost$mean/(10000*5)>=cost_constraint || prev_red$mean>1) {
        est_prev_red = -1
        # print(est_prev_red)
    }else {
        # print(total_cost)
        est_prev_red = prev_red$mean
    }
    # prev_red$mean * (total_cost$mean/(10000*5)>=cost_constraint)
    return(est_prev_red)
}

get_cost = function(X) {
    # total_cost = predict(x = X, object = cost_GP)
    # print(total_cost$mean/50000)
    selected_data_diff = cbind(selected_data[,1:5], rowSums(abs(selected_data[,1:5] - rep.row(X, nrow(selected_data)))), selected_data[6])
    ordered_data = selected_data_diff[order(selected_data_diff[,6]),]
    CM_cost = mean(ordered_data$CM_cost[1:10])
    non_CM_cost = calc_cost(X[1], X[2], X[3], X[4], X[5])
    total_cost = CM_cost + non_CM_cost
    # print(total_cost)
    return(total_cost/50000)
}

get_cost_GP = function(X) {
    total_cost = predict(x = X, object = cost_GP)
    # print(total_cost$mean/50000)
    # selected_data_diff = cbind(selected_data[,1:5], rowSums(abs(selected_data[,1:5] - rep.row(X, nrow(selected_data)))), selected_data[6])
    # ordered_data = selected_data_diff[order(selected_data_diff[,6]),]
    # CM_cost = mean(ordered_data$CM_cost[1:10])
    # non_CM_cost = calc_cost(X[1], X[2], X[3], X[4], X[5])
    # total_cost = CM_cost + non_CM_cost
    # # print(total_cost)
    return(total_cost$mean/50000)
}

get_p_red = function(X) {
    prev_red = predict(x = X, object = prev_red_GP)
    return((1-prev_red$mean))
}

get_p_red_ratio = function(X) {
    prev_red = predict(x = X, object = prev_red_GP)
    cost = get_cost(X) 
    return(prev_red$mean/cost)
}

train_dir = "~/MMC/elimination_modeling/results/GP/trained_models/prev_red_E14_RCD_after2/"
RCD_before_agg = read.table("MMC/elimination_modeling/results/results_tables/split_E14_RCD_after2/merged_agg.txt",
                            sep="\t", header = TRUE, as.is = TRUE)
RCD_before_agg$total_cost = RCD_before_agg$total_cost/(5*10000)
MDA_years = IRS_years = 1

file.names = dir(train_dir, pattern = "_cv.RData", full.names = TRUE)

ICER_df = NULL
prev_red_df = NULL
for(i in 1:length(file.names)){
    # load the prev red GP
    load(file.names[i])
    prev_red_GP = cv_result$GP_model
    # load the cost GP
    exp_name = (basename(file.names[i]))
    cost_file = str_replace(exp_name, "_cv.RData", "_cv_total_cost.RData")
    load(paste(train_dir, cost_file, sep=""))
    cost_GP = cv_result$GP_model
    param_ranges = rbind(c(0, 0.9), c(0, 0.9), c(0.0001798853, 0.5317248075), c(0, 50), c(0, 25))
    
    # For prevalence reduction under constraints
    for (cost_constraint in 5:10) {
        par = as.vector(lhs(1, param_ranges))
        
        # with Rsolnp
        n_par = 5
        LB = param_ranges[,1]
        UB = param_ranges[,2]
        eqLB = cost_constraint - 1
        eqUB = cost_constraint
        selected_data = RCD_before_agg[which(RCD_before_agg$EIR == cv_result$EIR & 
                                                 RCD_before_agg$Seasonality == cv_result$seasonality), c("Access", "MDA_cov", "IRS_cov", "RCD_nu", "RCD_iota", "CM_cost")]
        ans = gosolnp(pars  = NULL, fixed = NULL, fun = get_p_red, ineqfun = get_cost_GP,
                      ineqLB = eqLB, ineqUB = eqUB, LB = LB, UB = UB, control = list(outer.iter = 100, trace = 1),
                      distr = rep(1, length(LB)), distr.opt = list(), n.restarts = 2,
                      n.sim = 500, rseed = 443)
        est_p = predict( prev_red_GP, t(as.matrix(ans$pars)))
        est_cost = predict( cost_GP, t(as.matrix(ans$pars)))
        prev_red_df = rbind.data.frame(prev_red_df, cbind.data.frame(cost_constraint, cv_result$EIR, cv_result$seasonality, est_p$mean, est_cost$mean/50000,  t(ans$pars)))
        
        # result = metaOpt(get_p_red_ratio, optimType = "MIN", algorithm = "DA", numVar = 5, rangeVar = t(param_ranges))
        # result = optim(par, fn = get_p_red_ratio, method = "L-BFGS-B", lower = param_ranges[,1], upper = param_ranges[,2], hessian = FALSE, control = list(maxit=200000, fnscale=-1))
        
        # est_p = predict( prev_red_GP, t(as.matrix(result$par)))
        # est_cost = predict( cost_GP, t(as.matrix(result$par)))
        # # est_cost = calc_cost(result$par[1], result$par[2], result$par[3], result$par[4], result$par[5])
        # est_p = predict( prev_red_GP, (as.matrix(result$result)))
        # # est_cost = predict( cost_GP, (as.matrix(result$result)))
        # # cost_without_CM = calc_cost (result$result[1], result$result[2], result$result[3], result$result[4], result$result[5]) /50000
        # prev_red_df = rbind.data.frame(prev_red_df, cbind.data.frame(cost_constraint, cv_result$EIR, cv_result$seasonality, est_p$mean, est_cost$mean/50000,  t(result$par)))
    }
    print(i)
}


