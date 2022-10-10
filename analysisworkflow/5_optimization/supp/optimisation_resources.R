# Function returning the parameter to be optimized, 
get_param = function(x, GP_model, param_vec, param_name) {
    return(x)
}



# Function returning the predicted mean probability for non-inferiority
get_prob_noninf= function(x, GP_model, param_vec, param_name) {
  param_vec[which(names(param_vec)==param_name)] = x
  prob_noninf = predict(x = as.matrix(param_vec), GP_model)$mean
  if(prob_noninf>1) {
    prob_noninf = 1
  }
  return(prob_noninf)
}

# Function returning the predicted mean prevalence reduction - sd for a given input
get_prob_noninf_sd_minus = function(x, GP_model, param_vec, param_name) {
  param_vec[which(names(param_vec)==param_name)] = x
  pred = predict(x = as.matrix(param_vec), GP_model)
  prob_noninf = pred$mean - sqrt(pred$sd2 + pred$nugs)
  if(prob_noninf>1) {
    prob_noninf = 1
  }
  return(prob_noninf)
}

# Function returning the predicted mean prevalence reduction + sd for a given input
get_prob_noninf_sd_plus = function(x, GP_model, param_vec, param_name) {
  param_vec[which(names(param_vec)==param_name)] = x
  pred = predict(x = as.matrix(param_vec), GP_model)
  prob_noninf = pred$mean + sqrt(pred$sd2 + pred$nugs)
  if(prob_noninf>1) {
    prob_noninf = 1
  }
  return(prob_noninf)
}

# Function returning the predicted mean prevalence reduction - 2sd for a given input
get_prob_noninf_sd2_minus = function(x, GP_model, param_vec, param_name) {
  param_vec[which(names(param_vec)==param_name)] = x
  pred = predict(x = as.matrix(param_vec), GP_model)
  prob_noninf = pred$mean - sqrt(2*pred$sd2 + pred$nugs)
  if(prob_noninf>1) {
    prob_noninf = 1
  }
  return(prob_noninf)
}

# Function returning the predicted mean prevalence reduction + 2sd for a given input
get_prob_noninf_sd2_plus = function(x, GP_model, param_vec, param_name) {
  param_vec[which(names(param_vec)==param_name)] = x
  pred = predict(x = as.matrix(param_vec), GP_model)
  prob_noninf = pred$mean + sqrt(2*pred$sd2 + pred$nugs)
  if(prob_noninf>1) {
    prob_noninf = 1
  }
  return(prob_noninf)
}
