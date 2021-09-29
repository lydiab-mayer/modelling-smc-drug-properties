# Function returning the parameter to be optimized, 
get_param = function(x, GP_model, param_vec, param_name) {
  return(x)
}

# Function returning the predicted mean prevalence reduction for a given input
get_p_red = function(x, GP_model, param_vec, param_name) {
  param_vec[param_name] = x
  prev_red = predict(x = as.matrix(param_vec), GP_model)$mean
  if(prev_red>100) {
    prev_red = 100
  }
  return(prev_red)
}


# Function returning the predicted mean prevalence reduction - sd for a given input
get_p_red_sd_minus = function(x, GP_model, param_vec, param_name) {
  param_vec[param_name] = x
  prev_red = predict(x = as.matrix(param_vec), GP_model)
  prev_red = prev_red$mean - sqrt(prev_red$sd2 + prev_red$nugs)
  if(prev_red>100) {
    prev_red = 100
  }
  return(prev_red)
}

# Function returning the predicted mean prevalence reduction + sd for a given input
get_p_red_sd_plus = function(x, GP_model, param_vec, param_name) {
  param_vec[param_name] = x
  prev_red = predict(x = as.matrix(param_vec), GP_model)
  prev_red = prev_red$mean + sqrt(prev_red$sd2 + prev_red$nugs)
  if(prev_red>100) {
    prev_red = 100
  }
  return(prev_red)
}

# Function returning the predicted mean prevalence reduction - 2sd for a given input
get_p_red_sd2_minus = function(x, GP_model, param_vec, param_name) {
  param_vec[param_name] = x
  prev_red = predict(x = as.matrix(param_vec), GP_model)
  prev_red = prev_red$mean - sqrt(2*prev_red$sd2 + prev_red$nugs)
  if(prev_red>100) {
    prev_red = 100
  }
  return(prev_red)
}

# Function returning the predicted mean prevalence reduction + 2sd for a given input
get_p_red_sd2_plus = function(x, GP_model, param_vec, param_name) {
  param_vec[param_name] = x
  prev_red = predict(x = as.matrix(param_vec), GP_model)
  prev_red = prev_red$mean + sqrt(2*prev_red$sd2 + prev_red$nugs)
  if(prev_red>100) {
    prev_red = 100
  }
  return(prev_red)
}
