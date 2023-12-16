# ==================
# Lice model fitting
# ==================

fitLiceStream <- function(
    data = NA,
    order = list(),
    control = list(), 
    map = list(),
    silent = T
) {
  # Misc
  num_locs <- nrow(data$v)
  
  # Defined parameter based on init. Either from file or to warmstart
  parameters <- list(
    log_nu_AF = log(0.05),
    log_nu_MB = log(0.05),
    log_nu_ST = log(0.05),
    log_iota_I = 0,
    log_iota_H = 0,
    log_iota_N = 0,
    beta_T = 0,
    log_alpha = log(2.0),
    log_neg_rho = -1,
    log_neg_delta_Vmec = rep(-1, length.out = order$m_mec),
    log_neg_delta_Vmed = rep(-1, length.out = order$m_med),
    log_neg_delta_Vclf = rep(log(0.35 / (1e-8 + order$m_clf)), length.out = order$m_clf),
    log_c = log(1e-3),
    theta_scale = 0.5,
    theta_rate = -0.25
  )
  data$q <- order$q
  if (order$q == 0){
    map$theta_scale <- factor(NA)
    map$theta_rate <- factor(NA)
  }
  data$p <- order$p
  
  model <- MakeADFun(data = data, parameters = parameters, DLL = "nll_lice", map = map, silent = silent)
  
  # Minimize negative log-likelihood
  opt <- nlminb(model$par, model$fn, model$gr, control = control)
  rprt <- model$env$report() # Report of chosen estimates
  
  return(list(model = model, opt = opt, rprt = rprt))
}

# =======================
# Treatment model fitting
# =======================

fitTreat <- function(
    data = NA, 
    order = list(), 
    map = list(),
    parameters = list(),
    silent,
    control = list()
){
  # Define map and parameters as empty list and fill in accordingly
  
  # Misc
  nsites <- nrow(data$v)
  
  # parameters
  parameters = list(
    beta_0 = 0.0,
    beta_low = 0.0,
    phi_scale = 0.0,
    phi_rate = 0.0,
    theta_scale = 0.0,
    theta_rate = 0.0,
    # log_d = -3,
    gamma_AF = rep(0.0, length.out = order$r_AF),
    gamma_MB = rep(0.0, length.out = order$r_MB),
    gamma_ST = rep(0.0, length.out = order$r_ST),
    gamma_I_AF = 0.0,
    gamma_I_MB = 0.0,
    gamma_I_ST = 0.0,
    pi_AF = 1.0,
    pi_MB = 1.0,
    pi_ST = 1.0
  )
  data$p = order$p
  data$q = order$q
  
  # Fix parameters based on specification
  if (order$r_AF == 0){
    map$pi_AF = factor(NA)
    map$gamma_I_AF = factor(NA)
  }
  if (order$r_MB == 0){
    map$pi_MB = factor(NA)
    map$gamma_I_MB = factor(NA)
  }
  if (order$r_ST == 0){
    map$pi_ST = factor(NA)
    map$gamma_I_ST = factor(NA)
  }
  if (order$p == 0){
    map$phi_scale = factor(NA)
    map$phi_rate = factor(NA)
  }
  if (order$q == 0){
    map$theta_scale = factor(NA)
    map$theta_rate = factor(NA)
  }
  
  model <- MakeADFun(data = data, parameters = parameters, DLL = "nll_treat", map = map, silent = silent)
  
  # Minimize negative log-likelihood
  opt <- nlminb(model$par, model$fn, model$gr, control = control)
  rprt <- model$env$report() # Report of chosen estimates
  
  return(list(model = model, opt = opt, rprt = rprt))
}

# ===============
# Other utilities
# ===============

# Biological functions
development_time_stage <- function (temp) {
  return(exp(3.708396114525632 + -0.2322266272453008 * temp +  0.004496564435398739 * temp^2) / 7.0)
}
development_time_male <- function (temp) {
  return(exp(3.667550690993821 + -0.24720286112431838 * temp +  0.005024716188616198 * temp^2) / 7.0)
}
development_time_female <- function (temp) {
  return(exp(3.7492415380574373 + -0.21725039336628219 * temp +  0.003968412682181245 * temp^2) / 7.0)
}
development_time_larval <- function (temp){
  return(exp(3.290559408354516 + -0.03497494410938104 * temp + -0.0012845611962801293 * temp^2) / 7.0)
}
egg_count_normalised <- function (temp){
  return(exp(4.834458332510026 +  0.15614966081868187 * temp + -0.00769948277659029 * temp^2) / 300)
}

# Goodness-of-fit
gof <- function(res) {
  NLL <- res$opt$objective
  AIC <- 2 * (NLL + length(res$model$par))
  BIC <- 2 * NLL + length(res$model$par) * log(sum(res$model$env$data$A))
  gof <- list(
    NLL = NLL,
    AIC = AIC,
    BIC = BIC
  )  
  return(gof)  
}

# Split data
time_slice <- function (data_in, t_strt, t_stop, mdl_kind) {
  if (mdl_kind == "lice"){
    data_out <- list(
      y_AF = data_in$y_AF[, t_strt:t_stop],
      y_MB = data_in$y_MB[, t_strt:t_stop],
      y_ST = data_in$y_ST[, t_strt:t_stop],
      y_RN = data_in$y_RN[, t_strt:t_stop],
      y_RS = data_in$y_RS[, t_strt:t_stop],
      fish_counted = data_in$fish_counted[, t_strt:t_stop],
      A = data_in$A[, t_strt:t_stop],
      temperature = data_in$temperature[, t_strt:t_stop],
      dev_time_stage = data_in$dev_time_stage[, t_strt:t_stop],
      dev_time_ml = data_in$dev_time_ml[, t_strt:t_stop],
      dev_time_fm = data_in$dev_time_fm[, t_strt:t_stop],
      dev_time_larval = data_in$dev_time_larval[, t_strt:t_stop],
      new_hatches = data_in$new_hatches[, t_strt:t_stop],
      v_med = data_in$v_med[, t_strt:t_stop],
      v_mec = data_in$v_mec[, t_strt:t_stop],
      v_clf = data_in$v_clf[, t_strt:t_stop]
    )
  } else if (mdl_kind == "treat") {
    data_out <- list(
      v = data_in$v[, t_strt:t_stop],
      y_AF = data_in$y_AF[, t_strt:t_stop],
      y_MB = data_in$y_MB[, t_strt:t_stop],
      y_ST = data_in$y_ST[, t_strt:t_stop],
      A = data_in$A[, t_strt:t_stop],
      temperature = data_in$temperature[, t_strt:t_stop],
      low_limit = data_in$low_limit[, t_strt:t_stop],
      lice_limit = data_in$lice_limit[, t_strt:t_stop]
    )
  }
  return(data_out) 
}