{
  library(TMB)
  library(TMBhelper)
  library(tidyverse)
  library(magrittr)
  library(readr)
  library(jsonlite)
  library(boot)
  source("utils.R")
}

{ # Load data files
  load("data/report_dtl_PO03.RData")
  num_dates <- length(ds$dates)
  num_locs <- length(ds$locs)
  # These sites are excluded during estimation, but contribute to neighbor infestation
  excluded_locNos <- fromJSON("data/excluded_sites.json") 
  loc_indi <- setdiff(1:num_locs, match(excluded_locNos, ds$locs))
  t_strt <- 289
  t_stop <- 460
}

{ # Construct TMB dataset
  data <- list(
    y_AF = ds$lice_adult_female[loc_indi,],
    y_MB = ds$lice_mobile[loc_indi,],
    y_ST = ds$lice_stationary[loc_indi,],
    y_RN = ds$lice_recruits_neighbour[loc_indi,],
    y_RS = ds$lice_recruits_self[loc_indi,],
    fish_counted = ds$fish_counted_cage[loc_indi,] * ds$capacity[loc_indi,] / 790,
    A = 1 * ds$is_active[loc_indi,],
    temperature = ds$temperature[loc_indi,],
    dev_time_stage = development_time_stage(ds$temperature[loc_indi,]),
    dev_time_ml = development_time_male(ds$temperature[loc_indi,]),
    dev_time_fm = development_time_female(ds$temperature[loc_indi,]),
    dev_time_larval = development_time_larval(ds$temperature[loc_indi,]),
    new_hatches = egg_count_normalised(ds$temperature[loc_indi,]),
    v_med = ds$treat_medi[loc_indi,],
    v_mec = ds$treat_mech[loc_indi,],
    v_clf = ds$count_clf[loc_indi,] / ds$capacity[loc_indi,],
    v = 1 * (ds$treat_medi[loc_indi,] | ds$treat_mech[loc_indi,]),
    low_limit = 1 * (ds$lice_limit[loc_indi,] < 0.5),
    lice_limit = ds$lice_limit[loc_indi,]
  )
}

{ # Compile likelihood functions
  compile("nll_lice.cpp")
  dyn.load(dynlib("nll_lice"))
  compile("nll_treat.cpp")
  dyn.load(dynlib("nll_treat"))
}

# =====================
# Lice model estimation
# =====================

fit.lice <- fitLiceStream(
  data = time_slice(data, t_strt, t_stop, "lice"),
  order = list(q = 20, m_mec = 1, m_med = 5, m_clf = 5),
  map = list( # NB: We fix the effect of treatment fish after 2018.
    log_neg_delta_Vclf = factor(rep(NA, 5))), 
  control = list(trace = 1),
  silent = F)

# Goodness-of-fit
gof(fit.lice)

{ # Get estimates of standard deviations by Delta method
  sdrep.lice <- sdreport(fit.lice$model)
  summary(sdrep.lice, "report")
}

# Collect in-sample data from estimation (the last time index is out-of-sample)
df.lice <- data.frame( 
  date = rep(ds$dates[(t_strt + 20):(t_stop - 1)], each = length(loc_indi)), 
  locNo = rep(ds$locs[loc_indi], times = t_stop - t_strt - 20),
  mu_AF = c(fit.lice$rprt$mus_AF[,(1 + 20):(t_stop - t_strt)]),
  mu_MB = c(fit.lice$rprt$mus_MB[,(1 + 20):(t_stop - t_strt)]),
  mu_ST = c(fit.lice$rprt$mus_ST[,(1 + 20):(t_stop - t_strt)]),
  mu_bar_AF = c(fit.lice$rprt$mus_bar_AF[,(1 + 20):(t_stop - t_strt)]),
  mu_bar_MB = c(fit.lice$rprt$mus_bar_MB[,(1 + 20):(t_stop - t_strt)]),
  mu_bar_ST = c(fit.lice$rprt$mus_bar_ST[,(1 + 20):(t_stop - t_strt)]),
  surprise_AF = c(fit.lice$rprt$surprise_AF[,(1 + 20):(t_stop - t_strt)]),
  surprise_MB = c(fit.lice$rprt$surprise_MB[,(1 + 20):(t_stop - t_strt)]),
  surprise_ST = c(fit.lice$rprt$surprise_ST[,(1 + 20):(t_stop - t_strt)]),
  innovation_AF = c(fit.lice$rprt$innovation_AF[,(1 + 20):(t_stop - t_strt)]),
  innovation_MB = c(fit.lice$rprt$innovation_MB[,(1 + 20):(t_stop - t_strt)]),
  innovation_ST = c(fit.lice$rprt$innovation_ST[,(1 + 20):(t_stop - t_strt)]),
  lambda_N_AF = c(fit.lice$rprt$lambdas_N_AF[,(1 + 20):(t_stop - t_strt)]),
  lambda_N_MB = c(fit.lice$rprt$lambdas_N_MB[,(1 + 20):(t_stop - t_strt)]),
  lambda_N_ST = c(fit.lice$rprt$lambdas_N_ST[,(1 + 20):(t_stop - t_strt)]),
  lambda_H_AF = c(fit.lice$rprt$lambdas_H_AF[,(1 + 20):(t_stop - t_strt)]),
  lambda_H_MB = c(fit.lice$rprt$lambdas_H_MB[,(1 + 20):(t_stop - t_strt)]),
  lambda_H_ST = c(fit.lice$rprt$lambdas_H_ST[,(1 + 20):(t_stop - t_strt)]),
  lambda_I_AF = c(fit.lice$rprt$lambdas_I_AF[,(1 + 20):(t_stop - t_strt)]),
  lambda_I_MB = c(fit.lice$rprt$lambdas_I_MB[,(1 + 20):(t_stop - t_strt)]),
  lambda_I_ST = c(fit.lice$rprt$lambdas_I_ST[,(1 + 20):(t_stop - t_strt)]),
  nu_AF = fit.lice$rprt$nu_AF,
  nu_MB = fit.lice$rprt$nu_MB,
  nu_ST = fit.lice$rprt$nu_ST,
  y_AF = c(data$y_AF[,(t_strt + 20):(t_stop - 1)]),
  y_MB = c(data$y_MB[,(t_strt + 20):(t_stop - 1)]),
  y_ST = c(data$y_ST[,(t_strt + 20):(t_stop - 1)]),
  y_RN = c(data$y_RN[,(t_strt + 20):(t_stop - 1)]),
  y_RS = c(data$y_RS[,(t_strt + 20):(t_stop - 1)]),
  A = c(data$A[,(t_strt + 20):(t_stop - 1)]),
  n = c(data$fish_counted[,(t_strt + 20):(t_stop - 1)]))

# Example plot of results
df.lice %>% filter(locNo == 12019) %>% ggplot(aes(x = date)) +
  geom_line(aes(y = y_AF), color = "black") +
  geom_line(aes(y = mu_AF), color = "red")

df.lice %>% filter(locNo == 12019) %>% ggplot(aes(x = date)) +
  geom_line(aes(y = y_MB), color = "black") +
  geom_line(aes(y = mu_MB), color = "yellow")

df.lice %>% filter(locNo == 12019) %>% ggplot(aes(x = date)) +
  geom_line(aes(y = y_MB), color = "black") +
  geom_line(aes(y = mu_MB), color = "green")

# ==========================
# Treatment model estimation
# ==========================

fit.treat <- fitTreat(
  data = time_slice(data, t_strt, t_stop, "treat"),
  order = list(p = 0, q = 0, r_AF = 1, r_MB = 1, r_ST = 0),
  map = list(),
  control = list(trace = 0, eval.max = 3000, iter.max = 1500),
  silent = T)

# Goodness-of-fit
gof(fit.treat)

{ # Get estimates of standard deviations by Delta method
  sdrep.treat <- sdreport(fit.treat$model)
  summary(sdrep.treat, "report")
}

# Collect in-sample data from estimation (the last time index is out-of-sample)
df.treat <- data.frame(
  date = rep(ds$dates[(t_strt + 20):(t_stop - 1)], each = length(loc_indi)),
  locNo = rep(ds$locs[loc_indi], times = t_stop - t_strt - 20),
  treatment = 1 * c(fit.treat$model$env$data$v[,(1 + 20):(t_stop - t_strt)]),
  active = 1 * c(fit.treat$model$env$data$A[,(1 + 20):(t_stop - t_strt)]),
  lice_AF = c(fit.treat$model$env$data$y_AF[,(1 + 20):(t_stop - t_strt)]),
  lice_MB = c(fit.treat$model$env$data$y_MB[,(1 + 20):(t_stop - t_strt)]),
  lice_ST = c(fit.treat$model$env$data$y_ST[,(1 + 20):(t_stop - t_strt)]),
  prob_treatment = c(fit.treat$rprt$prob[,(1 + 20):(t_stop - t_strt)]),
  eta_surprise = c(fit.treat$rprt$eta_surprise[,(1 + 20):(t_stop - t_strt)]),
  eta_innovation = c(fit.treat$rprt$eta_innovation[,(1 + 20):(t_stop - t_strt)]),
  eta = c(fit.treat$rprt$eta[,(1 + 20):(t_stop - t_strt)]),
  eta_bar = c(fit.treat$rprt$eta_bars[,(1 + 20):(t_stop - t_strt)]))

# Example plot of results
ggplot(mapping = aes(x = date), data = df.treat %>% filter(locNo == 12019)) +
  geom_line(aes(y = prob_treatment), color = "red") +
  geom_vline(xintercept = df.treat$date[df.treat$treatment == 1 & df.treat$locNo == 12019], color = "black")
