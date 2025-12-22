#' Authors@R: person("Annice", "Najafi", email = "annicenajafi27@gmail.com")
#' Fall 2025
#' 
#' BIC vs states

suppressPackageStartupMessages({
  library(tidyverse)
  library(depmixS4)
})

set.seed(11)


params <- list(
  dt   = 0.02,
  T    = 20,
  tau_cmd      = 0.15,
  cmd_proc_std = 0.005,

  m = 1.0, c = 0.4, k = 4.0, plant_proc_std = 0.04,
  #sensors
  tau_sa = 0.08, sa_meas_std = 0.03, sa_bias = 0.0,
  tau_sb = 0.25, sb_meas_std = 0.02, sb_bias = 0.01,
  #fusion + performance
  alpha = 0.6,
  w1 = 1.0, w2 = 0.2,
  #ARX training split / orders
  train_frac = 0.5,
  na_cmd   = 1L, nb_cmd   = 1L,
  na_plant = 2L, nb_plant = 2L,
  na_sens  = 1L, nb_sens  = 1L
)


make_u_cmd <- function(n, amp = 1.0, dt = 0.02) {
  t <- (0:(n-1)) * dt
  u <- 0.6 * sin(2*pi*0.25*t) + 0.4 * sign(sin(2*pi*0.05*t))
  step_every <- as.integer(2.5/dt)
  for (i in seq(1, n, by = step_every)) u[i:n] <- u[i:n] + 0.2 * rnorm(1)
  amp * u
}

lpf_first_order <- function(u, tau, dt, proc_std = 0) {
  n <- length(u); y <- numeric(n)
  a <- exp(-dt / max(tau, 1e-6))
  for (k in 2:n) y[k] <- a*y[k-1] + (1 - a)*u[k-1] + proc_std*rnorm(1)
  y
}

plant_msd <- function(u, dt, m, c, k, proc_std = 0) {
  n <- length(u); x <- numeric(n); v <- numeric(n)
  for (kk in 2:n) {
    acc <- (u[kk-1] - c*v[kk-1] - k*x[kk-1]) / m + proc_std*rnorm(1)
    v[kk] <- v[kk-1] + dt*acc
    x[kk] <- x[kk-1] + dt*v[kk-1]
  }
  list(pos = x, vel = v)
}

sensor_first_order <- function(xpos, dt, tau, meas_std = 0, bias = 0) {
  y <- lpf_first_order(xpos, tau, dt, 0)
  y + bias + rnorm(length(y), sd = meas_std)
}

fit_arx <- function(y, u, na, nb) {
  stopifnot(length(y) == length(u))
  kmax <- length(y); p <- max(na, nb)
  if (kmax <= p) stop("ARX series too short.")
  rows <- vector("list", kmax - p); Tvec <- numeric(kmax - p)
  for (k in (p+1):kmax) {
    phi <- c(sapply(1:na, function(i) y[k-i]),
             sapply(1:nb, function(j) u[k-j]))
    rows[[k - p]] <- phi
    Tvec[k - p]  <- y[k]
  }
  Phi   <- do.call(rbind, rows)
  theta <- as.numeric(solve(t(Phi) %*% Phi, t(Phi) %*% Tvec))
  list(theta = theta, na = na, nb = nb)
}

sim_arx <- function(u, y_seed, model) {
  na <- model$na; nb <- model$nb; theta <- model$theta
  n <- length(u); yhat <- numeric(n)
  stopifnot(length(y_seed) >= na)
  yhat[1:na] <- tail(y_seed, na)
  p <- max(na, nb)
  for (k in (p+1):n) {
    phi <- c(sapply(1:na, function(i) yhat[k-i]),
             sapply(1:nb, function(j) u[k-j]))
    yhat[k] <- sum(theta * phi)
  }
  yhat
}


n <- as.integer(params$T / params$dt)
t <- (0:(n-1)) * params$dt
u_cmd <- make_u_cmd(n, dt = params$dt)


u_act_P <- lpf_first_order(u_cmd, params$tau_cmd, params$dt, params$cmd_proc_std)
pl      <- plant_msd(u_act_P, params$dt, params$m, params$c, params$k, params$plant_proc_std)
x_pos_P <- pl$pos; x_vel_P <- pl$vel
yA_P    <- sensor_first_order(x_pos_P, params$dt, params$tau_sa, params$sa_meas_std, params$sa_bias)
yB_P    <- sensor_first_order(x_pos_P, params$dt, params$tau_sb, params$sb_meas_std, params$sb_bias)
y_fused_P <- params$alpha*yA_P + (1 - params$alpha)*yB_P
z_P <- params$w1*x_pos_P + params$w2*x_vel_P


split <- as.integer(params$train_frac * n)
arx_cmd   <- fit_arx(u_act_P[1:split], u_cmd[1:split],     params$na_cmd,   params$nb_cmd)
arx_plant <- fit_arx(x_pos_P[1:split], u_act_P[1:split],   params$na_plant, params$nb_plant)
arx_sa    <- fit_arx(yA_P[1:split],    x_pos_P[1:split],   params$na_sens,  params$nb_sens)
arx_sb    <- fit_arx(yB_P[1:split],    x_pos_P[1:split],   params$na_sens,  params$nb_sens)


u_act_A <- sim_arx(u_cmd,   y_seed = u_act_P[1:params$na_cmd],   arx_cmd)
xpos_A  <- sim_arx(u_act_A, y_seed = x_pos_P[1:params$na_plant], arx_plant)
xvel_A  <- c(0, diff(xpos_A) / params$dt)
yA_A    <- sim_arx(xpos_A, y_seed = yA_P[1:params$na_sens],     arx_sa)
yB_A    <- sim_arx(xpos_A, y_seed = yB_P[1:params$na_sens],     arx_sb)
y_fused_A <- params$alpha*yA_A + (1 - params$alpha)*yB_A
z_A <- params$w1*xpos_A + params$w2*xvel_A


E <- tibble(
  e0  = u_act_A   - u_act_P,
  e1  = xpos_A    - x_pos_P,
  e2a = yA_A      - yA_P,
  e2b = yB_A      - yB_P,
  e3  = y_fused_A - y_fused_P,
  e4  = z_A       - z_P
)


Ez <- E %>% mutate(across(everything(), ~ as.numeric(scale(.x))))


fit_oneK <- function(K) {
  set.seed(123 + K)
  
  fams <- rep(list(gaussian()), 6)
  mod <- depmix(
    response = list(e0~1, e1~1, e2a~1, e2b~1, e3~1, e4~1),
    data     = as.data.frame(Ez),
    nstates  = K,
    family   = fams
  )
  fm <- fit(mod, verbose = FALSE)
  
  
  bic_val <- tryCatch(
    as.numeric(stats::BIC(fm)),
    error = function(e) {
      
      ll  <- as.numeric(stats::logLik(fm))          
      dfp <- attr(stats::logLik(fm), "df")          
      n   <- nrow(Ez)                                
      -2*ll + dfp*log(n)
    }
  )
  
  tibble(K = K, BIC = bic_val)
}


Ks <- 2:9
bic_tbl <- purrr::map_dfr(Ks, fit_oneK)

print(bic_tbl)


ggplot(bic_tbl, aes(K, BIC)) +
  geom_point(color="#547792", shape=8, stroke=1.5) +
  geom_line(color="#547792", size=2) +
  labs(title = "BIC vs number of HMM states (K)", x = "K", y = "BIC") +
  theme_minimal(base_size = 16)+
  theme(panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color = "black"),
        plot.title = element_text(hjust = 0.5, size = 20),
        strip.text = element_text(size = 12))

















