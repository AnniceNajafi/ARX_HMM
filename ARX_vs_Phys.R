#' Authors@R: person("Annice", "Najafi", email = "annicenajafi27@gmail.com")
#' Fall 2025
#' The following script performs physics to ARX comparison
#'


#Load relevant libraries

library(tidyverse)
library(ggplot2)

###Note: will add roxygen docs for functions later

CFG <- list(
  dt = 0.02,
  T  = 20.0,
  seed = 11L,
  
  #M0 (LPF actuator)
  tau_cmd = 0.15,
  cmd_proc_std = 0.005,
  
  #M1 (mass-spring-damper plant)
  m = 1.0, c = 0.4, k = 4.0,
  plant_proc_std = 0.04,
  
  #Sensors
  tau_sa = 0.08, sa_meas_std = 0.03, sa_bias = 0.0,
  tau_sb = 0.25, sb_meas_std = 0.02, sb_bias = 0.01,
  
  #Fusion + output
  alpha = 0.6,
  w1 = 1.0, w2 = 0.2,
  
  #ARX orders
  train_frac = 0.5,
  na_cmd = 1L, nb_cmd = 1L,
  na_plant = 2L, nb_plant = 2L,
  na_sens = 1L, nb_sens = 1L
)

set.seed(CFG$seed)


make_u_cmd <- function(n, amp = 1.0, dt = 0.02) {
  t <- (0:(n-1)) * dt
  u <- 0.6 * sin(2*pi*0.25*t) + 0.4 * sign(sin(2*pi*0.05*t))
  step_every <- as.integer(2.5/dt)
  for (i in seq(1, n, by = step_every)) {
    u[i:n] <- u[i:n] + 0.2 * rnorm(1)
  }
  amp * u
}

lpf_first_order <- function(u, tau, dt, proc_std = 0.0) {
  n <- length(u)
  y <- numeric(n)
  a <- exp(-dt / max(tau, 1e-6))
  for (k in 2:n) {
    y[k] <- a * y[k-1] + (1 - a) * u[k-1] + proc_std * rnorm(1)
  }
  y
}

plant_msd <- function(u, dt, m, c, k, proc_std = 0.0) {
  n <- length(u)
  x_pos <- numeric(n)
  x_vel <- numeric(n)
  for (kk in 2:n) {
    acc <- (u[kk-1] - c * x_vel[kk-1] - k * x_pos[kk-1]) / m + proc_std * rnorm(1)
    x_vel[kk] <- x_vel[kk-1] + dt * acc
    x_pos[kk] <- x_pos[kk-1] + dt * x_vel[kk-1]
  }
  list(pos = x_pos, vel = x_vel)
}

sensor_first_order <- function(xpos, dt, tau, meas_std = 0.0, bias = 0.0) {
  y <- lpf_first_order(xpos, tau, dt, proc_std = 0.0)
  y <- y + bias + rnorm(length(y), sd = meas_std)
  y
}


fit_arx <- function(y, u, na = 1L, nb = 1L) {
  stopifnot(length(y) == length(u))
  kmax <- length(y)
  p <- max(na, nb)
  if (kmax <= p) stop("Series too short for ARX given na/nb.")
  rows <- vector("list", kmax - p)
  Tvec <- numeric(kmax - p)
  for (k in (p+1):kmax) {
    phi <- c(
      sapply(1:na, function(i) y[k - i]),
      sapply(1:nb, function(j) u[k - j])
    )
    rows[[k - p]] <- phi
    Tvec[k - p] <- y[k]
  }
  Phi <- do.call(rbind, rows)
  theta <- as.numeric(solve(t(Phi) %*% Phi, t(Phi) %*% Tvec))
  list(theta = theta, na = na, nb = nb)
}

sim_arx <- function(u, y_seed, model) {
  na <- model$na; nb <- model$nb; theta <- model$theta
  n <- length(u)
  yhat <- numeric(n)
  stopifnot(length(y_seed) >= na)
  yhat[1:na] <- tail(y_seed, na)
  p <- max(na, nb)
  for (k in (p+1):n) {
    phi <- c(
      sapply(1:na, function(i) yhat[k - i]),
      sapply(1:nb, function(j) u[k - j])
    )
    yhat[k] <- sum(theta * phi)
  }
  yhat
}


n <- as.integer(CFG$T / CFG$dt)
t <- (0:(n-1)) * CFG$dt
u_cmd <- make_u_cmd(n, dt = CFG$dt)


u_act_P <- lpf_first_order(u_cmd, CFG$tau_cmd, CFG$dt, proc_std = CFG$cmd_proc_std)
plant_P <- plant_msd(u_act_P, CFG$dt, CFG$m, CFG$c, CFG$k, proc_std = CFG$plant_proc_std)
x_pos_P <- plant_P$pos
x_vel_P <- plant_P$vel
yA_P <- sensor_first_order(x_pos_P, CFG$dt, CFG$tau_sa, CFG$sa_meas_std, CFG$sa_bias)
yB_P <- sensor_first_order(x_pos_P, CFG$dt, CFG$tau_sb, CFG$sb_meas_std, CFG$sb_bias)
y_fused_P <- CFG$alpha * yA_P + (1 - CFG$alpha) * yB_P
z_P <- CFG$w1 * x_pos_P + CFG$w2 * x_vel_P


split <- as.integer(CFG$train_frac * n)
arx_cmd   <- fit_arx(u_act_P[1:split], u_cmd[1:split], na = CFG$na_cmd,   nb = CFG$nb_cmd)
arx_plant <- fit_arx(x_pos_P[1:split], u_act_P[1:split], na = CFG$na_plant, nb = CFG$nb_plant)
arx_sa    <- fit_arx(yA_P[1:split],    x_pos_P[1:split], na = CFG$na_sens,  nb = CFG$nb_sens)
arx_sb    <- fit_arx(yB_P[1:split],    x_pos_P[1:split], na = CFG$na_sens,  nb = CFG$nb_sens)

u_act_A <- sim_arx(u_cmd,     y_seed = u_act_P[1:CFG$na_cmd],   model = arx_cmd)
xpos_A  <- sim_arx(u_act_A,   y_seed = x_pos_P[1:CFG$na_plant], model = arx_plant)
xvel_A  <- c(0.0, diff(xpos_A) / CFG$dt)
yA_A    <- sim_arx(x_pos_P,   y_seed = yA_P[1:CFG$na_sens],     model = arx_sa)
yB_A    <- sim_arx(x_pos_P,   y_seed = yB_P[1:CFG$na_sens],     model = arx_sb)
y_fused_A <- CFG$alpha * yA_A + (1 - CFG$alpha) * yB_A
z_A <- CFG$w1 * xpos_A + CFG$w2 * xvel_A


df <- tibble(
  t = t,

  M0_u_act_P = u_act_P, M0_u_act_A = u_act_A,
  M1_x_pos_P = x_pos_P, M1_x_pos_A = xpos_A,
  M2a_yA_P = yA_P, M2a_yA_A = yA_A,
  M2b_yB_P = yB_P, M2b_yB_A = yB_A,
  M3_y_fused_P = y_fused_P, M3_y_fused_A = y_fused_A,
  M4_z_P = z_P, M4_z_A = z_A
)


pipeline_levels <- c("M0","M1","M2a","M2b","M3","M4")
modules <- tribble(
  ~module, ~truth_col,      ~approx_col,     ~label,
  "M0",    "M0_u_act_P",    "M0_u_act_A",    "Actuator output (M0)",
  "M1",    "M1_x_pos_P",    "M1_x_pos_A",    "Plant position (M1)",
  "M2a",   "M2a_yA_P",      "M2a_yA_A",      "Sensor A output (M2a)",
  "M2b",   "M2b_yB_P",      "M2b_yB_A",      "Sensor B output (M2b)",
  "M3",    "M3_y_fused_P",  "M3_y_fused_A",  "Fused measurement (M3)",
  "M4",    "M4_z_P",        "M4_z_A",        "Final output z (M4)"
) %>% mutate(module = factor(module, levels = pipeline_levels))



plot_df <- purrr::pmap_dfr(
  modules,
  function(module, truth_col, approx_col, label) {
    n <- nrow(df)
    mlev <- factor(module, levels = pipeline_levels)
    approx_tb <- tibble(t = df$t, module = mlev, label = label,
                        series = factor("approx", levels=c("approx","truth")),
                        value  = df[[approx_col]])
    truth_tb  <- tibble(t = df$t, module = mlev, label = label,
                        series = factor("truth", levels=c("approx","truth")),
                        value  = df[[truth_col]])
    bind_rows(approx_tb, truth_tb)
  }
)


col_truth  <- "#0C2B4E"
col_approx <- "#A1C2BD" 

plot_df <- plot_df %>%
  dplyr::mutate(series = factor(series, levels = c("truth","approx")))

p <- ggplot(plot_df, aes(t, value, color = series, linetype = series, group = series)) +
  geom_line(linewidth = 2) +
  facet_wrap(~ module + label, scales = "free_y", ncol = 2) +
  labs(title = "Physics (truth) vs ARX (approx) across modules",
       x = "Time (s)", y = "Signal", color = NULL, linetype = NULL) +
  scale_color_manual(values = c(truth = col_truth, approx = col_approx)) +
  scale_linetype_manual(values = c(truth = "solid", approx = "dashed")) +
  theme_minimal(base_size = 14) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color = "black"),
        plot.title = element_text(hjust = 0.5, size = 20),
        strip.text = element_text(size = 12))
print(p)


















