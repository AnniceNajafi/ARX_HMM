#' Authors@R: person("Annice", "Najafi", email = "annicenajafi27@gmail.com")
#' Fall 2025
#' Digital twin with ground-truth faults + HMM segmentation

#We inject faults according to a schedule - this is related to figure 8 of the article


suppressPackageStartupMessages({
  library(tidyverse)
  library(depmixS4)
  library(clue)        #for solve_LSAP (hungarian)
})

params <- list(
  dt   = 0.02,
  T    = 20,
  seed = 11L,
  #actuator LPF
  tau_cmd      = 0.15,
  cmd_proc_std = 0.005,
  #plant (mass-spring-damper)
  m = 1.0, c = 0.4, k = 4.0, plant_proc_std = 0.04,
  #sensors (nominal)
  tau_sa = 0.08, sa_meas_std = 0.03, sa_bias = 0.0,
  tau_sb = 0.25, sb_meas_std = 0.02, sb_bias = 0.01,
  #fusion + performance
  alpha = 0.6, w1 = 1.0, w2 = 0.2,
  #ARX (train on first half)
  train_frac = 0.5, na_cmd = 1L, nb_cmd = 1L,
  na_plant = 2L, nb_plant = 2L,
  na_sens  = 1L, nb_sens  = 1L,
  #HMM states, note that this is not necessarily optimal as we explain the discussion section
  K = 4
)

#Fault schedule (Nominal --> SensorNoisy --> DynamicsOff --> Drift)
schedule <- tibble::tribble(
  ~regime,        ~t_start, ~t_end,
  "Nominal",         0.0,      5.0,
  "SensorNoisy",     5.0,     10.0,
  "DynamicsOff",    10.0,     14.5,
  "Drift",          14.5,     params$T
)


make_u_cmd <- function(n, amp = 1.0, dt = 0.02) {
  t <- (0:(n-1)) * dt
  u <- 0.6 * sin(2*pi*0.25*t) + 0.4 * sign(sin(2*pi*0.05*t))
  step_every <- as.integer(2.5/dt)
  for (i in seq(1, n, by = step_every)) u[i:n] <- u[i:n] + 0.2 * rnorm(1)
  amp * u
}
lpf_first_order <- function(u, tau, dt, proc_std = 0.0) {
  n <- length(u); y <- numeric(n)
  a <- exp(-dt / max(tau, 1e-6))
  for (k in 2:n) y[k] <- a*y[k-1] + (1 - a)*u[k-1] + proc_std*rnorm(1)
  y
}
plant_msd <- function(u, dt, m, c, k, proc_std = 0.0) {
  n <- length(u); x <- numeric(n); v <- numeric(n)
  for (kk in 2:n) {
    acc <- (u[kk-1] - c*v[kk-1] - k*x[kk-1]) / m + proc_std*rnorm(1)
    v[kk] <- v[kk-1] + dt*acc
    x[kk] <- x[kk-1] + dt*v[kk-1]
  }
  list(pos = x, vel = v)
}
sensor_first_order <- function(xpos, dt, tau, meas_std = 0.0, bias = 0.0) {
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


simulate_with_faults <- function(params, schedule) {
  set.seed(params$seed)
  n <- as.integer(params$T/params$dt)
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
  xvel_A  <- c(0, diff(xpos_A)/params$dt)
  yA_A    <- sim_arx(xpos_A, y_seed = yA_P[1:params$na_sens],     arx_sa)
  yB_A    <- sim_arx(xpos_A, y_seed = yB_P[1:params$na_sens],     arx_sb)
  y_fused_A <- params$alpha*yA_A + (1 - params$alpha)*yB_A
  z_A <- params$w1*xpos_A + params$w2*xvel_A
  
  #Apply ground-truth regimes to create mismatches (residual structure)
  #SensorNoisy: inflate sensor noises
  #DynamicsOff: weaken plant dynamics (k -> small)
  #Drift: add slow bias to sensor A
  regime_by_t <- rep("Nominal", n)
  for (r in seq_len(nrow(schedule))) {
    idx <- which(t >= schedule$t_start[r] & t < schedule$t_end[r])
    regime_by_t[idx] <- schedule$regime[r]
  }
  
  #Copy physics signals to modify per regime for "approx" side before residuals
  x_pos_A_fault <- xpos_A
  yA_A_fault    <- yA_A
  yB_A_fault    <- yB_A
  
  #Inject effects on the ARX side to create residual patterns
  for (k in seq_len(n)) {
    reg <- regime_by_t[k]
    if (reg == "SensorNoisy") {
      yA_A_fault[k] <- yA_A_fault[k] + rnorm(1, 0, 0.06)  
      yB_A_fault[k] <- yB_A_fault[k] + rnorm(1, 0, 0.05)
      
      
    } else if (reg == "DynamicsOff") {
      
      x_pos_A_fault[k] <- 0.8 * x_pos_A_fault[k]
    } else if (reg == "Drift") {
   
      yA_A_fault[k] <- yA_A_fault[k] + 0.0007 * (k-1)
    }
  }
  y_fused_A_fault <- params$alpha*yA_A_fault + (1 - params$alpha)*yB_A_fault
  z_A_fault <- params$w1*x_pos_A_fault + params$w2*c(0, diff(x_pos_A_fault)/params$dt)
  
  
  E <- tibble::tibble(
    t   = t,
    e0  = u_act_A          - u_act_P,
    e1  = x_pos_A_fault    - x_pos_P,
    e2a = yA_A_fault       - yA_P,
    e2b = yB_A_fault       - yB_P,
    e3  = y_fused_A_fault  - y_fused_P,
    e4  = z_A_fault        - z_P
  )
  
  list(
    E = E, t = t, regime_t = regime_by_t,
    u_act_A = u_act_A, xpos_A = x_pos_A_fault, yA_A = yA_A_fault, yB_A = yB_A_fault,
    y_fused_A = y_fused_A_fault, z_A = z_A_fault,
    u_act_P = u_act_P, x_pos_P = x_pos_P, yA_P = yA_P, yB_P = yB_P,
    y_fused_P = y_fused_P, z_P = z_P
  )
}


sim <- simulate_with_faults(params, schedule)
E    <- sim$E
t    <- sim$t
gt   <- tibble::tibble(t = t, regime = sim$regime_t)

# Extract individual signals for use in TargetedIntervention.R
u_act_A   <- sim$u_act_A
xpos_A    <- sim$xpos_A
yA_A      <- sim$yA_A
yB_A      <- sim$yB_A
y_fused_A <- sim$y_fused_A
z_A       <- sim$z_A
u_act_P   <- sim$u_act_P
x_pos_P   <- sim$x_pos_P
yA_P      <- sim$yA_P
yB_P      <- sim$yB_P
y_fused_P <- sim$y_fused_P
z_P       <- sim$z_P


feat <- E %>%
  dplyr::transmute(
    t,
    e0 = e0, e1 = e1, e2a = e2a, e2b = e2b, e3 = e3, e4 = e4,
    e_abs_sum = abs(e0)+abs(e1)+abs(e2a)+abs(e2b)+abs(e3)+abs(e4),
    e_rms     = sqrt((e0^2+e1^2+e2a^2+e2b^2+e3^2+e4^2)/6),
    e1_diff   = c(0, diff(e1)),
    e3_diff   = c(0, diff(e3))
  ) %>%
  #z-scoring features
  dplyr::mutate(dplyr::across(-t, ~ as.numeric(scale(.x))))

#fitting HMMs here:

Ez <- feat %>% dplyr::select(-t) %>% as.data.frame()
fams <- rep(list(gaussian()), ncol(Ez))
mod  <- depmix(response = lapply(names(Ez), function(nm) as.formula(paste0(nm, "~ 1"))),
               data = Ez, nstates = params$K, family = fams)
fm   <- fit(mod, verbose = FALSE)
post <- tibble::as_tibble(posterior(fm, type = "smoothing"))

#Extract/normalize state posteriors γ_t(k)
s_cols <- grep("^(S|V)\\d+$", names(post), value = TRUE)
s_cols <- s_cols[order(as.integer(sub("^[SV]", "", s_cols)))]
gamma_long <- post %>%
  dplyr::mutate(time = dplyr::row_number()) %>%
  tidyr::pivot_longer(tidyselect::all_of(s_cols),
                      names_to = "state_col", values_to = "prob_raw") %>%
  dplyr::mutate(state = as.integer(sub("^[SV]", "", state_col))) %>%
  dplyr::group_by(time) %>%
  dplyr::mutate(prob = prob_raw / pmax(sum(prob_raw), 1e-12)) %>%
  dplyr::ungroup() %>%
  dplyr::select(time, state, prob)

#Viterbi path and contingency
vpath <- viterbi(fm)$state
tab   <- table(State = vpath, Regime = gt$regime)   


C <- as.matrix(tab)

nr <- nrow(C); nc <- ncol(C)
if (nr > nc) C <- cbind(C, matrix(0, nrow = nr, ncol = nr - nc))
if (nc > nr) C <- rbind(C, matrix(0, nrow = nc - nr, ncol = nc))


M <- max(C)
cost <- M - C             

###The hungarian part
perm_full <- clue::solve_LSAP(cost)   


state_names  <- rownames(C)
regime_names <- c(colnames(tab), rep("<pad>", ncol(C) - ncol(tab)))


map_regimes <- regime_names[perm_full][seq_len(nrow(tab))]
names(map_regimes) <- rownames(tab)  


gamma_mapped <- gamma_long %>%
  dplyr::mutate(regime = map_regimes[as.character(state)]) %>%
  dplyr::mutate(t = t[time])


p_gamma <- ggplot(gamma_mapped, aes(t, prob, color = regime)) +
  geom_rect(data = schedule,
            aes(xmin = t_start, xmax = t_end, ymin = -Inf, ymax = Inf, fill = regime),
            inherit.aes = FALSE, alpha = 0.10, color = NA) +
  geom_vline(data = schedule, aes(xintercept = t_start),
             linetype = "dashed", color = "gray40", alpha = 0.6) +
  geom_line(size=1) +
  labs(title = "HMM posteriors mapped to ground-truth regimes",
       x = "Time (s)", y = expression(gamma[t](regime)), color = "Regime", fill = "Regime") +
  theme_minimal()+
  scale_fill_manual(values=c('darkorange', 'darkred', 'navy', 'darkgreen'))+
  scale_color_manual(values=c('darkorange', 'darkred', 'navy', 'darkgreen'))
print(p_gamma)


feat_w_regime <- feat %>% dplyr::left_join(gt, by = "t")

diag_box <- feat_w_regime %>%
  tidyr::pivot_longer(cols = -c(t, regime), names_to = "feature", values_to = "val") %>%
  ggplot(aes(regime, val)) +
  geom_boxplot(outlier.alpha = 0.2, size=1) +
  facet_wrap(~ feature, scales = "free_y") +
  labs(title = "Feature separability by ground-truth regime",
       x = "Regime", y = "Feature value (z-scored)") +
  theme_minimal()
print(diag_box)















