#' Authors@R: person("Annice", "Najafi", email = "annicenajafi27@gmail.com")
#' Fall 2025
#' Script to run robusness analysis in the discussion section, figure 12

library(tidyverse)
library(depmixS4)



BASE <- list(
  dt   = 0.02, T = 20, seed = 11L,

  tau_cmd = 0.15, cmd_proc_std = 0.005,

  m = 1.0, c = 0.4, k = 4.0, plant_proc_std = 0.04,
 
  tau_sa = 0.08, sa_meas_std = 0.03, sa_bias = 0.0,
  tau_sb = 0.25, sb_meas_std = 0.02, sb_bias = 0.01,

  alpha = 0.6, w1 = 1.0, w2 = 0.2,
 
  train_frac = 0.5,
  na_cmd = 1L, nb_cmd = 1L,
  na_plant = 2L, nb_plant = 2L,
  na_sens = 1L, nb_sens = 1L,

  K = 4
)



make_u_cmd <- function(n, amp=1, dt=0.02){
  t <- (0:(n-1))*dt
  u <- 0.6*sin(2*pi*0.25*t) + 0.4*sign(sin(2*pi*0.05*t))
  step_every <- as.integer(2.5/dt)
  for(i in seq(1,n,by=step_every)) u[i:n] <- u[i:n] + 0.2*rnorm(1)
  amp*u
}
lpf_first_order <- function(u, tau, dt, proc_std=0){
  n <- length(u); y <- numeric(n); a <- exp(-dt/max(tau,1e-6))
  for(k in 2:n) y[k] <- a*y[k-1] + (1-a)*u[k-1] + proc_std*rnorm(1)
  y
}
plant_msd <- function(u, dt, m, c, k, proc_std=0){
  n <- length(u); x <- numeric(n); v <- numeric(n)
  for(i in 2:n){
    acc <- (u[i-1] - c*v[i-1] - k*x[i-1])/m + proc_std*rnorm(1)
    v[i] <- v[i-1] + dt*acc
    x[i] <- x[i-1] + dt*v[i-1]
  }
  list(pos=x, vel=v)
}
sensor_first_order <- function(xpos, dt, tau, meas_std=0, bias=0){
  y <- lpf_first_order(xpos, tau, dt, 0)
  y + bias + rnorm(length(y), sd=meas_std)
}

fit_arx <- function(y,u,na,nb){
  stopifnot(length(y)==length(u))
  kmax <- length(y); p <- max(na,nb); if(kmax<=p) stop("ARX: series too short")
  rows <- vector("list", kmax-p); Tvec <- numeric(kmax-p)
  for(k in (p+1):kmax){
    phi <- c(sapply(1:na, \(i) y[k-i]), sapply(1:nb, \(j) u[k-j]))
    rows[[k-p]] <- phi; Tvec[k-p] <- y[k]
  }
  Phi <- do.call(rbind, rows)
  theta <- as.numeric(solve(t(Phi)%*%Phi, t(Phi)%*%Tvec))
  list(theta=theta, na=na, nb=nb)
}
sim_arx <- function(u, y_seed, model){
  na <- model$na; nb <- model$nb; th <- model$theta
  n <- length(u); yhat <- numeric(n); stopifnot(length(y_seed)>=na)
  yhat[1:na] <- tail(y_seed, na); p <- max(na,nb)
  for(k in (p+1):n){
    phi <- c(sapply(1:na,\(i) yhat[k-i]), sapply(1:nb,\(j) u[k-j]))
    yhat[k] <- sum(th*phi)
  }
  yhat
}
ema_smooth <- function(x, a){ out <- numeric(length(x)); out[1]=x[1]; for(i in 2:length(x)) out[i] <- a*x[i]+(1-a)*out[i-1]; out }
rate_limit <- function(x, max_step){ out <- numeric(length(x)); out[1]=x[1]; for(i in 2:length(x)){ step <- x[i]-out[i-1]; step <- pmin(pmax(step,-max_step),max_step); out[i] <- out[i-1]+step }; out }

get_gamma <- function(fm){
  tb <- tibble::as_tibble(posterior(fm, type="smoothing"))
  s_cols <- grep("^(S|V)\\d+$", names(tb), value=TRUE)
  s_cols <- s_cols[order(as.integer(sub("^[SV]","",s_cols)))]
  tb %>%
    dplyr::mutate(time=dplyr::row_number()) %>%
    tidyr::pivot_longer(dplyr::all_of(s_cols), names_to="state_col", values_to="prob_raw") %>%
    dplyr::mutate(k = as.integer(sub("^[SV]","", state_col))) %>%
    dplyr::group_by(time) %>% dplyr::mutate(prob = prob_raw / pmax(sum(prob_raw),1e-12)) %>%
    dplyr::ungroup() %>% dplyr::select(time,k,prob)
}


run_one <- function(
    base = BASE,
    
    alpha_fuse = BASE$alpha,
    sa_mult    = 1.0, sb_mult = 1.0,
    na_plant   = BASE$na_plant, nb_plant = BASE$nb_plant,
    K          = BASE$K,
    kappa = 0.6, ema = 0.15, rmax = 5e-4
){
  set.seed(base$seed)
  n <- as.integer(base$T/base$dt); t <- (0:(n-1))*base$dt
  u_cmd <- make_u_cmd(n, dt=base$dt)
  

  u_act_P <- lpf_first_order(u_cmd, base$tau_cmd, base$dt, base$cmd_proc_std)
  pl <- plant_msd(u_act_P, base$dt, base$m, base$c, base$k, base$plant_proc_std)
  x_pos_P <- pl$pos; x_vel_P <- pl$vel
  yA_P <- sensor_first_order(x_pos_P, base$dt, base$tau_sa, base$sa_meas_std*sa_mult, base$sa_bias)
  yB_P <- sensor_first_order(x_pos_P, base$dt, base$tau_sb, base$sb_meas_std*sb_mult, base$sb_bias)
  y_fused_P <- alpha_fuse*yA_P + (1-alpha_fuse)*yB_P
  z_P <- base$w1*x_pos_P + base$w2*x_vel_P
  
 
  split <- as.integer(base$train_frac*n)
  arx_cmd   <- fit_arx(u_act_P[1:split], u_cmd[1:split], base$na_cmd, base$nb_cmd)
  arx_plant <- fit_arx(x_pos_P[1:split], u_act_P[1:split], na_plant, nb_plant)
  arx_sa    <- fit_arx(yA_P[1:split], x_pos_P[1:split], base$na_sens, base$nb_sens)
  arx_sb    <- fit_arx(yB_P[1:split], x_pos_P[1:split], base$na_sens, base$nb_sens)
  
  u_act_A <- sim_arx(u_cmd, y_seed=u_act_P[1:base$na_cmd], arx_cmd)
  xpos_A  <- sim_arx(u_act_A, y_seed=x_pos_P[1:na_plant], arx_plant)
  xvel_A  <- c(0, diff(xpos_A)/base$dt)
  yA_A    <- sim_arx(xpos_A, y_seed=yA_P[1:base$na_sens], arx_sa)
  yB_A    <- sim_arx(xpos_A, y_seed=yB_P[1:base$na_sens], arx_sb)
  y_fused_A <- alpha_fuse*yA_A + (1-alpha_fuse)*yB_A
  z_A <- base$w1*xpos_A + base$w2*xvel_A
  

  E <- tibble(
    t=t,
    e0 = u_act_A - u_act_P,
    e1 = xpos_A  - x_pos_P,
    e2a = yA_A - yA_P,
    e2b = yB_A - yB_P,
    e3 = y_fused_A - y_fused_P,
    e4 = z_A - z_P
  )
  

  Ez <- E %>% dplyr::mutate(dplyr::across(starts_with("e"), ~ as.numeric(scale(.x))))
  mod <- depmix(
    response = list(e0~1,e1~1,e2a~1,e2b~1,e3~1,e4~1),
    data = as.data.frame(Ez),
    nstates = K,
    family = replicate(6, gaussian(), simplify = FALSE)
  )
  fm <- fit(mod, verbose = FALSE)
  
 
  vpath <- viterbi(fm)$state
  st_means <- Ez %>%
    dplyr::mutate(state=vpath) %>%
    dplyr::group_by(state) %>%
    dplyr::summarise(mu_e1_z = mean(e1), .groups="drop") %>%
    dplyr::arrange(state)
  

  gamma <- get_gamma(fm)
  mu_k <- st_means$mu_e1_z
  bias_z <- gamma %>%
    dplyr::mutate(w = mu_k[k]) %>%
    dplyr::group_by(time) %>%
    dplyr::summarise(bias_z = sum(prob*w), .groups="drop") %>% dplyr::pull(bias_z)
  
 
  e1_sd <- sd(E$e1); e1_mean <- mean(E$e1)  
  bias_m <- kappa * (bias_z*e1_sd + e1_mean)
  bias_sm <- ema_smooth(bias_m, ema)
  bias_rl <- rate_limit(bias_sm, rmax)
  
  xpos_C <- xpos_A - bias_rl
  E_corr <- E %>% dplyr::mutate(e1 = xpos_C - (x_pos_P))  
  
 
  L2  <- function(df) sqrt(df$e0^2+df$e1^2+df$e2a^2+df$e2b^2+df$e3^2+df$e4^2)
  rmse <- function(x) sqrt(mean(x^2)); mae <- function(x) mean(abs(x))
  
  tibble(
    K=K, alpha_fuse=alpha_fuse, sa_mult=sa_mult, sb_mult=sb_mult,
    na_plant=na_plant, nb_plant=nb_plant,
    kappa=kappa, ema=ema, rmax=rmax,
    RMSE_e1_orig = rmse(E$e1), RMSE_e1_corr = rmse(E_corr$e1),
    MAE_e1_orig  = mae(E$e1),  MAE_e1_corr  = mae(E_corr$e1),
    L2_orig = mean(L2(E)),     L2_corr     = mean(L2(E_corr)),
    dRMSE_e1_pct = 100*(RMSE_e1_orig - RMSE_e1_corr)/RMSE_e1_orig,
    dMAE_e1_pct  = 100*(MAE_e1_orig  - MAE_e1_corr)/MAE_e1_orig,
    dL2_pct      = 100*(L2_orig - L2_corr)/L2_orig
  )
}


kappas <- c(0.0, 0.2, 0.4, 0.6, 0.8)
emas   <- c(0.05, 0.10, 0.15, 0.25, 0.35)
rmaxs  <- c(1e-4, 2e-4, 5e-4, 1e-3)



alphas <- c(0.2, 0.4, 0.6, 0.8)
noise_levels <- c(0.5, 1.0, 1.5)

sa2_grid <- expand.grid(alpha_fuse = alphas,
                        sa_mult = noise_levels,
                        sb_mult = noise_levels)

sa2 <- purrr::pmap_dfr(sa2_grid, \(alpha_fuse, sa_mult, sb_mult){
  run_one(alpha_fuse=alpha_fuse, sa_mult=sa_mult, sb_mult=sb_mult)
})

p_sa2 <- sa2 %>%
  mutate(sa_lab = paste0("Sensor A ×", sa_mult),
         sb_lab = paste0("Sensor B ×", sb_mult)) %>%
  ggplot(aes(alpha_fuse, dRMSE_e1_pct, group=alpha_fuse)) +
  geom_line(color="#547792", size=2) + geom_point(color="#547792", shape=8, stroke=1.5) +
  facet_grid(sa_lab ~ sb_lab) +
  scale_x_continuous(breaks = alphas) +
  labs(title="Robustness to fusion and noise",
       x=TeX("Fusion weight $\alpha$ (A vs B)"),
       y=TeX("$\\Delta$ RMSE_e1 (%)")) +
  theme_minimal(base_size = 16)+
  theme(panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color = "black"),
        plot.title = element_text(hjust = 0.5, size = 20),
        strip.text = element_text(size = 12))

print(p_sa2)







