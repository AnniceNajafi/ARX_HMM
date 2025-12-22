#' Authors@R: person("Annice", "Najafi", email = "annicenajafi27@gmail.com")
#' Fall 2025
#' This script performs targeted interventions using HMM results and plots the residuals
#' before and after HMM correction for demo

mu_state <- Ez %>%
  dplyr::mutate(state = vpath) %>%
  dplyr::group_by(state) %>%
  dplyr::summarise(dplyr::across(e0:e4, mean), .groups = "drop") %>%
  dplyr::arrange(state)                    # state = 1..K
K <- nrow(mu_state)


bias_z_from_gamma <- function(mu_k_vec, gamma_long) {
  #mu_k_vec: length-K vector of per-state means (z-units) for one residual
  gamma_long %>%
    dplyr::mutate(w = mu_k_vec[state_num]) %>%
    dplyr::group_by(time) %>%
    dplyr::summarise(bias_z = sum(prob * w), .groups = "drop") %>%
    dplyr::pull(bias_z)
}


ch_names <- c("e0","e1","e2a","e2b","e3","e4")


phys_mean <- sapply(ch_names, function(nm) mean(E[[nm]]))
phys_sd   <- sapply(ch_names, function(nm) sd(E[[nm]]))


mu_list <- list(
  e0  = mu_state$e0,
  e1  = mu_state$e1,
  e2a = mu_state$e2a,
  e2b = mu_state$e2b,
  e3  = mu_state$e3,
  e4  = mu_state$e4
)


bias_z <- lapply(ch_names, function(nm) bias_z_from_gamma(mu_list[[nm]], gamma_long))
names(bias_z) <- ch_names


bias_m <- lapply(ch_names, function(nm) {
  CFG$kappa * (bias_z[[nm]] * phys_sd[[nm]] + phys_mean[[nm]])
})
names(bias_m) <- ch_names


bias_corr <- lapply(ch_names, function(nm) {
  ema <- ema_smooth(bias_m[[nm]], CFG$ema_alpha)
  rate_limit(ema, CFG$r_max)
})
names(bias_corr) <- ch_names


u_act_C <- u_act_A - bias_corr$e0
xpos_C  <- xpos_A  - bias_corr$e1

xvel_C  <- c(0, diff(xpos_C) / CFG$dt)

yA_C    <- yA_A    - bias_corr$e2a
yB_C    <- yB_A    - bias_corr$e2b


y_fused_C <- CFG$alpha * yA_C + (1 - CFG$alpha) * yB_C
z_C       <- CFG$w1 * xpos_C + CFG$w2 * xvel_C


E_orig <- E %>% dplyr::transmute(t, e0,e1,e2a,e2b,e3,e4)

E_corr <- tibble::tibble(
  t  = E$t,
  e0 = u_act_C   - u_act_P,
  e1 = xpos_C    - x_pos_P,
  e2a= yA_C      - yA_P,
  e2b= yB_C      - yB_P,
  e3 = y_fused_C - y_fused_P,
  e4 = z_C       - z_P
)


metric_tbl <- function(Ea, Eb, label_a, label_b) {
  bind_rows(
    tibble::tibble(channel = ch_names,
                   set = label_a,
                   RMSE = sapply(ch_names, \(nm) sqrt(mean(Ea[[nm]]^2))),
                   MAE  = sapply(ch_names, \(nm) mean(abs(Ea[[nm]])))),
    tibble::tibble(channel = ch_names,
                   set = label_b,
                   RMSE = sapply(ch_names, \(nm) sqrt(mean(Eb[[nm]]^2))),
                   MAE  = sapply(ch_names, \(nm) mean(abs(Eb[[nm]]))))
  ) %>%
    tidyr::pivot_wider(names_from = set, values_from = c(RMSE, MAE)) %>%
    dplyr::mutate(
      dRMSE_pct = 100 * (RMSE_Original - RMSE_Corrected) / RMSE_Original,
      dMAE_pct  = 100 * (MAE_Original  - MAE_Corrected)  / MAE_Original
    )
}


metrics <- metric_tbl(E_orig, E_corr, "Original", "Corrected")

metrics_rounded <- metrics %>%
  dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 4)))

print(metrics_rounded, n = Inf)


ts_long <- bind_rows(
  E_orig %>% tidyr::pivot_longer(e0:e4, names_to = "channel", values_to = "err") %>% dplyr::mutate(set = "Original"),
  E_corr %>% tidyr::pivot_longer(e0:e4, names_to = "channel", values_to = "err") %>% dplyr::mutate(set = "Corrected")
)

lab_map <- c(
  e0="M0 (Actuator)", e1="M1 (Plant pos)", e2a="M2a (Sensor A)",
  e2b="M2b (Sensor B)", e3="M3 (Fusion)", e4="M4 (KPI)"
)

p_ts <- ggplot(ts_long, aes(t, err, color = set)) +
  geom_line(linewidth = 2) +
  facet_wrap(~ factor(channel, levels = ch_names, labels = lab_map[ch_names]),
             scales = "free_y", ncol = 2) +
  scale_color_manual(values = c(Original = "#561530", Corrected = "#050E3C")) +
  geom_hline(yintercept = 0, linetype="dashed", linewidth=2)+
  labs(title = "Per-module residuals under targeted intervention: Original vs HMM-corrected",
       x = "Time (s)", y = "Residual", color = NULL) +
  theme_minimal(base_size = 14) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color = "black"),
        plot.title = element_text(hjust = 0.5, size = 20),
        strip.text = element_text(size = 12))


dens_long <- ts_long %>%
  dplyr::group_by(channel, set) %>%
  dplyr::mutate(channel_lab = lab_map[channel[1]]) %>% dplyr::ungroup()

p_dens <- ggplot(dens_long, aes(x = err, fill = set)) +
  geom_density(alpha = 0.45) +
  facet_wrap(~ factor(channel, levels = ch_names, labels = lab_map[ch_names]), scales = "free", ncol = 2) +
  scale_fill_manual(values = c(Original = "#561530", Corrected = "#050E3C")) +
  labs(title = "Residual distributions by module",
       x = "Residual", y = "Density", fill = NULL) +
  theme_minimal(base_size = 14) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color = "black"),
        plot.title = element_text(hjust = 0.5, size = 20),
        strip.text = element_text(size = 12))


agg_long <- bind_rows(
  E_orig %>% dplyr::mutate(set = "Original"),
  E_corr %>% dplyr::mutate(set = "Corrected")
) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(L2 = sqrt(sum(c(e0,e1,e2a,e2b,e3,e4)^2))) %>%
  dplyr::ungroup()

p_agg <- ggplot(agg_long, aes(t, L2, color = set)) +
  geom_line(linewidth = 2) +
  scale_color_manual(values = c(Original = "#561530", Corrected = "#050E3C")) +
  labs(title = "Aggregate residual magnitude across modules",
       x = "Time (s)", y = "Magnitude", color = NULL) +
  geom_hline(yintercept = mean((agg_long%>%filter(set=="Original"))$L2), color = "#561530", linetype='dashed', linewidth=1.5)+
  geom_hline(yintercept = mean((agg_long%>%filter(set=="Corrected"))$L2), color = "#050E3C", linetype='dashed', linewidth=1.5)+
  theme_minimal(base_size = 14) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color = "black"),
        plot.title = element_text(hjust = 0.5, size = 20),
        strip.text = element_text(size = 12))


p_ts
p_dens
p_agg
