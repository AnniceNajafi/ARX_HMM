#' Authors@R: person("Annice", "Najafi", email = "annicenajafi27@gmail.com")
#' Fall 2025
#' This script performs targeted interventions using HMM results and plots the residuals
#' before and after HMM correction for demo
#' NOTE: Run AccuracyTesting.R first to create required variables

#Configuration for targeted intervention (uses params from AccuracyTesting.R)


CFG <- list(
  kappa     = 0.8,      #correction gain
  ema_alpha = 0.15,     #EMA smoothing factor
  r_max     = 0.05,     #rate limit max step
  dt        = params$dt,
  alpha     = params$alpha,
  w1        = params$w1,
  w2        = params$w2
)


ema_smooth <- function(x, a) {
  out <- numeric(length(x))
  out[1] <- x[1]
  for (i in 2:length(x)) out[i] <- a * x[i] + (1 - a) * out[i - 1]
  out
}

rate_limit <- function(x, max_step) {
  out <- numeric(length(x))
  out[1] <- x[1]
  for (i in 2:length(x)) {
    step <- x[i] - out[i - 1]
    step <- pmin(pmax(step, -max_step), max_step)
    out[i] <- out[i - 1] + step
  }
  out
}

K <- max(as.integer(gamma_long$state))


mu_state <- tibble::tibble(state = 1:K)
for (ch in c("e0", "e1", "e2a", "e2b", "e3", "e4")) {
  mu_state[[ch]] <- sapply(1:K, function(k) {
    gamma_k <- gamma_long %>%
      dplyr::filter(as.integer(state) == k) %>%
      dplyr::arrange(time) %>%
      dplyr::pull(prob)
    if (sum(gamma_k) < 1e-12) return(0)
    sum(gamma_k * E[[ch]]) / sum(gamma_k)
  })
}


bias_z_from_gamma <- function(mu_k_vec, gamma_long) {
  #mu_k_vec: length-K vector of per-state means (z-units) for one residual
  gamma_long %>%
    dplyr::mutate(w = mu_k_vec[as.integer(state)]) %>%
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
  CFG$kappa * bias_z[[nm]]  
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


E_orig <- E %>% dplyr::transmute(t, e0, e1, e2a, e2b, e3, e4)

E_corr <- tibble::tibble(
  t   = t,
  e0  = u_act_C   - u_act_P,
  e1  = xpos_C    - x_pos_P,
  e2a = yA_C      - yA_P,
  e2b = yB_C      - yB_P,
  e3  = y_fused_C - y_fused_P,
  e4  = z_C       - z_P
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
  E_corr %>% tidyr::pivot_longer(c(e0, e1, e2a, e2b, e3, e4), names_to = "channel", values_to = "err") %>% dplyr::mutate(set = "Corrected"),
  E_orig %>% tidyr::pivot_longer(c(e0, e1, e2a, e2b, e3, e4), names_to = "channel", values_to = "err") %>% dplyr::mutate(set = "Original")
)

lab_map <- c(
  e0="M0 (Actuator)", e1="M1 (Plant pos)", e2a="M2a (Sensor A)",
  e2b="M2b (Sensor B)", e3="M3 (Fusion)", e4="M4 (KPI)"
)


ts_long$set <- factor(ts_long$set, levels = c("Original", "Corrected"))

col_original  <- "#A5B68D"   
col_corrected <- "#273F4F"  

p_ts <- ggplot(ts_long, aes(x = t, y = err, color = set)) +
  geom_line(linewidth = 1.2, alpha = 0.85) +
  facet_wrap(~ factor(channel, levels = ch_names, labels = lab_map[ch_names]),
             scales = "free_y", ncol = 2) +
  scale_color_manual(values = c("Original" = col_original, "Corrected" = col_corrected)) +
  labs(title = "Per-Module Residuals: Original vs HMM-Corrected",
       x = "Time (s)", y = "Residual", color = NULL) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.major = element_line(color = "gray92", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    axis.line = element_line(color = "gray30", linewidth = 0.4),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold", color = "gray20"),
    strip.text = element_text(size = 11, face = "bold", color = "gray30"),
    strip.background = element_rect(fill = "gray95", color = NA),
    legend.position = "bottom",
    legend.text = element_text(size = 11),
    axis.title = element_text(size = 11, color = "gray30"),
    axis.text = element_text(size = 9, color = "gray40")
  )


dens_long <- ts_long %>%
  dplyr::group_by(channel, set) %>%
  dplyr::mutate(channel_lab = lab_map[channel[1]]) %>% dplyr::ungroup()


col_original_dark  <- "#6B7A5E"
col_corrected_dark <- "#1A2A34"

p_dens <- ggplot(dens_long, aes(x = err, fill = set, color = set)) +
  geom_density(alpha = 0.5, linewidth = 1.2) +
  facet_wrap(~ factor(channel, levels = ch_names, labels = lab_map[ch_names]), scales = "free", ncol = 2) +
  scale_fill_manual(values = c("Original" = col_original, "Corrected" = col_corrected)) +
  scale_color_manual(values = c("Original" = col_original_dark, "Corrected" = col_corrected_dark)) +
  labs(title = "Residual Distributions by Module",
       x = "Residual", y = "Density", fill = NULL, color = NULL) +
  theme_minimal(base_size = 13) +
  guides(color = "none") +
  theme(
    panel.grid.major = element_line(color = "gray92", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    axis.line = element_line(color = "gray30", linewidth = 0.4),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold", color = "gray20"),
    strip.text = element_text(size = 11, face = "bold", color = "gray30"),
    strip.background = element_rect(fill = "gray95", color = NA),
    legend.position = "bottom",
    legend.text = element_text(size = 11),
    axis.title = element_text(size = 11, color = "gray30"),
    axis.text = element_text(size = 9, color = "gray40")
  )


agg_long <- bind_rows(
  E_orig %>% dplyr::mutate(set = "Original"),
  E_corr %>% dplyr::mutate(set = "Corrected")
) %>%
  dplyr::rowwise() %>%
  dplyr::mutate(L2 = sqrt(sum(c(e0,e1,e2a,e2b,e3,e4)^2))) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(set = factor(set, levels = c("Original", "Corrected")))

p_agg <- ggplot(agg_long, aes(t, L2, color = set)) +
  geom_line(linewidth = 1.2, alpha = 0.85) +
  scale_color_manual(values = c("Original" = col_original, "Corrected" = col_corrected)) +
  geom_hline(yintercept = mean((agg_long %>% dplyr::filter(set == "Original"))$L2),
             color = col_original, linetype = "dashed", linewidth = 0.8, alpha = 0.8) +
  geom_hline(yintercept = mean((agg_long %>% dplyr::filter(set == "Corrected"))$L2),
             color = col_corrected, linetype = "dashed", linewidth = 0.8, alpha = 0.8) +
  labs(title = "Aggregate Residual Magnitude Across Modules",
       x = "Time (s)", y = "L2 Magnitude", color = NULL) +
  theme_minimal(base_size = 13) +
  theme(
    panel.grid.major = element_line(color = "gray92", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.background = element_rect(fill = "white", color = NA),
    plot.background = element_rect(fill = "white", color = NA),
    axis.line = element_line(color = "gray30", linewidth = 0.4),
    plot.title = element_text(hjust = 0.5, size = 16, face = "bold", color = "gray20"),
    legend.position = "bottom",
    legend.text = element_text(size = 11),
    axis.title = element_text(size = 11, color = "gray30"),
    axis.text = element_text(size = 9, color = "gray40")
  )


p_ts
p_dens
p_agg
