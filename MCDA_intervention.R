#' Authors@R: person("Annice", "Najafi", email = "annicenajafi27@gmail.com")
#' Fall 2025
#' 
#' 
#' Code for aMCDA based action plan for error propagation mitigation, figure 11 in the article


#Load relevant libraries
library(tidyverse)
library(RMCDA)   
library(zoo)     



params <- list(
  dt = 0.02,         
  roll_sec = 1.0,      #rolling window length for noise/drift proxies
  
  #TOPSIS weights
  weights = c(
    KPI_gain    = 0.35,
    Module_gain = 0.25,
    Drift_gain  = 0.15,
    Stability   = 0.10,
    Cost        = 0.15
  ),
  
  #Gate NoAction only when truly nominal/quiet
  gate_quantile = 0.20,  #use 20th percentile of kpi_n/agg_n as quiet
  gate_pNom     = 0.75,  #require nominal prob above this to allow NoAction
  
  #Anti-chatter switching
  min_dwell_sec = 0.50,
  switch_margin = 0.02
)


stopifnot(exists("E", inherits = TRUE))
E <- as_tibble(get("E", inherits = TRUE))

req_cols <- c("e0","e1","e2a","e2b","e3","e4")
miss <- setdiff(req_cols, names(E))
if (length(miss) > 0) {
  stop("E missing required columns: ", paste(miss, collapse = ", "),
       "\nAvailable: ", paste(names(E), collapse = ", "))
}


E <- E %>% mutate(time = row_number())


if (!("t" %in% names(E))) {
  message("E has no 't' column. Creating seconds using params$dt.")
  E <- E %>% mutate(t = (time - 1) * params$dt)
}

E <- E %>% arrange(time)


robust_scale <- function(x) {
  s <- mad(x, constant = 1, na.rm = TRUE)
  if (!is.finite(s) || s < 1e-12) s <- sd(x, na.rm = TRUE)
  if (!is.finite(s) || s < 1e-12) s <- 1
  x / s
}

roll_n <- max(5L, as.integer(round(params$roll_sec / params$dt)))

feat <- E %>%
  transmute(
    time = time,
    t    = t,
    
    
    kpi_mag    = abs(e4),
    plant_mag  = abs(e1),
    sensA_mag  = abs(e2a),
    sensB_mag  = abs(e2b),
    sensor_mag = abs(e2a) + abs(e2b),
    agg_L2     = sqrt(e0^2 + e1^2 + e2a^2 + e2b^2 + e3^2 + e4^2),
    

    sens_noise = zoo::rollapply(abs(e2a) + abs(e2b), roll_n, sd, fill = NA, align = "right"),
    plant_var  = zoo::rollapply(e1, roll_n, sd, fill = NA, align = "right"),
    
   
    drift_e2a  = zoo::rollapply(e2a, roll_n, mean, fill = NA, align = "right"),
    drift_e0   = zoo::rollapply(e0,  roll_n, mean, fill = NA, align = "right"),
    drift_e1   = zoo::rollapply(e1,  roll_n, mean, fill = NA, align = "right")
  ) %>%
  mutate(
    sens_noise = replace_na(sens_noise, 0),
    plant_var  = replace_na(plant_var, 0),
    drift_mag  = replace_na(abs(drift_e2a) + abs(drift_e0) + abs(drift_e1), 0),
    

    kpi_n    = robust_scale(kpi_mag),
    plant_n  = robust_scale(plant_mag),
    sensor_n = robust_scale(sensor_mag),
    noise_n  = robust_scale(sens_noise),
    drift_n  = robust_scale(drift_mag),
    agg_n    = robust_scale(agg_L2),
    
    A_dom = sensA_mag >= sensB_mag
  )


pcols <- c("p_Nominal","p_SensorNoisy","p_DynamicsOff","p_Drift")

have_gamma <- exists("gamma_mapped", inherits = TRUE)

if (have_gamma) {
  gm <- as_tibble(get("gamma_mapped", inherits = TRUE))
  
  
  if (!("time" %in% names(gm))) {

    gm <- gm %>% mutate(time = row_number())
  }
  
  stopifnot(all(c("time","regime","prob") %in% names(gm)))
  
  gm_wide <- gm %>%
    group_by(time, regime) %>%
    summarise(prob = sum(prob), .groups = "drop") %>%
    group_by(time) %>%
    mutate(prob = prob / pmax(sum(prob), 1e-12)) %>%
    ungroup() %>%
    pivot_wider(names_from = regime, values_from = prob, values_fill = 0) %>%
    rename_with(~paste0("p_", .x), -time)
  
  feat <- feat %>% left_join(gm_wide, by = "time")
}


if (!have_gamma) {
  proxy <- feat %>%
    transmute(
      time = time,
      s_nom   = exp(-pmax(0, agg_n)),
      s_sens  = pmax(0, noise_n) * pmax(0, sensor_n),
      s_dyn   = pmax(0, plant_n) * pmax(0, plant_var),
      s_drift = pmax(0, drift_n)
    ) %>%
    rowwise() %>%
    mutate(Z = s_nom + s_sens + s_dyn + s_drift + 1e-12,
           p_Nominal     = s_nom / Z,
           p_SensorNoisy = s_sens / Z,
           p_DynamicsOff = s_dyn / Z,
           p_Drift       = s_drift / Z) %>%
    ungroup() %>%
    select(time, all_of(pcols))
  
  feat <- feat %>% left_join(proxy, by = "time")
}


feat <- feat %>%
  mutate(across(all_of(pcols), ~ replace_na(.x, 0))) %>%
  rowwise() %>%
  mutate(ps = sum(c_across(all_of(pcols))),
         p_Nominal     = ifelse(ps < 1e-12, 1, p_Nominal     / ps),
         p_SensorNoisy = ifelse(ps < 1e-12, 0, p_SensorNoisy / ps),
         p_DynamicsOff = ifelse(ps < 1e-12, 0, p_DynamicsOff / ps),
         p_Drift       = ifelse(ps < 1e-12, 0, p_Drift       / ps)) %>%
  ungroup() %>%
  select(-ps)


actions <- tibble(
  action = c("DownweightSensorA","DownweightSensorB","IncreaseFiltering","ReidentifyPlant","BiasCorrect"),
  Cost = c(0.10, 0.10, 0.15, 0.40, 0.25),
  StabilityPenalty = c(0.08, 0.08, 0.05, 0.15, 0.12)
)


w <- params$weights
w <- w / sum(w)


action_scores_at_row <- function(f) {
  pN <- f$p_Nominal
  pS <- f$p_SensorNoisy
  pD <- f$p_DynamicsOff
  pR <- f$p_Drift
  
  kpi   <- f$kpi_n
  sens  <- f$sensor_n
  plant <- f$plant_n
  drift <- f$drift_n
  noise <- f$noise_n
  
  A_dom <- isTRUE(f$A_dom)
  
 
  tibble(
    action = actions$action,
    
    KPI_gain = c(
  
      pS * kpi * sens * ifelse(A_dom, 1.2, 0.6),
      pS * kpi * sens * ifelse(!A_dom, 1.2, 0.6),
      pS * kpi * noise,
      pD * kpi * plant,
      (pR + 0.5*pD) * kpi * (drift + 0.3*plant)
    ),
    
    Module_gain = c(
      pS * sens * ifelse(A_dom, 1.2, 0.6),
      pS * sens * ifelse(!A_dom, 1.2, 0.6),
      pS * sens * noise,
      pD * plant,
      (pR + 0.5*pD) * (drift + 0.3*plant)
    ),
    
    Drift_gain = c(
      pS * drift * ifelse(A_dom, 1.1, 0.7),
      pS * drift * ifelse(!A_dom, 1.1, 0.7),
      pS * drift * 0.6,
      pD * drift * 0.4,
      (pR + 0.2*pS) * drift
    )
  ) %>%
    left_join(actions, by = "action") %>%
    mutate(
      Stability = -StabilityPenalty
    ) %>%
    select(action, KPI_gain, Module_gain, Drift_gain, Stability, Cost)
}

rank_topsis <- function(score_tbl) {
  A <- score_tbl %>%
    column_to_rownames("action") %>%
    as.matrix()
  
  A[, "Cost"] <- -A[, "Cost"]
  
  s <- suppressWarnings(RMCDA::apply.TOPSIS(A, w = as.numeric(w[colnames(A)])))
  s <- as.numeric(s)
  s[!is.finite(s)] <- 0  
  
  tibble(action = rownames(A), score = s) %>%
    arrange(desc(score)) %>%
    mutate(rank = row_number())
}


thr_kpi <- as.numeric(quantile(feat$kpi_n, params$gate_quantile, na.rm = TRUE))
thr_agg <- as.numeric(quantile(feat$agg_n, params$gate_quantile, na.rm = TRUE))
thr_pN  <- params$gate_pNom

choose_action_at_row <- function(f) {
  if (f$kpi_n < thr_kpi && f$agg_n < thr_agg && f$p_Nominal > thr_pN) {
    return(tibble(action = "NoAction", score = NA_real_, rank = 1L))
  }
  sc <- action_scores_at_row(f)
  rank_topsis(sc)
}


rankings <- feat %>%
  group_by(time) %>%
  group_modify(~{
    f <- .x[1, ]
    out <- choose_action_at_row(f)
    out %>% mutate(t = f$t)
  }) %>%
  ungroup()

best_raw <- rankings %>% filter(rank == 1) %>% select(time, t, action, score)


min_dwell_n <- max(1L, as.integer(round(params$min_dwell_sec / params$dt)))
margin <- params$switch_margin

apply_hysteresis <- function(best_tbl, rankings_tbl) {
  best_tbl <- best_tbl %>% arrange(time)
  
  chosen <- character(nrow(best_tbl))
  chosen[1] <- best_tbl$action[1]
  hold <- 0L
  
  for (i in 2:nrow(best_tbl)) {
    if (hold < min_dwell_n) {
      chosen[i] <- chosen[i-1]
      hold <- hold + 1L
      next
    }
    
    time_now <- best_tbl$time[i]
    new_action <- best_tbl$action[i]
    new_score  <- best_tbl$score[i]
    if (!is.finite(new_score)) new_score <- -Inf
    
    prev_action <- chosen[i-1]
    prev_score <- rankings_tbl %>%
      filter(time == time_now, action == prev_action) %>%
      pull(score)
    prev_score <- if (length(prev_score) == 0) -Inf else prev_score[1]
    if (!is.finite(prev_score)) prev_score <- -Inf
    
    if (is.finite(new_score) && (new_score - prev_score) > margin) {
      chosen[i] <- new_action
      hold <- 0L
    } else {
      chosen[i] <- prev_action
      hold <- hold + 1L
    }
  }
  
  best_tbl %>% mutate(action_hyst = chosen)
}

best <- apply_hysteresis(best_raw, rankings)


p <- ggplot(best, aes(x = t, y = action_hyst)) +
  geom_point(size = 1.5, color="#3B4953", stroke=1.5) +
  labs(
    title = "TOPSIS-selected action over time",
    x = "Time (s)", y = "Selected action"
  ) +
  theme_minimal(base_size = 18) +
  theme(#panel.grid.major = element_blank(),
        #panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color = "black"),
        plot.title = element_text(hjust = 0.5, size = 20),
        strip.text = element_text(size = 12))
print(p)

cat("\nGate thresholds:\n")
cat("  thr_kpi =", thr_kpi, "\n")
cat("  thr_agg =", thr_agg, "\n")
cat("  thr_pNom=", thr_pN, "\n")

cat("\nAction counts (after hysteresis):\n")
print(best %>% count(action_hyst, sort = TRUE))

cat("\nNA count in p_Nominal: ", sum(is.na(feat$p_Nominal)), "\n")
cat("Range p_Nominal: ", paste(range(feat$p_Nominal, na.rm=TRUE), collapse=" .. "), "\n")


