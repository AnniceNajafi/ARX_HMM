#' Authors@R: person("Annice", "Najafi", email = "annicenajafi27@gmail.com")
#' Fall 2025
#' **you need to run the ARX_vs_phys script first to load the module functions

library(ggplot2)



##Generate random white noise for comparison:

set.seed(42)
n <- 1000
white_noise_data <- data.frame(
  time = 1:n,
  value = rnorm(n, mean = 0, sd = 1)
)


ggplot(white_noise_data, aes(x = time, y = value)) +
  geom_line(color = "#001F3D", linewidth = 2) +
  theme_minimal() +
  labs(
    title = "Random White Noise over Time",
    x = "Time Index",
    y = "Value", tag="A"
  )+
  theme_minimal(base_size = 14) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color = "black"),
        plot.title = element_text(hjust = 0.5, size = 20),
        strip.text = element_text(size = 12))->plt.A


######Calculate ACF values and Confidence Interval
acf_obj <- acf(white_noise_data$value, plot = FALSE)
acf_df <- data.frame(lag = as.numeric(acf_obj$lag), acf = as.numeric(acf_obj$acf))


ci <- 1.96 / sqrt(nrow(white_noise_data))


ggplot(acf_df, aes(x = lag, y = acf)) +

  geom_hline(yintercept = c(ci, -ci), linetype = "dashed", color = "darkred", size=2) +
  geom_hline(yintercept = 0, color = "#001F3D") + # Baseline
  

  geom_segment(aes(xend = lag, yend = 0), color = "#001F3D", linewidth = 2) +
  

  labs(
    title = "Autocorrelation",
    x = "Lag", 
    y = "ACF", 
    tag = "B"
  ) +
  

  theme_minimal(base_size = 14) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(color = "black"),
        plot.title = element_text(hjust = 0.5, size = 20),
        strip.text = element_text(size = 12)) -> plt.B

plt.A+plt.B->plt.white.noise




library(tidyverse)
library(patchwork)


all_errors <- tibble(
  t = df$t,
  e0_Actuator = df$M0_u_act_A - df$M0_u_act_P,
  e1_Plant_Pos = df$M1_x_pos_A - df$M1_x_pos_P,
  e2a_SensorA = df$M2a_yA_A - df$M2a_yA_P,
  e2b_SensorB = df$M2b_yB_A - df$M2b_yB_P,
  e3_Fused = df$M3_y_fused_A - df$M3_y_fused_P,
  e4_Output_Z = df$M4_z_A - df$M4_z_P
) %>%
  pivot_longer(-t, names_to = "Module", values_to = "Error")



plt_E <- ggplot(all_errors, aes(x = t, y = Error)) +
  geom_line(color = "#001F3D", linewidth = 1) +
  facet_wrap(~Module, scales = "free_y", ncol = 1) +
  labs(
    title = "Time Series of Errors (All Modules)",
    x = "Time (s)", y = "Error Value", tag = "C"
  ) +
  theme_minimal(base_size = 12) +
  theme(strip.text = element_text(face = "bold"),
        panel.grid.minor = element_blank())


acf_data <- all_errors %>%
  group_by(Module) %>%
  summarise(acf_vals = list({

    a <- acf(Error, plot = FALSE, lag.max = 30)
    data.frame(lag = as.numeric(a$lag), acf = as.numeric(a$acf))
  })) %>%
  unnest(acf_vals)


n_samples <- nrow(df)
ci_bound <- 1.96 / sqrt(n_samples)

plt_F <- ggplot(acf_data, aes(x = lag, y = acf)) +
  geom_hline(yintercept = c(ci_bound, -ci_bound), linetype = "dashed", color = "darkred") +
  geom_hline(yintercept = 0, color = "black") +
  geom_segment(aes(xend = lag, yend = 0), color = "#001F3D", linewidth = 1.5) +
  facet_wrap(~Module, ncol = 1) +
  labs(
    title = "Autocorrelation of Errors (All Modules)",
    x = "Lag", y = "ACF", tag = "D"
  ) +
  theme_minimal(base_size = 12) +
  theme(strip.text = element_text(face = "bold"),
        panel.grid.minor = element_blank())


print(plt_E)
print(plt_F)


(plt.A+plt.B)/(plt_E+plt_F)







