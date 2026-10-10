# -------------------------------------------------------------
# CLT Simulation: Cauchy Distribution
# -------------------------------------------------------------
library(ggplot2)
library(dplyr)
library(patchwork)

set.seed(8676)

B <- 10000L
n_values <- c(1, 5, 10, 15, 20, 25, 30, 40, 50, 75, 100, 150, 200)
output_dir <- file.path("..", "data", "cauchy")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Theoretical Parameters
theo_mu <- NA_real_
theo_sigma <- NA_real_
pop_shape <- "Continuous, symmetric, extremely heavy-tailed"

shape_statistics <- function(x) {
  center <- mean(x)
  deviations <- x - center
  m2 <- mean(deviations^2)
  if (!is.finite(m2) || m2 == 0) return(c(skewness = NA_real_, excess_kurtosis = NA_real_))
  c(skewness = mean(deviations^3) / m2^(3 / 2),
    excess_kurtosis = mean(deviations^4) / m2^2 - 3)
}

results_list <- list()
density_qq_list <- list()

for (n_val in n_values) {
  observations <- rcauchy(B * n_val, location = 0, scale = 1)
  sample_matrix <- matrix(observations, nrow = B, ncol = n_val)
  means <- rowMeans(sample_matrix)
  
  sim_mean <- mean(means)
  sim_se <- sd(means)
  stats <- shape_statistics(means)
  theo_se <- NA_real_
  
  sub_sample <- sample(means, min(5000L, B))
  shapiro_res <- shapiro.test(sub_sample)
  
  results_list[[length(results_list) + 1]] <- data.frame(
    n = n_val, theo_mean = theo_mu, sim_mean = sim_mean,
    theo_se = theo_se, sim_se = sim_se,
    skewness = unname(stats[["skewness"]]), excess_kurtosis = unname(stats[["excess_kurtosis"]]),
    shapiro_W = shapiro_res$statistic, shapiro_p = shapiro_res$p.value
  )
  
  z_means <- if (is.finite(sim_se) && sim_se > 0) (means - sim_mean) / sim_se else rep(NA_real_, B)
  density_qq_list[[length(density_qq_list) + 1]] <- data.frame(n = paste0("n = ", n_val), n_num = n_val, z = z_means)
}

summary_table <- do.call(rbind, results_list)
plot_data <- do.call(rbind, density_qq_list) %>% filter(is.finite(z))
plot_data$n <- factor(plot_data$n, levels = paste0("n = ", n_values))

write.csv(summary_table, file.path(output_dir, "cauchy_metrics.csv"), row.names = FALSE)

p_skew <- ggplot(summary_table, aes(x = n, y = skewness)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_line(color = "#E69F00", linewidth = 0.8) + geom_point(color = "#E69F00", size = 2) +
  labs(title = "Skewness (Fails to Converge)", x = "Sample Size (n)", y = "Skewness") + theme_bw()

p_kurt <- ggplot(summary_table, aes(x = n, y = excess_kurtosis)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
  geom_line(color = "#E69F00", linewidth = 0.8) + geom_point(color = "#E69F00", size = 2) +
  labs(title = "Kurtosis (Fails to Converge)", x = "Sample Size (n)", y = "Excess Kurtosis") + theme_bw()

p_dens <- ggplot(plot_data, aes(x = z)) +
  geom_density(fill = "#E69F00", alpha = 0.3) +
  stat_function(fun = dnorm, args = list(mean = 0, sd = 1), linetype = "dashed") +
  facet_wrap(~ n, scales = "free") +
  labs(title = "Empirical Density vs N(0,1)", x = "Z-score", y = "Density") + theme_bw()

p_qq <- ggplot(plot_data, aes(sample = z)) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
  stat_qq(color = "#E69F00", alpha = 0.3) +
  facet_wrap(~ n, scales = "free") +
  labs(title = "Normal Q-Q Plots", x = "Theoretical", y = "Sample") + theme_bw()

combined_plot <- (p_skew | p_kurt) / p_dens / p_qq + plot_layout(heights = c(1, 2, 2))
ggsave(file.path(output_dir, "cauchy_combined_plots.png"), plot = combined_plot, width = 14, height = 18, dpi = 200)