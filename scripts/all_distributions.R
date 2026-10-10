# Central Limit Theorem: skewness and kurtosis of sample means
#
# For each n, this script generates `B` independent samples of size n,
# calculates each sample mean, and then measures the shape of the resulting
# Monte Carlo distribution of sample means.
#
# The default settings are intended for a first run. Increase `B` and extend
# `n_values` after confirming that the script runs acceptably on your machine.
library(ggplot2)
library(dplyr)

set.seed(8676) #014118676 last 4 digits

# -----------------------------
# User-configurable settings
# -----------------------------
displays_plots <- TRUE
B <- 10000L
n_values <- c(1, 5, 10, 15, 20, 25, 30, 40, 50, 75, 100, 150, 200)

# Distribution parameters
binom_size <- 5L
exp_rate <- 1
binom_prob <- 0.10

output_dir <- file.path("..", "data", "all")

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# -----------------------------
# Moment calculations
# -----------------------------
# These are the moment-ratio definitions:
#   skewness = m3 / m2^(3/2)
#   excess kurtosis = m4 / m2^2 - 3
#
# The kurtosis plotted here is excess kurtosis, so a normal distribution has
# a target value of 0. Pearson kurtosis would be excess kurtosis + 3.
shape_statistics <- function(x) {
  center <- mean(x)
  deviations <- x - center
  m2 <- mean(deviations^2)
  m3 <- mean(deviations^3)
  m4 <- mean(deviations^4)
  
  if (!is.finite(m2) || m2 == 0) {
    return(c(skewness = NA_real_, excess_kurtosis = NA_real_))
  }
  
  c(
    skewness = m3 / m2^(3 / 2),
    excess_kurtosis = m4 / m2^2 - 3
  )
}

# -----------------------------
# Random sample generators
# -----------------------------
generators <- list(
  `Standard Normal` = function(size) {
    rnorm(size, mean = 0, sd = 1)
  },
  `Six-Sided Die` = function(size) {
    sample.int(6L, size = size, replace = TRUE)
  },
  Exponential = function(size) {
    rexp(size, rate = exp_rate)
  },
  `R-Binomial` = function(size) {
    rbinom(size, size = binom_size, prob = binom_prob)
  },
  Cauchy = function(size) {
    rcauchy(size, location = 0, scale = 1)
  }
)

# -----------------------------
# Monte Carlo simulation
# -----------------------------
simulate_distribution <- function(generator, distribution_name, n_values, B) {
  results <- vector("list", length(n_values))
  
  for (i in seq_along(n_values)) {
    n <- n_values[i]
    
    # Generate B samples of size n and calculate their means. A matrix keeps
    # the simulation vectorized while limiting memory use to one n at a time.
    observations <- generator(B * n)
    sample_matrix <- matrix(observations, nrow = B, ncol = n)
    sample_means <- rowMeans(sample_matrix)
    
    stats <- shape_statistics(sample_means)
    results[[i]] <- data.frame(
      distribution = distribution_name,
      n = n,
      skewness = unname(stats[["skewness"]]),
      excess_kurtosis = unname(stats[["excess_kurtosis"]])
    )
  }
  
  do.call(rbind, results)
}

simulation_results <- do.call(
  rbind,
  lapply(names(generators), function(distribution_name) {
    message("Simulating ", distribution_name)
    simulate_distribution(
      generator = generators[[distribution_name]],
      distribution_name = distribution_name,
      n_values = n_values,
      B = B
    )
  })
)

write.csv(
  simulation_results,
  file = file.path(output_dir, "shape_statistics.csv"),
  row.names = FALSE
)

saveRDS(
  simulation_results,
  file = file.path(output_dir, "shape_statistics.rds")
)

# -------------------------------------------------------------
# Theoretical vs. Empirical Convergence Analysis Table
# -------------------------------------------------------------
message("Computing theoretical vs. empirical metrics across sample sizes...")

# 1. Theoretical moments for each generator
theoretical_params <- list(
  `Standard Normal` = list(
    pop_shape = "Continuous, symmetric, bell-shaped",
    mu = 0,
    sigma = 1
  ),
  `Six-Sided Die` = list(
    pop_shape = "Discrete uniform on {1, 2, 3, 4, 5, 6}, symmetric",
    mu = 3.5,
    sigma = sqrt(35 / 12)
  ),
  `Exponential` = list(
    pop_shape = "Continuous, heavily right-skewed, strictly positive",
    mu = 1 / exp_rate,
    sigma = 1 / exp_rate
  ),
  `R-Binomial` = list(
    pop_shape = "Discrete, unimodal, moderately right-skewed (p = 0.3)",
    mu = binom_size * binom_prob,
    sigma = sqrt(binom_size * binom_prob * (1 - binom_prob))
  ),
  `Cauchy` = list(
    pop_shape = "Continuous, symmetric, extremely heavy-tailed",
    mu = NA_real_,
    sigma = NA_real_
  )
)

# 2. Compute empirical metrics and compare with theory
summary_rows <- list()

for (dist_name in names(generators)) {
  gen <- generators[[dist_name]]
  theo <- theoretical_params[[dist_name]]
  
  for (n_val in n_values) {
    # Generate B sample means
    obs <- gen(B * n_val)
    mat <- matrix(obs, nrow = B, ncol = n_val)
    means <- rowMeans(mat)
    
    # Empirical sample statistics
    sim_mean <- mean(means)
    sim_se <- sd(means)
    stats <- shape_statistics(means)
    
    # Theoretical standard error: SE = sigma / sqrt(n)
    theo_se <- if (is.finite(theo$sigma)) theo$sigma / sqrt(n_val) else NA_real_
    
    # Normality test: Shapiro-Wilk on a subsample (max 5,000 observations allowed)
    sub_sample <- sample(means, min(5000L, B))
    shapiro_res <- shapiro.test(sub_sample)
    
    summary_rows[[length(summary_rows) + 1]] <- data.frame(
      distribution     = dist_name,
      n                = n_val,
      pop_shape        = theo$pop_shape,
      theo_mean        = round(theo$mu, 4),
      sim_mean         = round(sim_mean, 4),
      mean_bias        = round(sim_mean - theo$mu, 4),
      theo_se          = round(theo_se, 4),
      sim_se           = round(sim_se, 4),
      se_ratio         = round(sim_se / theo_se, 4),
      skewness         = round(unname(stats[["skewness"]]), 4),
      excess_kurtosis  = round(unname(stats[["excess_kurtosis"]]), 4),
      shapiro_stat_W   = round(shapiro_res$statistic, 4),
      shapiro_p_value  = round(shapiro_res$p.value, 5),
      stringsAsFactors = FALSE
    )
  }
}

clt_summary_table <- do.call(rbind, summary_rows)

# Save full comparison table to disk
write.csv(
  clt_summary_table,
  file = file.path(output_dir, "clt_convergence_metrics_table.csv"),
  row.names = FALSE
)

message("Summary table written to: ", file.path(output_dir, "clt_convergence_metrics_table.csv"))

# -----------------------------
# Graphs (ggplot2 Faceted 2x3 Grid: Skewness & Kurtosis)
# -----------------------------
plot_shape_facet <- function(data, statistic, y_label, main_title) {
  ggplot(data, aes(x = n, y = .data[[statistic]], color = distribution)) +
    # Reference target line at 0 for normal convergence
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50", linewidth = 0.5) +
    # Trajectory lines and points
    geom_line(linewidth = 0.8) +
    geom_point(size = 1.8) +
    # 2 rows x 3 columns with parametric y-scales per panel
    facet_wrap(~ distribution, nrow = 2, ncol = 3, scales = "free_y") +
    scale_color_manual(values = c(
      "Standard Normal" = "#0072B2",
      "Six-Sided Die"   = "#D55E00",
      "Exponential"     = "#009E73",
      "R-Binomial"      = "#CC79A7",
      "Cauchy"          = "#E69F00"
    )) +
    labs(
      title = main_title,
      x = "Sample size n",
      y = y_label
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5, margin = margin(b = 12)),
      strip.text = element_text(face = "bold", size = 11),
      strip.background = element_rect(fill = "gray94", color = "gray70"),
      panel.grid.minor = element_blank()
    )
}

# 1. Build shape plots
p_skewness <- plot_shape_facet(
  data = simulation_results,
  statistic = "skewness",
  y_label = "Skewness of sample means",
  main_title = "Monte Carlo Skewness of Sample Means by Distribution"
)

p_kurtosis <- plot_shape_facet(
  data = simulation_results,
  statistic = "excess_kurtosis",
  y_label = "Excess kurtosis of sample means",
  main_title = "Monte Carlo Excess Kurtosis of Sample Means by Distribution"
)

# 2. Display on screen (if enabled in settings)
if (isTRUE(displays_plots)) {
  print(p_skewness)
  print(p_kurtosis)
}

# 3. Export PNG files (2x3 grid per metric)
ggsave(
  filename = file.path(output_dir, "skewness_by_n.png"),
  plot = p_skewness,
  width = 12,
  height = 8,
  dpi = 200
)

ggsave(
  filename = file.path(output_dir, "excess_kurtosis_by_n.png"),
  plot = p_kurtosis,
  width = 12,
  height = 8,
  dpi = 200
)

# 4. Export combined PDF (Page 1 = Skewness, Page 2 = Kurtosis)
pdf(file.path(output_dir, "shape_statistics_by_n.pdf"), width = 12, height = 8)
print(p_skewness)
print(p_kurtosis)
dev.off()

# ---------------------------------------------------
# Additional Graph: Standardized Density vs N(0,1)
# ---------------------------------------------------
message("Generating Density Overlays...")

# Key sample sizes to display progression toward normality
key_n <- c(1, 5, 15, 50, 100)

density_data <- do.call(rbind, lapply(names(generators), function(dist_name) {
  gen <- generators[[dist_name]]
  
  do.call(rbind, lapply(key_n, function(n_val) {
    obs <- gen(B * n_val)
    mat <- matrix(obs, nrow = B, ncol = n_val)
    means <- rowMeans(mat)
    
    # Standardize means to Z ~ N(0, 1)
    mean_val <- mean(means)
    sd_val <- sd(means)
    
    # Handle Cauchy or edge cases with non-finite variance gracefully
    z_means <- if (is.finite(sd_val) && sd_val > 0) (means - mean_val) / sd_val else rep(NA_real_, B)
    
    data.frame(
      distribution = dist_name,
      n = factor(paste0("n = ", n_val), levels = paste0("n = ", key_n)),
      z = z_means
    )
  }))
})) %>% filter(is.finite(z))

p_density <- ggplot(density_data, aes(x = z)) +
  # Empirical sample densities for different n
  geom_density(aes(color = n, fill = n), alpha = 0.08, linewidth = 0.8) +
  # Ground truth standard Gaussian bell curve
  stat_function(
    fun = dnorm,
    args = list(mean = 0, sd = 1),
    linetype = "dashed",
    color = "black",
    linewidth = 0.85
  ) +
  # 2x3 Grid per distribution
  facet_wrap(~ distribution, nrow = 2, ncol = 3, scales = "free_y") +
  scale_color_brewer(palette = "Set1") +
  scale_fill_brewer(palette = "Set1") +
  coord_cartesian(xlim = c(-4, 4)) +
  labs(
    title = "Empirical Densities of Standardized Sample Means vs. Normal N(0, 1)",
    subtitle = "Colored areas: sample mean distributions across n | Dashed black curve: standard Normal target",
    x = "Standardized Sample Mean (Z-score)",
    y = "Density",
    color = "Sample Size (n)",
    fill = "Sample Size (n)"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "gray30", margin = margin(b = 12)),
    strip.text = element_text(face = "bold", size = 11),
    strip.background = element_rect(fill = "gray95", color = "gray75"),
    panel.grid.minor = element_blank(),
    legend.position = "bottom"
  )

if (isTRUE(displays_plots)) {
  print(p_density)
}

ggsave(
  filename = file.path(output_dir, "density_vs_normal_2x3_grid.png"),
  plot = p_density,
  width = 12,
  height = 8,
  dpi = 200
)
# ---------------------------------------------------
# Additional Graph: Normal Q-Q Plots for Every n Value
# ---------------------------------------------------
message("Generating Normal Q-Q plots for each n value...")

# Create a dedicated subfolder to keep output files organized
qq_dir <- file.path(output_dir, "qq_plots_by_n")
dir.create(qq_dir, showWarnings = FALSE, recursive = TRUE)

# 1. Simulate standardized sample means across ALL n_values
all_qq_data <- do.call(rbind, lapply(names(generators), function(dist_name) {
  gen <- generators[[dist_name]]
  
  do.call(rbind, lapply(n_values, function(n_val) {
    obs <- gen(B * n_val)
    mat <- matrix(obs, nrow = B, ncol = n_val)
    means <- rowMeans(mat)
    
    mean_val <- mean(means)
    sd_val <- sd(means)
    
    # Standardize means to Z ~ N(0, 1)
    z_means <- if (is.finite(sd_val) && sd_val > 0) (means - mean_val) / sd_val else rep(NA_real_, B)
    
    data.frame(
      distribution = dist_name,
      n = n_val,
      z = z_means
    )
  }))
})) %>% filter(is.finite(z))

# 2. Open a multi-page PDF to store every n on its own page
pdf_file <- file.path(output_dir, "all_qq_plots_by_n.pdf")
pdf(pdf_file, width = 12, height = 8)

# 3. Loop over each n, generate the 2x3 grid, print, and save individual PNGs
for (n_val in n_values) {
  sub_qq <- all_qq_data %>% filter(n == n_val)
  
  p_qq_n <- ggplot(sub_qq, aes(sample = z, color = distribution)) +
    # Target 45-degree reference line
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black", linewidth = 0.75) +
    stat_qq(size = 0.9, alpha = 0.45) +
    # 2x3 Grid containing all distributions for this specific n
    facet_wrap(~ distribution, nrow = 2, ncol = 3, scales = "free") +
    scale_color_manual(values = c(
      "Standard Normal" = "#0072B2",
      "Six-Sided Die"   = "#D55E00",
      "Exponential"     = "#009E73",
      "R-Binomial"      = "#CC79A7",
      "Cauchy"          = "#E69F00"
    )) +
    coord_cartesian(xlim = c(-3.5, 3.5), ylim = c(-3.5, 3.5)) +
    labs(
      title = paste0("Normal Q-Q Convergence: Sample Size n = ", n_val),
      subtitle = "Points aligning with the dashed diagonal indicate convergence to normality",
      x = "Theoretical Normal Quantiles",
      y = "Sample Quantiles (Standardized)"
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5, color = "gray30", margin = margin(b = 12)),
      strip.text = element_text(face = "bold", size = 11),
      strip.background = element_rect(fill = "gray95", color = "gray75"),
      panel.grid.minor = element_blank()
    )
  
  # Display in interactive window if enabled
  if (isTRUE(displays_plots)) {
    print(p_qq_n)
  }
  
  # Save to the multi-page PDF
  print(p_qq_n)
  
  # Save standalone individual PNG formatted with 3-digit padding (e.g. qq_plot_n_005.png)
  png_filename <- sprintf("qq_plot_n_%03d.png", n_val)
  ggsave(
    filename = file.path(qq_dir, png_filename),
    plot = p_qq_n,
    width = 12,
    height = 8,
    dpi = 180
  )
}

# Close the multi-page PDF device
dev.off()

message("Done. Files written to: ", normalizePath(output_dir))