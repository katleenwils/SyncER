# OxCal-style agreement index based on overlap integrals
# This is the actual method OxCal uses, not chi-squared

# Helper function to create PDF from samples
create_pdf <- function(ages, bandwidth = NULL, n = 512) {
  if (is.null(bandwidth)) {
    dens <- density(ages, n = n)
  } else {
    dens <- density(ages, bw = bandwidth, n = n)
  }

  # Normalize
  dx <- diff(dens$x)[1]
  dens$y <- dens$y / (sum(dens$y) * dx)

  return(list(x = dens$x, y = dens$y))
}

# OxCal-style Combine: Based on overlap integral F ratios
oxcal_combine_pdfs <- function(list1, list2, bandwidth = NULL, n_bins = 512) {
  # Create PDFs from the age samples
  pdf1 <- create_pdf(list1, bandwidth, n = n_bins)
  pdf2 <- create_pdf(list2, bandwidth, n = n_bins)

  # Create common grid
  x_min <- min(c(pdf1$x, pdf2$x))
  x_max <- max(c(pdf1$x, pdf2$x))
  x_common <- seq(x_min, x_max, length.out = n_bins)
  dx <- diff(x_common)[1]

  # Interpolate both distributions to common grid
  y1 <- approx(pdf1$x, pdf1$y, xout = x_common, rule = 2)$y
  y2 <- approx(pdf2$x, pdf2$y, xout = x_common, rule = 2)$y

  # Replace NAs with 0
  y1[is.na(y1)] <- 0
  y2[is.na(y2)] <- 0

  # Renormalize
  y1 <- y1 / (sum(y1) * dx)
  y2 <- y2 / (sum(y2) * dx)

  # Combined distribution (product of PDFs - this acts as "posterior")
  y_combined <- y1 * y2
  y_combined <- y_combined / (sum(y_combined) * dx)

  # Calculate F ratios (overlap integrals)
  # F_i = ∫ likelihood × posterior / ∫ likelihood × likelihood
  F1 <- sum(y1 * y_combined) * dx / (sum(y1 * y1) * dx)
  F2 <- sum(y2 * y_combined) * dx / (sum(y2 * y2) * dx)

  # Overall F (product of individual F ratios)
  F_overall <- F1 * F2

  # Agreement index (OxCal formula)
  # A_comb = 100 × F^(1/√n)
  n_distributions <- 2
  A_comb <- 100 * F_overall^(1/sqrt(n_distributions))

  # Individual agreement indices
  A1 <- 100 * F1
  A2 <- 100 * F2

  # Calculate thresholds
  threshold_chi_sq <- 100 / sqrt(2 * n_distributions)  # Chi-squared at 95%
  threshold_standard <- 60  # OxCal standard

  return(list(
    F1 = F1,
    F2 = F2,
    F_overall = F_overall,
    A1 = A1,
    A2 = A2,
    A_comb = A_comb,
    n_distributions = n_distributions,
    threshold_chi_sq = threshold_chi_sq,
    threshold_standard = threshold_standard,
    passes_chi_sq = A_comb > threshold_chi_sq,
    passes_standard = A_comb > threshold_standard,
    pdf1 = list(x = x_common, y = y1),
    pdf2 = list(x = x_common, y = y2),
    combined_pdf = list(x = x_common, y = y_combined)
  ))
}

# OxCal-style Difference: Subtract PDFs
oxcal_difference_pdfs <- function(list1, list2, bandwidth = NULL, n_bins = 1024) {
  # Create PDFs from samples
  pdf1 <- create_pdf(list1, bandwidth, n = n_bins)
  pdf2 <- create_pdf(list2, bandwidth, n = n_bins)

  # Compute difference distribution via convolution
  # P(difference = d) = ∫ P1(x) * P2(x+d) dx

  x_min_diff <- min(pdf2$x) - max(pdf1$x)
  x_max_diff <- max(pdf2$x) - min(pdf1$x)
  x_diff <- seq(x_min_diff, x_max_diff, length.out = n_bins)

  y_diff <- numeric(length(x_diff))

  for (i in seq_along(x_diff)) {
    d <- x_diff[i]
    x_shifted <- pdf2$x - d
    y2_interp <- approx(x_shifted, pdf2$y, xout = pdf1$x, rule = 2)$y
    y2_interp[is.na(y2_interp)] <- 0
    y_diff[i] <- sum(pdf1$y * y2_interp) * diff(pdf1$x)[1]
  }

  # Normalize
  dx_diff <- diff(x_diff)[1]
  y_diff <- y_diff / (sum(y_diff) * dx_diff)

  # Calculate statistics from the difference PDF
  mean_diff <- sum(x_diff * y_diff) * dx_diff
  variance_diff <- sum((x_diff - mean_diff)^2 * y_diff) * dx_diff
  sd_diff <- sqrt(variance_diff)

  # Mode
  mode_idx <- which.max(y_diff)
  mode_diff <- x_diff[mode_idx]

  # Credible intervals
  cumsum_prob <- cumsum(y_diff * dx_diff)

  idx_16 <- which.min(abs(cumsum_prob - 0.16))
  idx_84 <- which.min(abs(cumsum_prob - 0.84))
  ci_68 <- c(x_diff[idx_16], x_diff[idx_84])

  idx_2.5 <- which.min(abs(cumsum_prob - 0.025))
  idx_97.5 <- which.min(abs(cumsum_prob - 0.975))
  ci_95 <- c(x_diff[idx_2.5], x_diff[idx_97.5])

  # Significance test
  overlaps_zero <- ci_95[1] < 0 & ci_95[2] > 0

  # Probability that list2 is younger
  prob_positive <- sum(y_diff[x_diff < 0]) * dx_diff

  return(list(
    mean_difference = mean_diff,
    sd_difference = sd_diff,
    mode_difference = mode_diff,
    ci_68 = ci_68,
    ci_95 = ci_95,
    overlaps_zero = overlaps_zero,
    significant = !overlaps_zero,
    prob_list2_younger = prob_positive,
    diff_pdf = list(x = x_diff, y = y_diff)
  ))
}

# Multiple distributions
oxcal_combine_multiple_pdfs <- function(age_lists, bandwidth = NULL, n_bins = 512) {
  n <- length(age_lists)

  # Create PDFs
  pdfs <- lapply(age_lists, function(ages) create_pdf(ages, bandwidth, n = n_bins))

  # Common grid
  all_x <- unlist(lapply(pdfs, function(p) p$x))
  x_common <- seq(min(all_x), max(all_x), length.out = n_bins)
  dx <- diff(x_common)[1]

  # Interpolate all to common grid
  y_list <- lapply(pdfs, function(pdf) {
    y <- approx(pdf$x, pdf$y, xout = x_common, rule = 2)$y
    y[is.na(y)] <- 0
    y / (sum(y) * dx)
  })

  # Combined (posterior) distribution - product of all
  y_combined <- Reduce("*", y_list)
  y_combined <- y_combined / (sum(y_combined) * dx)

  # Calculate F ratios for each distribution
  F_values <- sapply(y_list, function(y) {
    numerator <- sum(y * y_combined) * dx
    denominator <- sum(y * y) * dx
    numerator / denominator
  })

  # Overall F
  F_overall <- prod(F_values)

  # Agreement indices
  A_individual <- 100 * F_values
  A_comb <- 100 * F_overall^(1/sqrt(n))

  # Calculate thresholds based on n
  threshold_chi_sq <- 100 / sqrt(2 * n)
  threshold_standard <- 60

  return(list(
    F_values = F_values,
    F_overall = F_overall,
    A_individual = A_individual,
    A_comb = A_comb,
    n_distributions = n,
    threshold_chi_sq = threshold_chi_sq,
    threshold_standard = threshold_standard,
    passes_chi_sq = A_comb > threshold_chi_sq,
    passes_standard = A_comb > threshold_standard,
    combined_pdf = list(x = x_common, y = y_combined),
    individual_pdfs = lapply(1:n, function(i) list(x = x_common, y = y_list[[i]]))
  ))
}

# Parnell et al. (2008) method - Simple age differences
parnell_age_differences <- function(age_list) {
  n <- length(age_list)
  core_names <- names(age_list)
  n_samples <- length(age_list[[1]])  # Assuming all have same number of samples

  pairwise_results <- list()

  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      pair_name <- paste(core_names[i], "vs", core_names[j])

      # Calculate age differences for each Monte Carlo sample
      age_diffs <- age_list[[j]] - age_list[[i]]

      # Calculate statistics
      mean_diff <- mean(age_diffs)
      sd_diff <- sd(age_diffs)
      median_diff <- median(age_diffs)

      # 95% HDR (Highest Posterior Density Region)
      ci_95 <- quantile(age_diffs, probs = c(0.025, 0.975))

      # Probability that j is younger than i (negative difference)
      prob_j_younger <- sum(age_diffs < 0) / length(age_diffs)

      # Test if zero is in 95% HDR (synchroneity test)
      overlaps_zero <- ci_95[1] < 0 & ci_95[2] > 0

      pairwise_results[[pair_name]] <- list(
        age_diffs = age_diffs,
        mean_diff = mean_diff,
        median_diff = median_diff,
        sd_diff = sd_diff,
        ci_95 = ci_95,
        prob_j_younger = prob_j_younger,
        overlaps_zero = overlaps_zero,
        significant = !overlaps_zero
      )
    }
  }

  return(pairwise_results)
}

# Calculate pairwise differences for all cores (OxCal Difference method)
calculate_pairwise_differences <- function(age_list) {
  n <- length(age_list)
  core_names <- names(age_list)

  diff_results <- list()

  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      pair_name <- paste(core_names[i], "vs", core_names[j])
      diff_results[[pair_name]] <- oxcal_difference_pdfs(age_list[[i]], age_list[[j]])
    }
  }

  return(diff_results)
}