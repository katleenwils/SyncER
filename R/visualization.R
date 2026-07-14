#' Create Comparison Histogram Visualizations
#'
#' Creates publication-quality comparison visualizations showing age probability distributions
#' and log-ratio analysis for two records.
#'
#' @param Amin Numeric minimum age range for record A.
#' @param Amax Numeric maximum age range for record A.
#' @param Bmin Numeric minimum age range for record B.
#' @param Bmax Numeric maximum age range for record B.
#' @param ABlr Numeric vector of log-ratio values (log(A/B)).
#' @param sheet1_name Character string identifying record A.
#' @param sheet2_name Character string identifying record B.
#' @param column_name Character string identifying the event/horizon being compared.
#' @param col1 Numeric vector of posterior age samples for record 1 (default: NULL).
#' @param col2 Numeric vector of posterior age samples for record 2 (default: NULL).
#' @param age_diff_log_bounds Two-element numeric vector specifying log-space age difference
#'   thresholds (shown as green lines).
#' @param conf_log_bounds Two-element numeric vector specifying confidence interval bounds
#'   in log-space (shown as blue shading).
#' @param conf_level Numeric value specifying the confidence level used (e.g., 0.95).
#' @param variant1 Optional character string specifying exact horizon name in record 1
#'   (for generic horizon groups) (default: NULL).
#' @param variant2 Optional character string specifying exact horizon name in record 2
#'   (for generic horizon groups) (default: NULL).
#' @param plot_opts Named list of plot layout options:
#'   \itemize{
#'     \item \code{n_bins}: histogram resolution (default: \code{500})
#'     \item \code{base_size}: ggplot base font size in pt (default: \code{10})
#'     \item \code{label_size}: record label text size (default: \code{3})
#'     \item \code{annotation_size}: statistical annotation text size (default: \code{3})
#'     \item \code{minor_break_by}: x-axis minor tick spacing in cal yrs (default: \code{100})
#'     \item \code{x_n_breaks}: approximate number of major x-axis breaks (default: \code{10})
#'     \item \code{age_axis_label}: age distribution panel's x-axis label
#'           (default: \code{"Age (cal yrs BP)"})
#'   }
#'
#' @return No return value. Creates a composite plot with two panels displayed in the current
#'   graphics device.
#'
#' @details Generates a two-panel visualization: (1) histogram of log-ratios with fitted normal
#'   distribution, confidence intervals, and age difference thresholds; (2) age probability
#'   distributions for both records shown as histograms with frequency gradients.
#'
#' @export
create_visualization <- function(Amin,
                                 Amax,
                                 Bmin,
                                 Bmax,
                                 ABlr,
                                 sheet1_name,
                                 sheet2_name,
                                 column_name,
                                 col1 = NULL,
                                 col2 = NULL,
                                 age_diff_log_bounds = NULL,
                                 conf_log_bounds  = NULL,
                                 conf_level=0.95,
                                 variant1 = NULL,
                                 variant2 = NULL,
                                 plot_opts = list()) {

  opts <- modifyList(
    list(
      n_bins          = 500,                # histogram resolution
      base_size       = 10,                 # ggplot base font size (pt)
      label_size      = 3,                  # record label text size
      annotation_size = 3,                  # statistical annotation text size
      minor_break_by  = 100,                # x-axis minor tick spacing (cal yrs)
      x_n_breaks      = 10,                 # approximate number of major x-axis breaks
      age_axis_label  = "Age (cal yrs BP)"  # age distribution panel's x-axis label
    ),
    plot_opts
  )

  n_bins <- opts$n_bins

  # Create bin edges for both records
  Abins <- seq(Amin, Amax, length.out = n_bins + 1)
  Bbins <- seq(Bmin, Bmax, length.out = n_bins + 1)

  # Calculate histograms (counts per bin)
  A_counts <- hist(col1, breaks = Abins, plot = FALSE)$counts
  B_counts <- hist(col2, breaks = Bbins, plot = FALSE)$counts

  # Find maximum count for color scaling
  max_count <- max(c(A_counts, B_counts), na.rm = TRUE)

  # Pre-compute bin lengths
  n_Abins <- length(Abins)
  n_Bbins <- length(Bbins)

  # Create data frames for histogram rectangles
  A_df <- data.frame(
    xmin = Abins[-n_Abins],
    xmax = Abins[-1],
    ymin = 1.75,  # Y-position for record A
    ymax = 2.25,
    fill = A_counts
  )

  B_df <- data.frame(
    xmin = Bbins[-n_Bbins],
    xmax = Bbins[-1],
    ymin = 0.75,  # Y-position for record B
    ymax = 1.25,
    fill = B_counts
  )

  # Helper function to safely check variant labels
  # Returns empty string if variant is NULL or NA, otherwise returns formatted label
  format_variant <- function(variant, column_name) {
    if (is.null(variant) || is.na(variant) || variant == column_name) {
      return("")
    } else {
      return(paste0(" (", variant, ")"))
    }
  }

  # PLOT 1: Age distributions for both records
  p1 <- ggplot2::ggplot() +
    # Draw histogram rectangles with gradient fill
    ggplot2::geom_rect(data = A_df, ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill),
              color = NA) +
    ggplot2::geom_rect(data = B_df, ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill),
              color = NA) +
    # Draw borders around distributions
    ggplot2::geom_rect(ggplot2::aes(xmin = Amin, xmax = Amax, ymin = 1.75, ymax = 2.25),
              fill = NA, color="black") +
    ggplot2::geom_rect(ggplot2::aes(xmin = Bmin, xmax = Bmax, ymin = 0.75, ymax = 1.25),
              fill = NA, color="black") +
    # Color scale (white to black based on frequency)
    ggplot2::scale_fill_gradient(low = "white", high = "black", limits = c(0, max_count),
                        na.value = "white") +
    # Add record labels with variant info (safely formatted)
    ggplot2::annotate("text", x = Inf, y = 2.35,
             label = paste0("Age PDF: ", sheet1_name, format_variant(variant1, column_name)),
             hjust = -0.1, vjust = 0, size = 3) +
    ggplot2::annotate("text", x = Inf, y = 1.35,
             label = paste0("Age PDF: ", sheet2_name, format_variant(variant2, column_name)),
             hjust = -0.1, vjust = 0, size = 3) +
    ggplot2::xlab(opts$age_axis_label) +
    ggplot2::ylim(0, 3) +
    ggplot2::theme_minimal() +
    ggplot2::theme(
      panel.grid = ggplot2::element_blank(),
      axis.title.y = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_blank(),
      axis.ticks.y = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(color = "black"),
      panel.border = ggplot2::element_rect(color = "black", fill = NA),
      legend.position = "none"
    ) +
    ggplot2::scale_x_continuous(
      breaks = pretty(c(Amin, Amax, Bmin, Bmax), n = opts$x_n_breaks),
      minor_breaks = seq(min(c(Amin, Bmin)), max(c(Amax, Bmax)), by = opts$minor_break_by),
      trans = "reverse"
    )

  # Fit normal distribution to log-ratios
  fit_mean <- mean(ABlr, na.rm = TRUE)
  fit_sd <- sd(ABlr, na.rm = TRUE)
  fit_median <- median(ABlr, na.rm = TRUE)

  # Create fitted normal density curve
  x_vals <- seq(min(ABlr, na.rm = TRUE), max(ABlr, na.rm = TRUE), length.out = 1000)
  y_vals <- dnorm(x_vals, mean = fit_mean, sd = fit_sd)

  # Scale density to match histogram counts
  range_ABlr <- max(ABlr, na.rm = TRUE) - min(ABlr, na.rm = TRUE)
  binwidth <- range_ABlr / 50
  scale_factor <- length(ABlr) * binwidth
  y_vals_scaled <- y_vals * scale_factor
  fit_df <- data.frame(x = x_vals, y = y_vals_scaled)

  # PLOT 2: Log-ratio histogram with statistical overlays
  df <- data.frame(ABlr = ABlr)
  p2 <- ggplot2::ggplot(df, ggplot2::aes(x = ABlr)) +
    # Histogram with conditional coloring (blue = within CI, light blue = outside)
    ggplot2::geom_histogram(
      binwidth = binwidth,
      ggplot2::aes(fill = (ggplot2::after_stat(xmin) >= conf_log_bounds[1] &
                    ggplot2::after_stat(xmax) <= conf_log_bounds[2])),
      color = "black"
    ) +
    ggplot2::scale_fill_manual(values = c("TRUE" = "blue", "FALSE" = "lightblue"), guide = "none") +
    # Overlay fitted normal distribution
    ggplot2::geom_line(data = fit_df, ggplot2::aes(x = x, y = y), color = "blue", linewidth = 1) +
    # Add mean line
    ggplot2::geom_vline(xintercept = fit_mean, linetype = "dashed", color = "black", linewidth = 0.8) +
    ggplot2::geom_vline(xintercept = fit_median, linetype = "dotted", color = "darkred", linewidth = 0.8) +
    # Add age difference threshold lines (red)
    ggplot2::geom_vline(xintercept = age_diff_log_bounds[1], linetype = "dashed",
               color = "green", linewidth = 0.9) +
    ggplot2::geom_vline(xintercept = age_diff_log_bounds[2], linetype = "dashed",
               color = "green", linewidth = 0.9) +
    # Annotate statistical markers
    ggplot2::annotate("text", x = fit_mean, y = Inf, label = "Mean",
             vjust = -0.5, hjust = 0.5, color = "black", size = 3) +
    ggplot2::annotate("text", x = fit_median, y = Inf, label = "Median",
             vjust = -1.8, hjust = 0.5, color = "darkred", size = 3) +
    ggplot2::annotate("text", x = conf_log_bounds[2], y = Inf,
             label = paste0("+", round(conf_level*100), "% CI"),
             vjust = -0.5, hjust = 0, color = "blue", size = 3) +
    ggplot2::annotate("text", x = conf_log_bounds[1], y = Inf,
             label = paste0("-", round(conf_level*100), "% CI"),
             vjust = -0.5, hjust = 1, color = "blue", size = 3) +
    # Add statistics text box
    ggplot2::annotate(
      "text",
      x = max(ABlr, na.rm = TRUE),
      y = max(y_vals_scaled) * 1.05,
      label = paste0("Mean = ", round(fit_mean, 2),
                     "\nMedian = ", round(fit_median, 2),
                     "\nSD = ", round(fit_sd, 2)),
      hjust = 1, vjust = 1,
      size = opts$annotation_size + 0.5, color = "black"
    ) +
    ggplot2::theme_minimal(base_size = opts$base_size) +
    ggplot2::theme(
      aspect.ratio = 1,
      panel.grid = ggplot2::element_blank(),
      plot.margin = ggplot2::margin(5, 5, 5, 5),
      axis.line = ggplot2::element_line(color = "black"),
      panel.border = ggplot2::element_rect(color = "black", fill = NA),
      plot.title = ggplot2::element_text(hjust = 0.5),
      axis.title.y = ggplot2::element_blank()
    ) +
    ggplot2::xlab(paste0("Log(", sheet1_name, "/", sheet2_name, ")")) +
    ggplot2::ylab("Frequency")

  invisible(list(log_ratio = p2, age_dist = p1))
}

#' Plot Synchronicity Comparison Figures
#'
#' Produces per-horizon PDF files and console visualizations from a
#' \code{compute_synchronicity()} result.
#'
#' @param synchro_result The list returned by \code{compute_synchronicity()}.
#' @param output_dir Character string specifying the directory the PDF files will be
#'   saved to (default: \code{syncer_output_dir()}, i.e. the \code{SyncER_outputs}
#'   folder in the working directory).
#' @param offset Numeric offset correction value applied to ages before plotting
#'   (default: \code{bp_datum()}, i.e. the current year minus 1950).
#' @param synced Character string appended to PDF file names (default: "").
#' @param fig_width Numeric width of output PDF figures in inches (default: 10).
#' @param fig_height Numeric height of output PDF figures in inches (default: 8).
#' @inheritParams create_visualization
#'
#' @return Invisibly returns \code{NULL}.
#'
#' @export
plot_synchronicity <- function(synchro_result,
                               output_dir = syncer_output_dir(),
                               offset     = bp_datum(),
                               synced     = "",
                               fig_width  = 10,
                               fig_height = 8,
                               plot_opts  = list()) {

  viz_data <- synchro_result$viz_data
  if (is.null(viz_data) || nrow(viz_data) == 0) return(invisible(NULL))

  cat("Printing comparison PDFs...\n")
  out_dir <- output_dir
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  horizons <- unique(viz_data$Event)
  for (horizon in horizons) {
    horizon_rows <- viz_data[viz_data$Event == horizon, ]

    horizon_plots <- vector("list", nrow(horizon_rows))
    for (row_idx in seq_len(nrow(horizon_rows))) {
      row   <- horizon_rows[row_idx, ]
      plots <- create_visualization(
        Amin = row$Amin - offset,
        Amax = row$Amax - offset,
        Bmin = row$Bmin - offset,
        Bmax = row$Bmax - offset,
        ABlr = unlist(row$ABlr),
        sheet1_name = row$record1,
        sheet2_name = row$record2,
        column_name = row$Event,
        col1 = unlist(row$col1) - offset,
        col2 = unlist(row$col2) - offset,
        age_diff_log_bounds = c(row$age_diff_log_bounds_low, row$age_diff_log_bounds_high),
        conf_log_bounds     = c(row$conf_log_bounds_low,     row$conf_log_bounds_high),
        conf_level  = row$conf_level,
        variant1    = row$variant1,
        variant2    = row$variant2,
        plot_opts   = plot_opts
      )
      horizon_plots[[row_idx]] <- plots
      # Print to console outside any PDF device
      print(plots$log_ratio)
      print(plots$age_dist)
    }

    pdf_path <- file.path(out_dir, paste0(horizon, synced, ".pdf"))
    pdf(file = pdf_path, width = fig_width, height = fig_height)
    for (plots in horizon_plots) {
      print(plots$log_ratio)
      print(plots$age_dist)
    }
    dev.off()
  }

  invisible(NULL)
}

#' Plot Bayesian Age Combination Results
#'
#' Draws the Bayesian posterior combination figure for one horizon group from
#' the raw data returned by \code{compute_synchronized_ages()$bayesian_plot_data}.
#' Prints to the active graphics device (console) and, when \code{output_dir} is
#' supplied, also saves a PDF.
#'
#' @param pd Named list for a single horizon, as stored in
#'   \code{compute_synchronized_ages()$bayesian_plot_data}. Must contain elements
#'   \code{samples_list_shifted}, \code{combined_pdf_x}, \code{combined_pdf_vals},
#'   \code{valid_records}, \code{mu_comb}, \code{sigma_comb}, \code{group_name},
#'   and \code{bayes_opts}.
#' @param output_dir Character string; directory the PDF is saved to. Pass
#'   \code{NULL} to skip PDF output.
#'
#' @return Invisibly returns \code{NULL}.
#'
#' @export
plot_bayesian_combination <- function(pd, output_dir = NULL) {

  samples_list_shifted <- pd$samples_list_shifted
  combined_pdf_x       <- pd$combined_pdf_x
  combined_pdf_vals    <- pd$combined_pdf_vals
  valid_records        <- pd$valid_records
  mu_comb              <- pd$mu_comb
  sigma_comb           <- pd$sigma_comb
  group_name           <- pd$group_name
  bayes_opts           <- pd$bayes_opts

  all_samples <- unlist(samples_list_shifted)
  xrange <- range(all_samples,
                  mu_comb + c(-bayes_opts$plot_range_sigma * sigma_comb,
                               bayes_opts$plot_range_sigma * sigma_comb))

  max_density <- max(
    vapply(samples_list_shifted, function(s) max(density(s)$y), numeric(1)),
    max(combined_pdf_vals)
  )

  draw_plot <- function() {
    plot(0, 0, type = "n",
         xlim = xrange, ylim = c(0, max_density),
         xlab = "Age (cal yrs BP)", ylab = "Density",
         main = paste("Bayesian Combination for", group_name))

    cols <- rainbow(length(samples_list_shifted))
    for (i in seq_along(samples_list_shifted)) {
      d <- density(samples_list_shifted[[i]])
      lines(d$x, d$y, col = cols[i], lwd = bayes_opts$posterior_lwd)
    }
    lines(combined_pdf_x, combined_pdf_vals, col = "black", lwd = bayes_opts$combined_lwd)
    legend(bayes_opts$legend_pos,
           legend = c(valid_records, "Bayesian combined"),
           col    = c(cols, "black"),
           lwd    = c(rep(bayes_opts$posterior_lwd, length(cols)), bayes_opts$combined_lwd))
  }

  draw_plot()

  if (!is.null(output_dir)) {
    dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
    pdf_path <- file.path(output_dir,
                          paste0(group_name, "_Bayesian_combination.pdf"))
    pdf(pdf_path, width = bayes_opts$fig_width, height = bayes_opts$fig_height)
    draw_plot()
    dev.off()
    message("Saved Bayesian plot: ", pdf_path)
  }

  invisible(NULL)
}
