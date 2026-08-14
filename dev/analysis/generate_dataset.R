# Set the seed for reproducibility
set.seed(125)

required_packages <- c("ggplot2", "dplyr", "gridExtra", "cowplot", "rintcal")
missing_packages  <- required_packages[!required_packages %in% installed.packages()[, "Package"]]
if (length(missing_packages) > 0) install.packages(missing_packages)

library(ggplot2)
library(dplyr)
library(gridExtra)
library(cowplot)

# All output (core_data/ CSVs, core_plots.pdf) is written to out_dir. By default
# this is the current working directory; change it to write elsewhere.
out_dir <- getwd()
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Load the IntCal20 calibration curve (Reimer et al. 2020) via rintcal.
intcal_data <- rintcal::ccurve()
colnames(intcal_data) <- c("cal_BP", "C14_age", "sigma")

# Generate reproducible sedimentation rates for each core
set.seed(126)
base_sedrates <- runif(5, min = 0.01, max = 0.1)
names(base_sedrates) <- paste0("core", 1:5)

# -------------------------------
# Define per-core sedimentation rate segments
# Each segment: list(depth_start, depth_end, sedrate_multiplier)
# depth_start and depth_end are proportions of total core length (0 to 1)
# sedrate_multiplier scales the base sedrate for that segment
# -------------------------------
sedrate_segments_list <- list(
  core1 = data.frame(
    prop_start = c(0.0, 0.3),
    prop_end   = c(0.3, 1.0),
    multiplier = c(1.0, 1/1.8)  # below 30%: slower sedimentation
  ),
  core2 = data.frame(
    prop_start = c(0.0),
    prop_end   = c(1.0),
    multiplier = c(1.0)
  ),
  core3 = data.frame(
    prop_start = c(0.0, 0.2, 0.6),
    prop_end   = c(0.2, 0.6, 1.0),
    multiplier = c(1.0, 1/0.5, 1/1.3)
  ),
  core4 = data.frame(
    prop_start = c(0.0),
    prop_end   = c(1.0),
    multiplier = c(1.0)
  ),
  core5 = data.frame(
    prop_start = c(0.0, 0.7),
    prop_end   = c(0.7, 1.0),
    multiplier = c(1.0, 1/2.0)
  )
)

# -------------------------------
# Function to convert age to depth given segments
# Segments are defined in depth-proportion space, but we build
# the age-depth relationship by computing cumulative age per segment
# -------------------------------
build_age_depth_model <- function(base_sedrate, segments, max_depth) {
  segments$depth_start <- segments$prop_start * max_depth
  segments$depth_end   <- segments$prop_end   * max_depth
  segments$sedrate     <- base_sedrate * segments$multiplier

  segments$age_start <- 0

  if (nrow(segments) > 1) {
    for (i in 2:nrow(segments)) {
      seg_depth <- segments$depth_end[i-1] - segments$depth_start[i-1]
      segments$age_start[i] <- segments$age_start[i-1] + seg_depth / segments$sedrate[i-1]
    }
  }

  return(segments)
}

depth_to_age <- function(depth, adm) {
  # For each depth value, find which segment it falls in and compute age
  sapply(depth, function(d) {
    seg <- adm[d >= adm$depth_start & d <= adm$depth_end, ]
    if (nrow(seg) == 0) seg <- adm[nrow(adm), ]  # clamp to last segment
    seg <- seg[1, ]
    seg$age_start + (d - seg$depth_start) / seg$sedrate
  })
}

age_to_depth <- function(age, adm) {
  # Compute the age at end of each segment
  adm$age_end <- sapply(1:nrow(adm), function(i) {
    adm$age_start[i] + (adm$depth_end[i] - adm$depth_start[i]) / adm$sedrate[i]
  })

  sapply(age, function(a) {
    seg <- adm[a >= adm$age_start & a <= adm$age_end, ]
    if (nrow(seg) == 0) seg <- adm[nrow(adm), ]
    seg <- seg[1, ]
    seg$depth_start + (a - seg$age_start) * seg$sedrate
  })
}

# -------------------------------
# Generate non-synchronous depths (same as before)
# -------------------------------
generate_positive_non_synchros <- function(n, mean_age, sd_age, sedrate, seed) {
  set.seed(seed)
  age_increments <- abs(rnorm(n, mean = mean_age, sd = sd_age))
  depth_cum <- cumsum(age_increments * sedrate)
  return(round(depth_cum))
}

cores <- lapply(1:5, function(i) {
  set.seed(200 + i)
  num_non_synchros <- sample(5:20, 1)
  sort(generate_positive_non_synchros(num_non_synchros,
                                      mean_age = 500,
                                      sd_age = 200,
                                      sedrate = base_sedrates[i],
                                      seed = 300 + i))
})
names(cores) <- paste0("core", 1:5)

core_depths <- cores
core_data_list <- list()

synchro_recurrence <- 300
chrono_recurrence  <- 1000

# -------------------------------
# Generate all events for a core using the age-depth model
# -------------------------------
generate_core_data <- function(core_index, core_name, non_synchro_depths,
                               base_sedrate, segments,
                               chrono_recurrence, synchro_recurrence) {
  if (length(non_synchro_depths) == 0) return(NULL)

  set.seed(400 + core_index)

  max_depth <- non_synchro_depths[length(non_synchro_depths)]

  # Build age-depth model for this core
  adm <- build_age_depth_model(base_sedrate, segments, max_depth)

  # Non-synchronous events: depth is fixed, age derived from ADM
  non_synchro_ages <- round(depth_to_age(non_synchro_depths, adm))
  non_synchros <- data.frame(
    depth = non_synchro_depths,
    event = "non-synchronous",
    age   = non_synchro_ages
  )

  # Max age of the core
  max_age <- round(depth_to_age(max_depth, adm))

  # isochrons: defined by age, depth derived from ADM
  max_chrono_index <- ceiling((max_age - 175) / chrono_recurrence)
  chrono_ages  <- 175 + (0:(max_chrono_index - 1)) * chrono_recurrence
  chrono_depths <- round(age_to_depth(chrono_ages, adm))
  chronos <- data.frame(
    depth = chrono_depths,
    event = paste0("isochron", seq_along(chrono_depths)),
    age   = chrono_ages
  )

  # Synchronous events: defined by age, depth derived from ADM
  max_synchro_index <- ceiling((max_age - 100) / synchro_recurrence) + 1
  synchro_ages  <- 100 + (0:(max_synchro_index - 1)) * synchro_recurrence
  synchro_depths <- round(age_to_depth(synchro_ages, adm))
  synchros <- data.frame(
    depth = synchro_depths,
    event = "synchronous",
    age   = synchro_ages
  )

  # Samples: evenly spaced in age, depth derived from ADM
  sediment_years <- ceiling(max_age / 500) * 500
  set.seed(400 + core_index)
  sample_ages  <- round(sort(runif(sediment_years / 500, 0, sediment_years)))
  sample_depths <- round(age_to_depth(sample_ages, adm))
  sample_events <- data.frame(
    depth = sample_depths,
    event = "sample",
    age   = sample_ages
  )

  # Combine and sort
  all_events <- rbind(non_synchros, chronos, synchros, sample_events)
  all_events <- all_events[order(all_events$depth), ]
  n <- nrow(all_events)

  all_events$relative_depth <- seq(1, n)
  all_events$calendar_age   <- 2025 - all_events$age
  all_events$BP             <- 1950 - all_events$calendar_age

  all_events$C14_age <- sapply(all_events$BP, function(bp_val) {
    matched_row <- intcal_data[which.min(abs(intcal_data$cal_BP - bp_val)), ]
    return(matched_row$C14_age)
  })

  set.seed(600 + core_index)
  all_events$C14_error <- round(runif(n, min = 20, max = 50))

  all_events <- all_events[, c("relative_depth", "depth", "age", "calendar_age",
                               "BP", "C14_age", "C14_error", "event")]

  return(list(data = all_events, max_depth = max_depth, adm = adm))
}

# -------------------------------
# Run and save
# -------------------------------
core_data_dir <- file.path(out_dir, "core_data")
dir.create(core_data_dir, showWarnings = FALSE, recursive = TRUE)
max_depths <- list()

for (i in seq_along(core_depths)) {
  core_name   <- names(core_depths)[i]
  core_output <- generate_core_data(
    core_index        = i,
    core_name         = core_name,
    non_synchro_depths = core_depths[[i]],
    base_sedrate      = base_sedrates[core_name],
    segments          = sedrate_segments_list[[core_name]],
    chrono_recurrence = chrono_recurrence,
    synchro_recurrence = synchro_recurrence
  )

  df <- core_output$data
  max_depths[[core_name]] <- max(df$depth)
  core_data_list[[core_name]] <- df

  if (!is.null(df)) {
    write.csv(df, file.path(core_data_dir, paste0(core_name, ".csv")), row.names = FALSE)
  }
}

cat("Data successfully saved to", core_data_dir, "\n")

# -------------------------------
# Plotting (unchanged)
# -------------------------------
csv_files <- list.files(core_data_dir, pattern = "\\.csv$", full.names = TRUE)
sheets <- tools::file_path_sans_ext(basename(csv_files))
core_data_list <- lapply(csv_files, read.csv)
names(core_data_list) <- sheets

event_colors <- c(
  "non-synchronous" = "#769897",
  "isochron"   = "#E51A4B",
  "synchronous"     = "#FFF7B2",
  "sample"          = "#393185",
  "sediment"        = "burlywood"
)

prepare_core_data <- function(core_data, min_height = 0.5) {
  core_data <- core_data %>% arrange(depth)
  core_data$legend_event <- sapply(core_data$event, function(ev) {
    if (grepl("isochron", ev, ignore.case = TRUE)) return("isochron")
    return(ev)
  })
  core_data$legend_event <- factor(core_data$legend_event,
                                   levels = c("non-synchronous","isochron","synchronous","sample","sediment"))
  event_layers <- data.frame(
    top          = core_data$depth + min_height/2,
    bottom       = core_data$depth - min_height/2,
    legend_event = core_data$legend_event
  )
  sediment_intervals <- data.frame(
    top          = 0,
    bottom       = max(core_data$depth),
    legend_event = factor("sediment", levels = levels(core_data$legend_event))
  )
  full_core <- bind_rows(sediment_intervals, event_layers) %>% arrange(top)
  return(full_core)
}

plot_core <- function(core_name, core_data) {
  prepared_data <- prepare_core_data(core_data)
  ggplot(prepared_data, aes(x = 1, ymin = bottom, ymax = top, fill = legend_event, color = legend_event)) +
    geom_rect(aes(xmin = 0.8, xmax = 1.2), linewidth = 2) +
    scale_fill_manual(values = event_colors) +
    scale_color_manual(values = event_colors) +
    scale_y_reverse() +
    theme_minimal() +
    theme(axis.title.x = element_blank(),
          axis.text.x  = element_blank(),
          axis.ticks.x = element_blank()) +
    labs(y = "Depth (cm)", title = core_name)
}

pdf(file.path(out_dir, "core_plots.pdf"), width = 8, height = 20)
core_plots <- lapply(names(core_data_list), function(core) {
  plot_core(core, core_data_list[[core]]) +
    theme(legend.position  = "bottom",
          legend.direction = "horizontal",
          legend.box       = "horizontal",
          legend.title     = element_blank(),
          legend.text      = element_text(size = 12),
          legend.spacing.x = unit(0.5, 'cm')) +
    guides(color = guide_legend(nrow = 1, byrow = TRUE))
})

get_legend <- function(my_plot) {
  tmp <- ggplot_gtable(ggplot_build(my_plot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  return(tmp$grobs[[leg]])
}

legend       <- get_legend(core_plots[[1]])
core_plots   <- lapply(core_plots, function(p) p + theme(legend.position = "none"))
plot_grid_combined <- plot_grid(plotlist = core_plots, ncol = 3)
final_plot   <- plot_grid(plot_grid_combined, legend, ncol = 1, rel_heights = c(1, 0.1))
print(final_plot)
dev.off()
cat("Plots successfully saved to", file.path(out_dir, "core_plots.pdf"), "\n")
