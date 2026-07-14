#' Add Adjusted Horizons to Bacon CSV Files
#'
#' Adds synchronized horizon ages to existing Bacon input CSV files, creating updated
#' versions for re-running age-depth models.
#'
#' @param parent_folder Character string specifying the root directory containing Bacon
#'   input folders.
#' @param record_data Named list of record metadata (output from \code{load_excel_data()}).
#' @param adjusted_ages Output from \code{set_to_zero()} containing synchronized ages.
#' @param non_synchro Logical indicating whether to print detailed information about added
#'   samples (default: FALSE).
#' @param overwrite Logical indicating whether to overwrite existing files (TRUE) or create
#'   new folders with _synced suffix (FALSE, default).
#' @param synced Logical indicating whether input folders already have _synced suffix (TRUE)
#'   or are original folders (FALSE, default).
#' @param update_records Character string suffix identifying folders to update when
#'   \code{synced = TRUE} (default: "_synced").
#' @param age_record Character string suffix to strip from record names in \code{adjusted_ages}
#'   for matching (default: "_synced").
#'
#' @return No return value. Modifies or creates CSV files and prints progress messages.
#'
#' @details When \code{synced = FALSE} (default):
#'   - Searches for folders matching record names from \code{record_data}
#'   - If \code{overwrite = TRUE}: modifies CSV files in original folders
#'   - If \code{overwrite = FALSE}: creates new folders with "_synced" suffix
#'
#'   When \code{synced = TRUE}:
#'   - Searches for folders with \code{update_records} suffix
#'   - If \code{overwrite = TRUE}: modifies CSV files in _synced folders
#'   - If \code{overwrite = FALSE}: creates new folders with incremental numeric suffixes
#'     (e.g., "core1_synced_1", "core1_synced_2")
#'
#' @export
add_adjusted_horizons_to_csv <- function(
    parent_folder,
    record_data,
    adjusted_ages,
    non_synchro = FALSE,
    overwrite = FALSE,
    synced = FALSE,
    update_records = "_synced",
    age_record = "_synced") {

  # Snapshot next_n for ALL folders BEFORE any writing begins
  if (synced && !overwrite) {
    all_top_dirs <- list.dirs(parent_folder, recursive = FALSE, full.names = TRUE)
    pattern      <- paste0(update_records, "$")
    target_folders <- all_top_dirs[grepl(pattern, basename(all_top_dirs))]

    if (length(target_folders) == 0) {
      cat("WARNING: No folders found matching pattern '",
          paste0(update_records, "$"), "'\n", sep = "")
      return(invisible(NULL))
    }

    # Pre-calculate next_n for each target folder NOW before any folders are created
    folder_next_n <- setNames(vapply(target_folders, function(folder) {
      base_name    <- sub(paste0(update_records, "$"), "", basename(folder))
      all_siblings <- basename(list.dirs(dirname(folder), recursive = FALSE))
      suffixes <- suppressWarnings(
        as.integer(sub(paste0("^", base_name, update_records, "_"), "",
                       all_siblings[grepl(paste0("^", base_name, update_records, "_[0-9]+$"),
                                          all_siblings)]))
      )
      if (length(suffixes) == 0) 1L else max(suffixes, na.rm = TRUE) + 1L
    }, integer(1)), target_folders)

  } else if (!synced) {
    all_top_dirs   <- list.dirs(parent_folder, recursive = FALSE, full.names = TRUE)
    target_folders <- character(0)
    for (record_name in names(record_data)) {
      matched <- all_top_dirs[basename(all_top_dirs) == record_name]
      if (length(matched) > 0) target_folders <- c(target_folders, matched[1])
    }
    if (length(target_folders) == 0) {
      cat("WARNING: No folders found matching record names from record_data\n")
      return(invisible(NULL))
    }
    folder_next_n <- NULL
  } else {
    # synced + overwrite
    all_top_dirs   <- list.dirs(parent_folder, recursive = FALSE, full.names = TRUE)
    pattern        <- paste0(update_records, "$")
    target_folders <- all_top_dirs[grepl(pattern, basename(all_top_dirs))]
    folder_next_n  <- NULL
  }

  # Process each folder
  for (folder in target_folders) {

    csv_files <- list.files(folder, pattern = "\\.csv$", full.names = TRUE)

    if (length(csv_files) == 0) {
      cat("WARNING: No CSV files found in", folder, "\n")
      next
    }

    # Determine new folder path using pre-calculated next_n (but don't create yet)
    if (!overwrite) {
      if (synced) {
        base_name  <- sub(paste0(update_records, "$"), "", basename(folder))
        next_n     <- folder_next_n[[folder]]
        new_folder <- file.path(dirname(folder),
                                paste0(base_name, update_records, "_", next_n))
      } else {
        first_csv  <- tools::file_path_sans_ext(basename(csv_files[1]))
        core_name  <- sub("_C14$", "", first_csv)
        new_folder <- file.path(dirname(folder), paste0(core_name, update_records))
      }
    }

    # Determine whether a _C14 file exists in this folder
    has_c14_file <- any(grepl("_C14\\.csv$", csv_files))

    # Collect outputs in memory first
    outputs      <- list()
    any_modified <- FALSE

    for (file in csv_files) {

      core_name_full <- tools::file_path_sans_ext(basename(file))

      # Determine core names
      if (synced) {
        core_name_synced <- sub("_C14$", "", core_name_full)
        core_name        <- sub(paste0(update_records, "$"), "", core_name_synced)
      } else {
        core_name        <- basename(folder)
        core_name_synced <- core_name
      }

      if (!core_name %in% names(record_data)) {
        cat("WARNING: Core name '", core_name, "' not found in record_data. Skipping.\n", sep = "")
        next
      }

      # Decide whether this file gets adjusted horizons or is copied as-is
      is_c14_file   <- grepl("_C14\\.csv$", file)
      should_modify <- (!has_c14_file) || is_c14_file

      if (should_modify) {

        df <- read.csv(file, stringsAsFactors = FALSE)

        if (!"labID" %in% colnames(df)) {
          df$labID <- NA_integer_
          df <- df[c("labID", setdiff(names(df), "labID"))]
        }
        df$labID <- as.integer(df$labID)

        added_rows_list <- list()
        assigned_labIDs <- integer(0)
        labID_counter   <- if (all(is.na(df$labID))) 1L else max(df$labID, na.rm = TRUE) + 1L

        for (horizon_name in names(adjusted_ages)) {
          adj_df             <- adjusted_ages[[horizon_name]]
          pattern_age        <- paste0(age_record, ".*$")
          adj_df$record_base <- sub(pattern_age, "", adj_df$record)

          if (core_name %in% adj_df$record_base) {
            indices <- which(adj_df$record_base == core_name)

            for (idx in indices) {
              adj_age     <- adj_df$adjusted_age[idx]
              adj_error   <- adj_df$adjusted_error[idx]
              core_depths <- record_data[[core_name]]
              match_row   <- core_depths[core_depths$event == horizon_name, , drop = FALSE]

              if (nrow(match_row) == 0) {
                variant_pattern <- paste0("^", horizon_name, "[a-zA-Z0-9]*$")
                match_row <- core_depths[grepl(variant_pattern, core_depths$event), , drop = FALSE]
              }

              if (nrow(match_row) == 0) {
                warning("No matching depth for horizon '", horizon_name,
                        "' in record '", core_name, "'")
                next
              }

              depth_val <- match_row$depth[1]

              new_row <- as.data.frame(matrix(NA, nrow = 1, ncol = ncol(df)))
              colnames(new_row) <- colnames(df)

              new_row$labID <- labID_counter
              if ("age"   %in% names(df)) new_row$age   <- adj_age
              if ("error" %in% names(df)) new_row$error <- adj_error
              if ("depth" %in% names(df)) new_row$depth <- depth_val
              if ("cc"    %in% names(df)) new_row$cc    <- 0
              if ("d.R"   %in% names(df)) new_row$d.R   <- 0
              if ("d.STD" %in% names(df)) new_row$d.STD <- 0
              if ("t.a"   %in% names(df)) new_row$t.a   <- 33
              if ("t.b"   %in% names(df)) new_row$t.b   <- 34

              added_rows_list[[length(added_rows_list) + 1]] <- new_row
              assigned_labIDs <- c(assigned_labIDs, labID_counter)
              labID_counter   <- labID_counter + 1L
            }
          }
        }

        if (length(added_rows_list) > 0) {
          added_block <- do.call(rbind, added_rows_list)
          for (nm in intersect(names(df), names(added_block))) {
            if (is.integer(df[[nm]]))      added_block[[nm]] <- as.integer(added_block[[nm]])
            else if (is.numeric(df[[nm]])) added_block[[nm]] <- as.numeric(added_block[[nm]])
            else                           added_block[[nm]] <- as.character(added_block[[nm]])
          }
          df_out       <- rbind(df, added_block)
          any_modified <- TRUE
        } else {
          df_out <- df
          cat("No new rows added for", core_name, "- no matching horizons found\n")
        }

        if ("depth" %in% colnames(df_out)) {
          df_out$depth <- as.numeric(df_out$depth)
          df_out       <- df_out[order(df_out$depth, na.last = TRUE), ]
          rownames(df_out) <- NULL
        }

        if (overwrite) {
          output_file <- file
        } else {
          new_basename <- if (synced) {
            if (is_c14_file) {
              paste0(core_name_synced, "_", next_n, "_C14",  ".csv")
            } else {
              paste0(core_name_synced, "_", next_n, ".csv")
            }
          } else {
            paste0(core_name, update_records, ".csv")
          }
          output_file <- file.path(new_folder, new_basename)
        }

        outputs[[length(outputs) + 1]] <- list(
          type            = "write",
          output_file     = output_file,
          df_out          = df_out,
          core_name       = core_name,
          assigned_labIDs = assigned_labIDs
        )

      } else {
        # Pb210 file: copy as-is
        if (!overwrite) {
          pb_basename <- if (synced) {
            paste0(sub("\\.csv$", "", basename(file)), "_", next_n, ".csv")
          } else {
            basename(file)
          }
          output_file <- file.path(new_folder, pb_basename)
          outputs[[length(outputs) + 1]] <- list(
            type        = "copy",
            source_file = file,
            output_file = output_file
          )
        }
      }
    }

    # Only write files if something was actually modified
    if (!any_modified) {
      cat("No changes for", basename(folder), "- skipping folder creation\n\n")
      next
    }

    # Create folder and write all outputs
    if (!overwrite) {
      dir.create(new_folder, showWarnings = FALSE, recursive = TRUE)
    }

    for (out in outputs) {
      if (out$type == "write") {
        write.csv(out$df_out, out$output_file, row.names = FALSE)
        cat("Updated CSV for", out$core_name, "-> saved to:", out$output_file, "\n")

        if (non_synchro && length(out$assigned_labIDs) > 0) {
          df_out       <- out$df_out
          df_out$labID <- as.integer(df_out$labID)
          positions    <- vapply(out$assigned_labIDs, function(id) {
            pos <- which(df_out$labID == id)
            if (length(pos) == 0) NA_integer_ else pos[1]
          }, integer(1))
          for (k in seq_along(out$assigned_labIDs)) {
            cat(sprintf("Added sample: labID = %d at final row = %s in %s\n",
                        out$assigned_labIDs[k], positions[k],
                        basename(out$output_file)))
          }
        }

      } else if (out$type == "copy") {
        file.copy(out$source_file, out$output_file, overwrite = TRUE)
        cat("Copied Pb210 file as-is:", basename(out$source_file),
            "->", out$output_file, "\n")
      }
      cat("\n")
    }
  }
}
