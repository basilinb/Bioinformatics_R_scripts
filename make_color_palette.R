# ── Dependencies ────────────────────────────────────────────────────────────
# install.packages(c("paletteer", "pals", "dplyr", "colorspace", "farver"))

#' Global color registry for stable variable→color assignment across plots
#' Stored in the session environment; optionally persisted to disk.
.color_registry <- new.env(parent = emptyenv())

# Known colorblind-safe qualitative palettes (curated list)
# Sources: colorbrewer2.org, Bang Wong (Nature Methods 2011), Okabe-Ito
.colorblind_safe_palettes <- c(
  "RColorBrewer::Set2",
  "RColorBrewer::Paired",
  "RColorBrewer::Dark2",
  "ggthemes::colorblind",
  "ggthemes::few_Dark",
  "ggthemes::few_Medium",
  "ggsci::category10_d3",
  "ochRe::ngunnawal",
  "colorBlindness::paletteMartin",
  "colorBlindness::Blue2DarkRed12Steps"
)

# ── Scoring helper ───────────────────────────────────────────────────────────

#' Score a palette by perceptual distinctiveness (mean min CIEDE2000 distance)
#' Higher = more visually distinct colors = better
.score_palette_distinctiveness <- function(hex_colors) {
  # Convert hex → LAB color space for perceptual distance
  lab <- farver::decode_colour(hex_colors, to = "lab")

  if (nrow(lab) < 2) {
    return(0)
  }

  # Pairwise CIEDE2000 distances
  dist_matrix <- farver::compare_colour(
    lab,
    lab,
    from_space = "lab",
    metric = "cie2000"
  )
  diag(dist_matrix) <- NA

  # Mean of each color's nearest neighbor distance
  mean(apply(dist_matrix, 1, min, na.rm = TRUE))
}

# ── Registry helpers ─────────────────────────────────────────────────────────

#' Save the registry to disk for cross-session persistence
#' @param path File path for the RDS cache (default: ~/.color_registry.rds)
save_color_registry <- function(path = "~/.color_registry.rds") {
  registry_list <- as.list(.color_registry)
  saveRDS(registry_list, path)
  message(sprintf(
    "Color registry saved to %s (%d entries)",
    path,
    length(registry_list)
  ))
  invisible(registry_list)
}

#' Load a previously saved registry from disk
#' @param path File path for the RDS cache
load_color_registry <- function(path = "~/.color_registry.rds") {
  if (!file.exists(path)) {
    message("No registry file found at: ", path)
    return(invisible(NULL))
  }
  registry_list <- readRDS(path)
  for (nm in names(registry_list)) {
    assign(nm, registry_list[[nm]], envir = .color_registry)
  }
  message(sprintf(
    "Color registry loaded from %s (%d entries)",
    path,
    length(registry_list)
  ))
  invisible(registry_list)
}

#' View all registered variable sets
view_color_registry <- function() {
  keys <- ls(.color_registry)
  if (length(keys) == 0) {
    message("Registry is empty.")
    return(invisible(NULL))
  }
  for (k in keys) {
    colors <- get(k, envir = .color_registry)
    cat(sprintf("\n[%s]\n", k))
    print(colors)
  }
  invisible(as.list(.color_registry))
}

#' Remove a specific entry from the registry (forces re-assignment next call)
forget_color_registry <- function(key) {
  if (exists(key, envir = .color_registry)) {
    rm(list = key, envir = .color_registry)
    message("Removed '", key, "' from registry.")
  } else {
    message("Key '", key, "' not found in registry.")
  }
}

# ── Main function ─────────────────────────────────────────────────────────────

#' Make a stable, scored named color palette for categorical variables
#'
#' @param variables       Character/factor vector. Unique values are extracted
#'                        in first-appearance order (or factor level order).
#' @param palette         Optional "package::palette" string to pin a specific
#'                        palette. Bypasses scoring. NULL = auto-select.
#' @param colorblind_safe Prefer colorblind-safe palettes? Default TRUE.
#'                        If TRUE and no safe palette covers n, falls back to
#'                        full candidate pool with a warning.
#' @param registry_key    Optional string key for the color registry. If
#'                        supplied, previously assigned colors are reused and
#'                        only new variables get new colors. Use the same key
#'                        across batches/plots for stable assignment.
#' @param verbose         Print which palette was chosen. Default TRUE.
#'
#' @return Named character vector of hex colors, one per unique variable value.
#'
#' @examples
#' # Basic auto-selection
#' make_color_palette(c("Control", "24hr", "48hr", "96hr", "10Day"))
#'
#' # Stable across batches — same variables always get same colors
#' cols <- make_color_palette(timepoints, registry_key = "timepoint")
#'
#' # Pin a specific palette
#' make_color_palette(timepoints, palette = "wesanderson::AsteroidCity1")
#'
#' # ggplot2 usage
#' ggplot(df, aes(x, y, color = TimePoint)) +
#'   geom_point() +
#'   scale_color_manual(values = make_color_palette(df$TimePoint,
#'                                                  registry_key = "timepoint"))

make_color_palette <- function(
  variables,
  palette = NULL,
  colorblind_safe = TRUE,
  registry_key = NULL,
  verbose = TRUE
) {
  stopifnot(is.character(variables) || is.factor(variables))

  # ── 1. Extract unique levels ───────────────────────────────────────────────
  if (is.factor(variables)) {
    # Respect factor levels (e.g. from arrange + factor())
    levels_ordered <- levels(droplevels(as.factor(variables)))
    levels_ordered <- levels_ordered[
      levels_ordered %in% as.character(variables)
    ]
  } else {
    levels_ordered <- unique(as.character(variables))
  }
  n <- length(levels_ordered)

  # ── 2. Hard cap ────────────────────────────────────────────────────────────
  if (n > 22) {
    stop(sprintf(
      paste0(
        "You have %d unique values. Maximum supported is 22 (pals::kelly).\n",
        "Consider: collapsing rare categories, using a continuous scale,\n",
        "or faceting your plot."
      ),
      n
    ))
  }

  # ── 3. Registry lookup — reuse existing assignments ───────────────────────
  new_levels <- levels_ordered # will be narrowed if registry hits exist
  existing_colors <- character(0)

  if (!is.null(registry_key)) {
    reg_key_clean <- gsub("[^a-zA-Z0-9_]", "_", registry_key)

    if (exists(reg_key_clean, envir = .color_registry)) {
      cached <- get(reg_key_clean, envir = .color_registry)
      already_assigned <- intersect(levels_ordered, names(cached))

      if (length(already_assigned) > 0) {
        existing_colors <- cached[already_assigned]
        new_levels <- setdiff(levels_ordered, already_assigned)

        if (verbose && length(new_levels) == 0) {
          message(
            "All colors retrieved from registry (key: '",
            registry_key,
            "')"
          )
          return(existing_colors[levels_ordered])
        }
        if (verbose) {
          message(sprintf(
            "Registry hit (key: '%s'): %d reused, %d new variables to assign.",
            registry_key,
            length(already_assigned),
            length(new_levels)
          ))
        }
      }
    }
  }

  n_new <- length(new_levels)

  # ── 4. Explicit palette ────────────────────────────────────────────────────
  if (!is.null(palette)) {
    if (!grepl("::", palette)) {
      stop('palette must be "package::name".')
    }
    colors_raw <- tryCatch(
      as.character(paletteer::paletteer_d(palette)),
      error = function(e) stop("Palette not found: ", palette, "\n", e$message)
    )
    if (length(colors_raw) < n_new) {
      stop(sprintf(
        "Palette '%s' has %d colors but %d are needed.",
        palette,
        length(colors_raw),
        n_new
      ))
    }
    new_colors <- setNames(colors_raw[seq_len(n_new)], new_levels)
    if (verbose) message("Using palette: ", palette)
  } else {
    # ── 5. pals::kelly fallback for > 10 new variables ─────────────────────
    if (n_new > 10) {
      warning(sprintf(
        paste0(
          "%d new unique values exceeds the qualitative palette limit (10).\n",
          "Falling back to pals::kelly (22 colors). ",
          "Consider reviewing category structure."
        ),
        n_new
      ))
      kelly_colors <- pals::kelly(22)[-1] # drop first color (white/near-white)
      # Avoid reusing colors already in the registry
      used_colors <- unname(existing_colors)
      avail_colors <- kelly_colors[!kelly_colors %in% used_colors]
      new_colors <- setNames(avail_colors[seq_len(n_new)], new_levels)
      if (verbose) message("Using palette: pals::kelly")
    } else {
      # ── 6. Score and select best qualitative palette ──────────────────────

      candidates <- paletteer::palettes_d_names |>
        dplyr::filter(type == "qualitative", length >= n_new)

      if (nrow(candidates) == 0) {
        stop(sprintf("No qualitative palette found with >= %d colors.", n_new))
      }

      # 6a. Prefer colorblind-safe palettes if requested
      if (colorblind_safe) {
        safe_key <- paste0(candidates$package, "::", candidates$palette)
        safe_candidates <- candidates[safe_key %in% .colorblind_safe_palettes, ]

        if (nrow(safe_candidates) == 0) {
          warning(paste0(
            "No colorblind-safe palette covers ",
            n_new,
            " colors. ",
            "Falling back to full palette pool."
          ))
        } else {
          candidates <- safe_candidates
        }
      }

      # 6b. Score each candidate by perceptual distinctiveness
      # (only score first n_new colors — what we'd actually use)
      if (verbose) {
        message(sprintf(
          "Scoring %d candidate palettes...",
          nrow(candidates)
        ))
      }

      candidates$score <- vapply(
        seq_len(nrow(candidates)),
        function(i) {
          pal_str <- sprintf(
            "%s::%s",
            candidates$package[i],
            candidates$palette[i]
          )
          tryCatch(
            {
              cols <- as.character(paletteer::paletteer_d(pal_str))[seq_len(
                n_new
              )]
              .score_palette_distinctiveness(cols)
            },
            error = function(e) -1
          )
        },
        numeric(1)
      )

      # 6c. Penalize oversized palettes (prefer palette length close to n_new)
      # Small penalty: subtract 0.2 * (extra colors / 10), keeps scoring primary
      candidates$score <- candidates$score -
        0.2 * pmax(0, candidates$length - n_new) / 10

      best <- candidates[which.max(candidates$score), ]
      palette_str <- sprintf("%s::%s", best$package, best$palette)

      colors_raw <- as.character(paletteer::paletteer_d(palette_str))
      new_colors <- setNames(colors_raw[seq_len(n_new)], new_levels)

      if (verbose) {
        message(sprintf(
          "Selected palette: %s (score: %.2f, %d of %d colors used)",
          palette_str,
          best$score,
          n_new,
          best$length
        ))
      }
    }
  }

  # ── 7. Merge new assignments with registry cache ────────────────────────────
  all_colors <- c(existing_colors, new_colors)[levels_ordered]

  # Update registry
  if (!is.null(registry_key)) {
    reg_key_clean <- gsub("[^a-zA-Z0-9_]", "_", registry_key)
    updated_cache <- if (exists(reg_key_clean, envir = .color_registry)) {
      existing_cache <- get(reg_key_clean, envir = .color_registry)
      # Merge: existing entries take priority, new ones added
      c(
        existing_cache,
        new_colors[!names(new_colors) %in% names(existing_cache)]
      )
    } else {
      all_colors
    }
    assign(reg_key_clean, updated_cache, envir = .color_registry)
  }

  all_colors
}

# ── ggplot2 scale wrappers ───────────────────────────────────────────────────

#' ggplot2 color scale using make_color_palette
#' Drop-in replacement for scale_color_manual
scale_color_auto <- function(
  variables,
  palette = NULL,
  colorblind_safe = TRUE,
  registry_key = NULL,
  ...
) {
  ggplot2::scale_color_manual(
    values = make_color_palette(
      variables,
      palette = palette,
      colorblind_safe = colorblind_safe,
      registry_key = registry_key,
      verbose = FALSE
    ),
    ...
  )
}

#' ggplot2 fill scale using make_color_palette
scale_fill_auto <- function(
  variables,
  palette = NULL,
  colorblind_safe = TRUE,
  registry_key = NULL,
  ...
) {
  ggplot2::scale_fill_manual(
    values = make_color_palette(
      variables,
      palette = palette,
      colorblind_safe = colorblind_safe,
      registry_key = registry_key,
      verbose = FALSE
    ),
    ...
  )
}
