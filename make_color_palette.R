# ── Dependencies ────────────────────────────────────────────────────────────
# install.packages(c("paletteer", "pals", "dplyr", "farver", "ggplot2"))
# Colorblind-safe candidates additionally live in RColorBrewer, ggthemes,
# ggsci, colorBlindness and ochRe; missing ones are skipped during scoring.

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

# Ceiling (CIEDE2000 ΔE) on how much credit a palette earns for being far from
# colors used by other variables. Without a cap, a palette that is wildly
# distant from the avoided colors could outscore one that is internally more
# distinct. ~40 ΔE is already "obviously different" to the eye.
.AVOID_DIST_CAP <- 40

# Largest number of unique values we can color: pals::kelly(22) minus its
# first (near-white) entry.
.MAX_LEVELS <- 21

# ── Color utilities ──────────────────────────────────────────────────────────

#' Validate and canonicalize colors to uppercase hex
#'
#' Everything that gets set-compared (registry entries, avoid pools, palette
#' colors) must pass through here first, or "black" and "#000000" read as
#' different colors.
#'
#' @param x        Character vector of color names or hex codes (NULL allowed).
#' @param arg_name Name to use in error messages.
#' @return Uppercase "#RRGGBB" (or "#RRGGBBAA" when partly transparent).
.normalize_colors <- function(x, arg_name = "colors") {
  if (is.null(x) || length(x) == 0) {
    return(character(0))
  }
  if (!is.character(x)) {
    stop(sprintf(
      "`%s` must be a character vector of color names or hex codes.",
      arg_name
    ))
  }

  ok <- vapply(
    x,
    function(col) {
      !is.na(col) &&
        isTRUE(tryCatch(
          {
            grDevices::col2rgb(col, alpha = TRUE)
            TRUE
          },
          error = function(e) FALSE
        ))
    },
    logical(1),
    USE.NAMES = FALSE
  )

  if (any(!ok)) {
    stop(sprintf(
      paste0(
        "`%s` contains %d value(s) R cannot interpret as a color: %s\n",
        'Use a hex code ("#1B9E77") or a name from colors().'
      ),
      arg_name,
      sum(!ok),
      paste0('"', x[!ok], '"', collapse = ", ")
    ))
  }

  rgba <- grDevices::col2rgb(x, alpha = TRUE)
  out <- ifelse(
    rgba["alpha", ] == 255L,
    grDevices::rgb(
      rgba["red", ],
      rgba["green", ],
      rgba["blue", ],
      maxColorValue = 255
    ),
    grDevices::rgb(
      rgba["red", ],
      rgba["green", ],
      rgba["blue", ],
      rgba["alpha", ],
      maxColorValue = 255
    )
  )
  toupper(unname(out))
}

#' Colors from `pool` that are not in the (already normalized) `avoid` set
.free_colors <- function(pool, avoid) {
  if (length(avoid) == 0) {
    return(pool)
  }
  pool[!.normalize_colors(pool) %in% avoid]
}

#' Greedy farthest-point pick — last resort when the pool cannot supply
#' `n` colors that are unused elsewhere.
#'
#' Repeatedly takes the color whose nearest neighbor (among the avoided colors
#' and the colors already picked) is furthest away. Deterministic.
.farthest_point_pick <- function(pool, n, avoid) {
  pool_lab <- farver::decode_colour(pool, to = "lab")

  min_d <- if (length(avoid) > 0) {
    apply(
      farver::compare_colour(
        pool_lab,
        farver::decode_colour(avoid, to = "lab"),
        from_space = "lab",
        method = "cie2000"
      ),
      1,
      min
    )
  } else {
    rep(Inf, nrow(pool_lab))
  }

  picked <- integer(0)
  for (i in seq_len(n)) {
    candidates <- setdiff(seq_len(nrow(pool_lab)), picked)
    best <- candidates[which.max(min_d[candidates])]
    picked <- c(picked, best)
    min_d <- pmin(
      min_d,
      farver::compare_colour(
        pool_lab,
        pool_lab[best, , drop = FALSE],
        from_space = "lab",
        method = "cie2000"
      )[, 1]
    )
  }
  pool[picked]
}

#' Take `n` colors out of `pool`, steering away from `avoid`
#'
#' With no avoid set this is just the first `n` colors (the historical
#' behavior). With one, avoided colors are dropped and the palette's own order
#' is preserved; only if too few remain do we fall back to a distance-based
#' pick and warn.
.pick_colors <- function(pool, n, avoid = character(0), verbose = TRUE) {
  if (n <= 0) {
    return(character(0))
  }
  if (n > length(pool)) {
    stop(sprintf(
      "Internal error: asked for %d colors from a pool of %d.",
      n,
      length(pool)
    ))
  }
  if (length(avoid) == 0) {
    return(pool[seq_len(n)])
  }

  free <- .free_colors(pool, avoid)
  if (length(free) >= n) {
    return(free[seq_len(n)])
  }

  warning(sprintf(
    paste0(
      "Only %d of %d palette colors are unused by the avoided variables, ",
      "but %d are needed.\n",
      "Falling back to the most perceptually distant colors available — ",
      "some may repeat colors used elsewhere.\n",
      "Consider passing a larger `palette`, or `append_colors`."
    ),
    length(free),
    length(pool),
    n
  ))
  .farthest_point_pick(pool, n, avoid)
}

# ── Scoring helper ───────────────────────────────────────────────────────────

#' Score a palette by perceptual distinctiveness (mean min CIEDE2000 distance)
#' Higher = more visually distinct colors = better
#'
#' @param hex_colors Colors that would actually be used.
#' @param avoid      Optional normalized colors already used by other
#'                   variables. When supplied, a color sitting on top of one of
#'                   these is penalized exactly as hard as one sitting on top of
#'                   another color in the same palette. With `avoid = NULL` this
#'                   is the original mean-nearest-neighbor score, unchanged.
.score_palette_distinctiveness <- function(hex_colors, avoid = NULL) {
  # Convert hex → LAB color space for perceptual distance
  lab <- farver::decode_colour(hex_colors, to = "lab")

  no_avoid <- is.null(avoid) || length(avoid) == 0
  if (nrow(lab) == 0 || (nrow(lab) < 2 && no_avoid)) {
    return(0)
  }

  # Each color's nearest neighbor within the palette
  dist_within <- if (nrow(lab) >= 2) {
    dist_matrix <- farver::compare_colour(
      lab,
      lab,
      from_space = "lab",
      method = "cie2000"
    )
    diag(dist_matrix) <- NA
    apply(dist_matrix, 1, min, na.rm = TRUE)
  } else {
    rep(Inf, nrow(lab))
  }

  if (no_avoid) {
    return(mean(dist_within))
  }

  # Each color's nearest neighbor among colors used by other variables
  dist_avoid <- apply(
    farver::compare_colour(
      lab,
      farver::decode_colour(avoid, to = "lab"),
      from_space = "lab",
      method = "cie2000"
    ),
    1,
    min
  )

  mean(pmin(dist_within, pmin(dist_avoid, .AVOID_DIST_CAP)))
}

# ── Registry helpers ─────────────────────────────────────────────────────────

#' Canonical form of a registry key (single source of truth)
.sanitize_registry_key <- function(key) {
  if (!is.character(key) || length(key) != 1 || is.na(key)) {
    stop("`registry_key` must be a single string.")
  }
  gsub("[^a-zA-Z0-9_]", "_", key)
}

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
#'
#' To re-color only some levels of a key while leaving the rest alone, use
#' `make_color_palette(..., registry_overwrite = TRUE)` instead.
forget_color_registry <- function(key) {
  key_clean <- .sanitize_registry_key(key)
  if (exists(key_clean, envir = .color_registry)) {
    rm(list = key_clean, envir = .color_registry)
    message("Removed '", key, "' from registry.")
  } else {
    message("Key '", key, "' not found in registry.")
  }
}

#' Gather colors to steer away from
#'
#' @param avoid_keys   Registry key(s) whose colors are off limits, or "all"
#'                     for every key other than `own_key`.
#' @param avoid_colors Explicit colors to steer away from.
#' @param own_key      Sanitized key of the variable being assigned (never
#'                     avoided against itself here — the caller handles its own
#'                     cached colors).
#' @return Normalized, deduplicated character vector.
.collect_avoid_colors <- function(
  avoid_keys = NULL,
  avoid_colors = NULL,
  own_key = NULL
) {
  from_registry <- character(0)

  if (!is.null(avoid_keys)) {
    if (identical(as.character(avoid_keys), "all")) {
      keys <- setdiff(ls(.color_registry), own_key)
    } else {
      keys <- unique(vapply(
        as.character(avoid_keys),
        .sanitize_registry_key,
        character(1),
        USE.NAMES = FALSE
      ))
      found <- vapply(keys, exists, logical(1), envir = .color_registry)
      if (any(!found)) {
        warning(sprintf(
          paste0(
            "avoid_keys not found in the registry: %s\n",
            "Nothing will be avoided for these. Check the spelling, or call ",
            "load_color_registry() first."
          ),
          paste0("'", keys[!found], "'", collapse = ", ")
        ))
      }
      keys <- setdiff(keys[found], own_key)
    }

    for (k in keys) {
      from_registry <- c(
        from_registry,
        unname(get(k, envir = .color_registry))
      )
    }
  }

  unique(c(
    .normalize_colors(from_registry, "registry colors"),
    .normalize_colors(avoid_colors, "avoid_colors")
  ))
}

#' Resolve the `palette` / `append_colors` arguments into a color pool
#'
#' @return List with `colors` (NULL when auto-selection is wanted), `appended`
#'         (normalized `append_colors`) and `label` for messages.
.resolve_palette <- function(palette = NULL, append_colors = NULL) {
  appended <- .normalize_colors(append_colors, "append_colors")

  if (is.null(palette)) {
    return(list(colors = NULL, appended = appended, label = NULL))
  }
  if (!is.character(palette)) {
    stop(paste0(
      '`palette` must be a "package::name" string (e.g. ',
      '"wesanderson::AsteroidCity1") or a character vector of colors.'
    ))
  }

  if (length(palette) == 1 && grepl("::", palette)) {
    colors_raw <- tryCatch(
      as.character(paletteer::paletteer_d(palette)),
      error = function(e) stop("Palette not found: ", palette, "\n", e$message)
    )
    label <- palette
  } else {
    colors_raw <- .normalize_colors(palette, "palette")
    label <- sprintf("custom palette (%d colors)", length(colors_raw))
  }

  if (length(appended) > 0) {
    label <- sprintf("%s + %d appended", label, length(appended))
  }

  # `colors` is the palette's own contribution only — appended colors are
  # reserved for the trailing levels by the caller, so drop any the palette
  # already contains rather than offering them twice.
  list(
    colors = unique(colors_raw),
    appended = appended,
    label = label
  )
}

# ── Main function ─────────────────────────────────────────────────────────────

#' Make a stable, scored named color palette for categorical variables
#'
#' @param variables       Character/factor vector. Unique values are extracted
#'                        in first-appearance order (or factor level order).
#' @param palette         Optional palette to use instead of auto-selection.
#'                        Either a "package::palette" string (resolved through
#'                        paletteer) or a character vector of colors / hex codes
#'                        of your own. NULL = auto-select by scoring.
#' @param append_colors   Optional colors appended to the end of the palette
#'                        pool, e.g. c("black", "grey") for control / reference
#'                        levels. Works with a pinned palette (extends it) and
#'                        with auto-selection (that many fewer colors are taken
#'                        from the scored palette, and the appended colors are
#'                        themselves avoided during scoring).
#' @param colorblind_safe Prefer colorblind-safe palettes? Default TRUE.
#'                        If TRUE and no safe palette covers n, falls back to
#'                        full candidate pool with a warning.
#' @param registry_key    Optional string key for the color registry. If
#'                        supplied, previously assigned colors are reused and
#'                        only new variables get new colors. Use the same key
#'                        across batches/plots for stable assignment.
#' @param avoid_keys      Optional registry key(s) whose colors this variable
#'                        should steer away from, so several variables in one
#'                        document do not share colors. Use "all" for every
#'                        other key currently in the registry. Avoidance is
#'                        folded into palette scoring, so it influences which
#'                        palette is chosen — not just which colors are taken
#'                        from it.
#' @param avoid_colors    Optional explicit colors to steer away from, in
#'                        addition to `avoid_keys`.
#' @param registry_overwrite Reassign the levels passed in this call even if
#'                        they are already in the registry? Default FALSE
#'                        (cached colors win). Other levels stored under the
#'                        key are left untouched — use forget_color_registry()
#'                        to reset a key entirely.
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
#' # Several variables in one figure, none sharing colors
#' tp_cols <- make_color_palette(df$TimePoint, registry_key = "timepoint")
#' sex_cols <- make_color_palette(
#'   df$Sex,
#'   registry_key = "sex",
#'   avoid_keys = "timepoint"
#' )
#' trt_cols <- make_color_palette(
#'   df$Treatment,
#'   registry_key = "treatment",
#'   avoid_keys = "all"
#' )
#'
#' # Pin a specific palette, extended with colors of your own
#' make_color_palette(
#'   timepoints,
#'   palette = "wesanderson::AsteroidCity1",
#'   append_colors = c("black", "grey")
#' )
#'
#' # Bring your own palette entirely
#' make_color_palette(
#'   df$Treatment,
#'   palette = c("#E1BD6D", "#EABE94", "black", "grey40"),
#'   registry_key = "treatment"
#' )
#'
#' # Re-color just these levels, keeping the rest of the key intact
#' make_color_palette(
#'   df$Treatment,
#'   palette = "RColorBrewer::Dark2",
#'   registry_key = "treatment",
#'   registry_overwrite = TRUE
#' )
#'
#' # ggplot2 usage
#' ggplot(df, aes(x, y, color = TimePoint)) +
#'   geom_point() +
#'   scale_color_manual(values = make_color_palette(df$TimePoint,
#'                                                  registry_key = "timepoint"))

make_color_palette <- function(
  variables,
  palette = NULL,
  append_colors = NULL,
  colorblind_safe = TRUE,
  registry_key = NULL,
  avoid_keys = NULL,
  avoid_colors = NULL,
  registry_overwrite = FALSE,
  verbose = TRUE
) {
  stopifnot(is.character(variables) || is.factor(variables))
  stopifnot(is.logical(registry_overwrite), length(registry_overwrite) == 1)

  if (registry_overwrite && is.null(registry_key)) {
    warning(
      "registry_overwrite = TRUE has no effect without a registry_key; ",
      "nothing is cached to overwrite."
    )
  }

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
  if (n > .MAX_LEVELS) {
    stop(sprintf(
      paste0(
        "You have %d unique values. Maximum supported is %d ",
        "(pals::kelly, minus its near-white first color).\n",
        "Consider: collapsing rare categories, using a continuous scale,\n",
        "or faceting your plot."
      ),
      n,
      .MAX_LEVELS
    ))
  }

  # ── 3. Registry lookup — reuse existing assignments ───────────────────────
  new_levels <- levels_ordered # will be narrowed if registry hits exist
  existing_colors <- character(0) # reused colors, keyed by level
  cached <- character(0) # everything already stored under this key
  reg_key_clean <- if (is.null(registry_key)) {
    NULL
  } else {
    .sanitize_registry_key(registry_key)
  }

  if (!is.null(reg_key_clean) && exists(reg_key_clean, envir = .color_registry)) {
    cached <- get(reg_key_clean, envir = .color_registry)
    already_assigned <- intersect(levels_ordered, names(cached))

    if (registry_overwrite) {
      if (verbose && length(already_assigned) > 0) {
        message(sprintf(
          "registry_overwrite: reassigning %d level(s) under key '%s' (%d other level(s) kept).",
          length(already_assigned),
          registry_key,
          length(setdiff(names(cached), levels_ordered))
        ))
      }
    } else if (length(already_assigned) > 0) {
      existing_colors <- cached[already_assigned]
      new_levels <- setdiff(levels_ordered, already_assigned)

      if (length(new_levels) == 0) {
        if (verbose) {
          message(
            "All colors retrieved from registry (key: '",
            registry_key,
            "')"
          )
        }
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

  n_new <- length(new_levels)

  # ── 4. Build the set of colors to steer away from ─────────────────────────
  # Other variables (opt-in) plus every color already stored under this key,
  # so incrementally adding a level cannot duplicate one of its siblings.
  resolved <- .resolve_palette(palette, append_colors)
  appended <- resolved$appended

  avoid <- unique(c(
    .collect_avoid_colors(avoid_keys, avoid_colors, own_key = reg_key_clean),
    .normalize_colors(unname(cached), "registry colors"),
    appended
  ))

  if (verbose && length(avoid) > 0) {
    message(sprintf("Steering away from %d already-used color(s).", length(avoid)))
  }

  # Appended colors are reserved for the tail of the assignment, so the
  # auto-selected palette only needs to cover the rest.
  n_pal <- max(0, n_new - length(appended))

  # ── 5. append_colors alone covers every level ─────────────────────────────
  if (n_pal == 0) {
    new_colors <- setNames(appended[seq_len(n_new)], new_levels)
    if (verbose) {
      message(sprintf("Using append_colors only (%d colors).", n_new))
    }
  } else if (!is.null(resolved$colors)) {
    # ── 6. Explicit palette (pinned string or your own vector) ──────────────
    pool <- resolved$colors
    if (length(pool) < n_pal) {
      stop(sprintf(
        paste0(
          "Palette %s supplies %d colors but %d are needed.\n",
          "Add more with `append_colors`, or pick a longer palette."
        ),
        resolved$label,
        length(pool) + length(appended),
        n_new
      ))
    }
    new_colors <- setNames(
      c(.pick_colors(pool, n_pal, avoid, verbose), appended),
      new_levels
    )
    if (verbose) message("Using palette: ", resolved$label)
  } else {
    # ── 7. pals::kelly fallback for > 10 new variables ─────────────────────
    if (n_pal > 10) {
      warning(sprintf(
        paste0(
          "%d new unique values exceeds the qualitative palette limit (10).\n",
          "Falling back to pals::kelly (%d colors). ",
          "Consider reviewing category structure."
        ),
        n_pal,
        .MAX_LEVELS
      ))
      kelly_colors <- pals::kelly(22)[-1] # drop first color (white/near-white)
      new_colors <- setNames(
        c(.pick_colors(kelly_colors, n_pal, avoid, verbose), appended),
        new_levels
      )
      if (verbose) message("Using palette: pals::kelly")
    } else {
      # ── 8. Score and select best qualitative palette ──────────────────────

      candidates <- paletteer::palettes_d_names |>
        # `type` and `length` are columns of palettes_d_names, resolved by the
        # data mask (not base::length). `.data$` cannot be used here: this file
        # is sourced, not a package, so the pronoun has no import.
        dplyr::filter(type == "qualitative", length >= n_pal)

      if (nrow(candidates) == 0) {
        stop(sprintf("No qualitative palette found with >= %d colors.", n_pal))
      }

      # 8a. Prefer colorblind-safe palettes if requested
      if (colorblind_safe) {
        safe_key <- paste0(candidates$package, "::", candidates$palette)
        safe_candidates <- candidates[safe_key %in% .colorblind_safe_palettes, ]

        if (nrow(safe_candidates) == 0) {
          warning(paste0(
            "No colorblind-safe palette covers ",
            n_pal,
            " colors. ",
            "Falling back to full palette pool."
          ))
        } else {
          candidates <- safe_candidates
        }
      }

      # 8b. Score each candidate by perceptual distinctiveness — of the colors
      # we would actually use, i.e. after dropping any that are already taken
      # by an avoided variable. A candidate that cannot supply n_pal free
      # colors scores -1 (deprioritized, but still usable as a last resort).
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
              all_cols <- as.character(paletteer::paletteer_d(pal_str))
              free <- .free_colors(all_cols, avoid)
              if (length(free) < n_pal) {
                -1
              } else {
                .score_palette_distinctiveness(free[seq_len(n_pal)], avoid)
              }
            },
            error = function(e) -1
          )
        },
        numeric(1)
      )

      # 8c. Penalize oversized palettes (prefer palette length close to n_pal)
      # Small penalty: subtract 0.2 * (extra colors / 10), keeps scoring primary
      candidates$score <- candidates$score -
        0.2 * pmax(0, candidates$length - n_pal) / 10

      best <- candidates[which.max(candidates$score), ]
      palette_str <- sprintf("%s::%s", best$package, best$palette)

      colors_raw <- as.character(paletteer::paletteer_d(palette_str))
      new_colors <- setNames(
        c(.pick_colors(colors_raw, n_pal, avoid, verbose), appended),
        new_levels
      )

      if (verbose) {
        message(sprintf(
          "Selected palette: %s (score: %.2f, %d of %d colors used)",
          palette_str,
          best$score,
          n_pal,
          best$length
        ))
      }
    }
  }

  # ── 9. Merge new assignments with registry cache ───────────────────────────
  all_colors <- c(existing_colors, new_colors)[levels_ordered]

  # Update registry
  if (!is.null(reg_key_clean)) {
    updated_cache <- if (length(cached) > 0) {
      if (registry_overwrite) {
        # New assignments win for the levels in this call; levels stored under
        # the key but not passed in here carry through untouched.
        c(new_colors, cached[!names(cached) %in% names(new_colors)])
      } else {
        # Merge: existing entries take priority, new ones added
        c(cached, new_colors[!names(new_colors) %in% names(cached)])
      }
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
  append_colors = NULL,
  colorblind_safe = TRUE,
  registry_key = NULL,
  avoid_keys = NULL,
  avoid_colors = NULL,
  registry_overwrite = FALSE,
  ...
) {
  ggplot2::scale_color_manual(
    values = make_color_palette(
      variables,
      palette = palette,
      append_colors = append_colors,
      colorblind_safe = colorblind_safe,
      registry_key = registry_key,
      avoid_keys = avoid_keys,
      avoid_colors = avoid_colors,
      registry_overwrite = registry_overwrite,
      verbose = FALSE
    ),
    ...
  )
}

#' ggplot2 fill scale using make_color_palette
scale_fill_auto <- function(
  variables,
  palette = NULL,
  append_colors = NULL,
  colorblind_safe = TRUE,
  registry_key = NULL,
  avoid_keys = NULL,
  avoid_colors = NULL,
  registry_overwrite = FALSE,
  ...
) {
  ggplot2::scale_fill_manual(
    values = make_color_palette(
      variables,
      palette = palette,
      append_colors = append_colors,
      colorblind_safe = colorblind_safe,
      registry_key = registry_key,
      avoid_keys = avoid_keys,
      avoid_colors = avoid_colors,
      registry_overwrite = registry_overwrite,
      verbose = FALSE
    ),
    ...
  )
}
