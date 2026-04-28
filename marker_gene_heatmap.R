#' Marker gene / gene-set heatmap with optional group splitting and within-group clustering
#'
#' Accepts either:
#'   (a) an EList (`elist`) — original behaviour, or
#'   (b) a wide module tibble (`module_df`) where the first column (`geneSet`) holds
#'       row labels and all remaining columns are sample IDs (libids).
#'
#' @param genes        For EList mode: character vector of HGNC symbols to plot.
#'                     For module_df mode: optional character vector of geneSet names
#'                     to subset; NULL plots all rows.
#' @param elist        An EList (e.g. from voom) with $E, $genes (needs `hgnc_symbol`,
#'                     `geneName`), and $targets (needs `libid`).
#' @param module_df    Wide tibble/data.frame: first column = geneSet (row labels),
#'                     remaining columns = samples (column names are libids).
#' @param targets      Metadata data.frame required when using `module_df`. Must
#'                     contain a column matching the sample-ID column names of
#'                     `module_df` (see `libid_col`).
#' @param libid_col    Name of the sample-ID column in `targets`. Default "libid".
#' @param split_by     Column name in targets to split columns by (e.g.
#'                     "asthma_status"). NULL → no split.
#' @param split_levels Optional character vector giving the desired factor order for
#'                     `split_by`. Inferred from data if NULL.
#' @param n_clusters   How many clusters to cut per group.
#'                     - NULL          : no within-group clustering
#'                     - single int    : same k for every group
#'                     - named int vec : per-group k, e.g. c(Asthma = 4, Non_Asthma = 3)
#' @param anno_vars    Named list of extra columns in targets to add as column
#'                     annotations, e.g. list(T2 = "t2_1sd", Eosinophils = "EOS_cat").
#' @param anno_colors  Named list of color vectors for annotations (passed to
#'                     `HeatmapAnnotation(col = ...)`).
#' @param cluster_colors  Optional named character vector of colours for the
#'                        Donor_Clusters annotation. Auto-generated if NULL.
#' @param hclust_method   Method passed to `hclust()`. Default "ward.D2".
#' @param scale_rows   Logical. Apply row-wise z-score? Default TRUE. Set FALSE if
#'                     values are already on a comparable scale.
#' @param row_title    Title for the row axis. Default "Marker genes".
#' @param output_path  Full file path to save (PDF or PNG detected from extension).
#'                     NULL → no file written.
#' @param width,height Plot dimensions in inches. Defaults 9 × 5.
#'
#' @return Invisibly: list(ht, col_meta, cluster_vec, cluster_colors).
marker_gene_heatmap <- function(
  genes = NULL,
  elist = NULL,
  module_df = NULL,
  targets = NULL,
  libid_col = "libid",
  split_by = NULL,
  split_levels = NULL,
  n_clusters = NULL,
  anno_vars = list(),
  anno_colors = list(),
  cluster_colors = NULL,
  hclust_method = "ward.D2",
  scale_rows = TRUE,
  row_title = "Marker genes",
  output_path = NULL,
  width = 9,
  height = 5
) {
  library(ComplexHeatmap)
  library(circlize)

  if (is.null(elist) == is.null(module_df)) {
    stop("Supply exactly one of `elist` or `module_df`.")
  }

  # ══════════════════════════════════════════════════════════════════════════
  # PATH A — EList
  # ══════════════════════════════════════════════════════════════════════════
  if (!is.null(elist)) {
    if (is.null(genes) || length(genes) == 0) {
      stop("`genes` must be a non-empty character vector when using `elist`.")
    }

    gene_map <- elist$genes |>
      dplyr::filter(hgnc_symbol %in% genes) |>
      dplyr::distinct(geneName, hgnc_symbol)

    if (nrow(gene_map) == 0) {
      stop("None of the requested genes were found in elist$genes$hgnc_symbol.")
    }

    expr_raw <- elist$E[gene_map$geneName, , drop = FALSE]

    if (anyDuplicated(gene_map$hgnc_symbol)) {
      idx_list <- split(seq_len(nrow(expr_raw)), gene_map$hgnc_symbol)
      expr_mat <- t(sapply(idx_list, function(i) {
        colMeans(expr_raw[i, , drop = FALSE])
      }))
    } else {
      expr_mat <- expr_raw
      rownames(expr_mat) <- gene_map$hgnc_symbol
    }

    meta <- elist$targets
    join_col <- "libid" # EList targets always uses "libid"

    # ══════════════════════════════════════════════════════════════════════════
    # PATH B — Wide module tibble
    # ══════════════════════════════════════════════════════════════════════════
  } else {
    if (is.null(targets)) {
      stop("`targets` must be supplied when using `module_df`.")
    }
    if (!libid_col %in% colnames(targets)) {
      stop("`targets` has no column named '", libid_col, "' (see `libid_col`).")
    }

    gs_col <- colnames(module_df)[1]
    samp_cols <- setdiff(colnames(module_df), gs_col)

    if (!is.null(genes)) {
      keep <- module_df[[gs_col]] %in% genes
      if (!any(keep)) {
        stop(
          "None of the values in `genes` were found in the '",
          gs_col,
          "' column."
        )
      }
      module_df <- module_df[keep, , drop = FALSE]
    }

    expr_mat <- as.matrix(module_df[, samp_cols, drop = FALSE])
    rownames(expr_mat) <- module_df[[gs_col]]
    storage.mode(expr_mat) <- "double"

    meta <- targets
    join_col <- libid_col # ← use as-is; no renaming, no collision risk
  }

  # ══════════════════════════════════════════════════════════════════════════
  # SHARED PATH
  # ══════════════════════════════════════════════════════════════════════════

  # ── 3. Optional row-wise z-score ──────────────────────────────────────────
  if (scale_rows) {
    expr_z <- t(scale(t(expr_mat)))
    expr_z[is.na(expr_z)] <- 0
  } else {
    expr_z <- expr_mat
  }

  # ── 4. Align metadata to matrix columns ───────────────────────────────────
  col_meta <- meta[
    match(colnames(expr_z), meta[[join_col]]),
    ,
    drop = FALSE
  ]

  # ── 5. Column split factor ────────────────────────────────────────────────
  col_split <- NULL
  if (!is.null(split_by)) {
    if (!split_by %in% colnames(col_meta)) {
      stop("`split_by` column '", split_by, "' not found in targets.")
    }
    grp_vals <- as.character(col_meta[[split_by]])
    lvls <- split_levels %||% unique(grp_vals[!is.na(grp_vals)])
    col_split <- factor(grp_vals, levels = lvls)
  }

  # ── 6. Within-group clustering ────────────────────────────────────────────
  cluster_vec <- NULL
  col_order_vec <- NULL
  do_cluster <- !is.null(n_clusters) && !is.null(col_split)

  if (do_cluster) {
    cluster_vec <- rep(NA_character_, ncol(expr_z))
    names(cluster_vec) <- colnames(expr_z)

    col_order_list <- vector("list", nlevels(col_split))
    names(col_order_list) <- levels(col_split)

    for (grp in levels(col_split)) {
      idx <- which(col_split == grp)

      k <- if (length(n_clusters) > 1 && grp %in% names(n_clusters)) {
        n_clusters[[grp]]
      } else {
        n_clusters[[1]]
      }

      if (length(idx) >= 2) {
        mat_sub <- expr_z[, idx, drop = FALSE]
        hc <- hclust(dist(t(mat_sub)), method = hclust_method)
        cl <- cutree(hc, k = min(k, length(idx)))
        cluster_vec[idx] <- paste0(grp, "_C", cl[colnames(mat_sub)])
        col_order_list[[grp]] <- colnames(mat_sub)[hc$order]
      } else {
        cluster_vec[idx] <- paste0(grp, "_C1")
        col_order_list[[grp]] <- colnames(expr_z)[idx]
      }
    }

    col_order_vec <- unlist(lapply(col_order_list, function(x) {
      match(x, colnames(expr_z))
    }))

    cluster_levels <- unlist(lapply(levels(col_split), function(grp) {
      k <- if (length(n_clusters) > 1 && grp %in% names(n_clusters)) {
        n_clusters[[grp]]
      } else {
        n_clusters[[1]]
      }
      paste0(grp, "_C", seq_len(k))
    }))
    col_meta$Donor_Clusters <- factor(cluster_vec, levels = cluster_levels)

    if (is.null(cluster_colors)) {
      pal <- colorRampPalette(c(
        "#1d1d1d",
        "#ebce2b",
        "#702c8c",
        "#db6917",
        "#96cde6",
        "#ba1c30",
        "#c0bd7f",
        "#7f7e80",
        "#5fa641",
        "#d485b2",
        "#4277b6",
        "#df8461",
        "#463397",
        "#e1a11a",
        "#91218c",
        "#e8e948",
        "#7e1510",
        "#92ae31",
        "#6f340d",
        "#d32b1e",
        "#2b3514"
      ))
      cluster_colors <- setNames(pal(length(cluster_levels)), cluster_levels)
    }
  }

  # ── 7. Build column annotation ────────────────────────────────────────────
  anno_list <- list()
  color_list <- anno_colors

  if (!is.null(col_split)) {
    anno_list[[split_by]] <- col_split
  }

  for (nm in names(anno_vars)) {
    col_nm <- anno_vars[[nm]]
    if (!col_nm %in% colnames(col_meta)) {
      warning("anno_vars column '", col_nm, "' not found in targets — skipped.")
      next
    }
    anno_list[[nm]] <- col_meta[[col_nm]]
  }

  if (do_cluster) {
    anno_list[["Donor_Clusters"]] <- col_meta$Donor_Clusters
    color_list[["Donor_Clusters"]] <- cluster_colors
  }

  ha <- do.call(
    HeatmapAnnotation,
    c(anno_list, list(col = color_list, annotation_name_side = "left"))
  )

  # ── 8. Draw heatmap ───────────────────────────────────────────────────────
  z_label <- if (scale_rows) "z-score" else "value"

  ht <- Heatmap(
    expr_z,
    name = z_label,
    top_annotation = ha,
    show_row_names = TRUE,
    show_column_names = FALSE,
    cluster_rows = TRUE,
    cluster_columns = !do_cluster,
    column_split = col_split,
    column_order = col_order_vec,
    cluster_column_slices = !is.null(col_split),
    use_raster = TRUE,
    row_title = row_title,
    heatmap_legend_param = list(title = z_label)
  )

  grid::grid.newpage()
  draw(ht, heatmap_legend_side = "right")

  # ── 9. Optional file output ───────────────────────────────────────────────
  if (!is.null(output_path)) {
    ext <- tolower(tools::file_ext(output_path))
    if (ext == "pdf") {
      pdf(output_path, width = width, height = height)
    } else if (ext %in% c("png", "tiff", "jpeg", "jpg")) {
      png(output_path, width = width, height = height, units = "in", res = 300)
    } else {
      warning("Unrecognised extension '", ext, "' — writing as PDF.")
      pdf(output_path, width = width, height = height)
    }
    draw(ht, heatmap_legend_side = "right")
    invisible(dev.off())
    message("Saved: ", output_path)
  }

  invisible(list(
    ht = ht,
    col_meta = col_meta,
    cluster_vec = cluster_vec,
    cluster_colors = cluster_colors
  ))
}

`%||%` <- function(a, b) if (!is.null(a)) a else b
