volcano_funcxn <- function(
  limma_output,
  coeff = NULL,
  fdr = 0.01,
  n_top_gene = 10,
  add_genes_int = NULL,
  gene_col = "hgnc_symbol"
) {
  required_cols <- c("gene", "volcano_comp", "estimate", "FDR", gene_col)
  missing_cols <- setdiff(required_cols, names(limma_output))

  if (length(missing_cols) > 0) {
    stop(
      "Missing required column(s): ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  required_pkgs <- c("dplyr", "ggplot2", "ggrepel", "tibble")
  missing_pkgs <- required_pkgs[
    !vapply(required_pkgs, requireNamespace, logical(1), quietly = TRUE)
  ]

  if (length(missing_pkgs) > 0) {
    stop(
      "Install required package(s): ",
      paste(missing_pkgs, collapse = ", "),
      call. = FALSE
    )
  }

  if (is.null(coeff)) {
    coeff <- unique(limma_output$volcano_comp)
  }

  options(ggrepel.max.overlaps = Inf)

  make_label_df <- function(df) {
    df |>
      dplyr::mutate(delabel = .data[[gene_col]]) |>
      dplyr::select(gene, delabel)
  }

  vol_plt_list <- vector("list", length(coeff))
  names(vol_plt_list) <- coeff

  for (i in coeff) {
    pval_coef <- limma_output |>
      dplyr::filter(volcano_comp == i)

    up_reg <- pval_coef |>
      dplyr::filter(estimate > 0, FDR < fdr) |>
      dplyr::arrange(FDR) |>
      utils::head(n_top_gene) |>
      make_label_df()

    down_reg <- pval_coef |>
      dplyr::filter(estimate < 0, FDR < fdr) |>
      dplyr::arrange(FDR) |>
      utils::head(n_top_gene) |>
      make_label_df()

    add_gene_reg <- if (is.null(add_genes_int)) {
      tibble::tibble(gene = character(), delabel = character())
    } else {
      pval_coef |>
        dplyr::filter(.data[[gene_col]] %in% add_genes_int) |>
        dplyr::arrange(FDR) |>
        utils::head(n_top_gene) |>
        make_label_df()
    }

    label_df <- dplyr::bind_rows(up_reg, down_reg, add_gene_reg) |>
      dplyr::distinct(gene, .keep_all = TRUE)

    pval_plot <- pval_coef |>
      dplyr::left_join(label_df, by = "gene") |>
      dplyr::mutate(
        FDR = dplyr::if_else(FDR <= 0, .Machine$double.xmin, FDR),
        FC.group = dplyr::case_when(
          FDR <= fdr & estimate > 0 ~ "up",
          FDR <= fdr & estimate < 0 ~ "down",
          !is.na(delabel) ~ "Gene of Interest(NS)",
          TRUE ~ "Not significant"
        ),
        FC.group = factor(
          FC.group,
          levels = c("Not significant", "Gene of Interest(NS)", "down", "up")
        )
      ) |>
      dplyr::arrange(FC.group)

    cols <- c(
      "up" = "red",
      "down" = "blue",
      "Gene of Interest(NS)" = "steelblue4",
      "Not significant" = "grey"
    )

    plt_vol <- ggplot2::ggplot(
      data = pval_plot,
      ggplot2::aes(x = estimate, y = -log10(FDR), color = FC.group)
    ) +
      ggplot2::geom_point(size = 3, alpha = 0.7) +
      ggplot2::theme_classic() +
      ggrepel::geom_text_repel(
        ggplot2::aes(label = delabel),
        color = "black",
        na.rm = TRUE
      ) +
      ggplot2::scale_color_manual(values = cols, drop = FALSE) +
      ggplot2::geom_vline(xintercept = 0, color = "black", linetype = 2) +
      ggplot2::geom_hline(
        yintercept = -log10(fdr),
        color = "black",
        linetype = 2
      ) +
      ggplot2::labs(title = i, x = "estimate", y = "-log10(FDR)", color = NULL)

    vol_plt_list[[i]] <- plt_vol
  }

  vol_plt_list
}
