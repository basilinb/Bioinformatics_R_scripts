#' Fit mediator and outcome mixed-effects models across multiple genes
#'
#' For each gene, fits two \code{lmerTest::lmer} models used in causal
#' mediation analysis and extracts fixed-effect coefficients into a tidy
#' data frame with FDR-corrected p-values.
#'
#' @details
#' Two models are fit per gene:
#' \itemize{
#'   \item \strong{Mediator model}: \code{mediator ~ treatment + covariates + (1 | random_effect)}
#'     Tests whether treatment changes gene expression (the mediator).
#'   \item \strong{Outcome model}: \code{outcome ~ treatment + mediator + covariates + (1 | random_effect)}
#'     Tests whether the outcome depends on treatment and gene expression,
#'     adjusting for one another.
#' }
#' Both models are fit with \code{lmerTest::lmer()} to obtain Satterthwaite
#' p-values from \code{summary()}.
#'
#' FDR correction (Benjamini-Hochberg) is applied within each
#' \code{model × term} combination so that, e.g., the treatment coefficient
#' in the mediator model is corrected across genes separately from the
#' treatment coefficient in the outcome model.
#'
#' @param dat Data frame containing gene expression, outcome, treatment, and
#'   covariate columns. Must contain one row per gene-sample combination.
#' @param genes Character vector of gene names to analyse. If \code{NULL}
#'   (default), all unique values in \code{gene_col} are used.
#' @param treat_col String. Column name for the treatment variable.
#' @param mediator_col String. Column name for gene expression (the mediator).
#' @param outcome_col String. Column name for the outcome variable.
#' @param covariates Character vector of covariate column names to adjust for.
#' @param random_effect String. Column name for the random effect grouping
#'   variable. Default \code{"donorId"}.
#' @param gene_col String. Column name identifying genes. Default
#'   \code{"hgnc_symbol"}.
#' @param ref_level String. Reference (control) level of the treatment factor.
#' @param treat_levels Character vector of treatment levels to include.
#'   Combined with \code{ref_level} to set factor levels.
#' @param min_obs Integer. Minimum observations required per gene after
#'   dropping NAs. Default 10.
#'
#' @return A tibble with one row per gene × model × fixed-effect term,
#'   containing:
#' \describe{
#'   \item{gene}{Gene symbol}
#'   \item{model}{Either \code{"mediator"} or \code{"outcome"}}
#'   \item{model_formula}{The formula used (as a string)}
#'   \item{term}{Fixed-effect term name from the model summary}
#'   \item{estimate}{Coefficient estimate}
#'   \item{std_error}{Standard error}
#'   \item{df}{Satterthwaite degrees of freedom}
#'   \item{t_value}{t-statistic}
#'   \item{p_value}{Satterthwaite p-value}
#'   \item{n_obs}{Number of complete observations used}
#'   \item{fdr}{FDR-corrected p-value (within model × term)}
#' }
#'
#' @examples
#' # Single gene
#' models_irf1 <- run_mediation_models(
#'   dat          = dat_Meta_w_clinical_ALI_Batch_5_RV16_Baseline_postinfection_viral_load,
#'   genes        = "IRF1",
#'   treat_col    = "Treatment_Infection_alt",
#'   mediator_col = "geneExp",
#'   outcome_col  = "log10_viral_load_24hr",
#'   covariates   = c("sex", "median_cv_coverage"),
#'   ref_level    = "No_Treatment_RV16",
#'   treat_levels = c("IL-13_Low_RV16", "IL-13_High_RV16")
#' )
#'
#' # Multiple genes
#' my_genes <- c("IRF1", "OAS1", "MX1", "IFIT1")
#' models_all <- run_mediation_models(
#'   dat          = dat_Meta_w_clinical_ALI_Batch_5_RV16_Baseline_postinfection_viral_load,
#'   genes        = my_genes,
#'   treat_col    = "Treatment_Infection_alt",
#'   mediator_col = "geneExp",
#'   outcome_col  = "log10_viral_load_24hr",
#'   covariates   = c("sex", "median_cv_coverage"),
#'   ref_level    = "No_Treatment_RV16",
#'   treat_levels = c("IL-13_Low_RV16", "IL-13_High_RV16")
#' )
#'
#' # Filter to just treatment coefficients from the mediator model
#' models_all %>%
#'   filter(model == "mediator", grepl("Treatment", term)) %>%
#'   select(gene, term, estimate, p_value, fdr)
#'
#' # Filter to gene expression coefficient in the outcome model
#' # (the "b path" — does the mediator predict the outcome?)
#' models_all %>%
#'   filter(model == "outcome", term == "geneExp") %>%
#'   select(gene, estimate, p_value, fdr) %>%
#'   arrange(p_value)
#'
#' @export
run_mediation_models <- function(
  dat,
  genes = NULL,
  treat_col = NULL,
  mediator_col = NULL,
  outcome_col = NULL,
  covariates = NULL,
  random_effect = "donorId",
  gene_col = "hgnc_symbol",
  ref_level = NULL,
  treat_levels = NULL,
  min_obs = 10
) {
  # --- Input validation -------------------------------------------------------
  if (is.null(treat_col)) {
    stop("`treat_col` must be provided (e.g., 'Treatment_Infection_alt')")
  }
  if (is.null(mediator_col)) {
    stop("`mediator_col` must be provided (e.g., 'geneExp')")
  }
  if (is.null(outcome_col)) {
    stop("`outcome_col` must be provided (e.g., 'log10_viral_load_24hr')")
  }
  if (is.null(covariates)) {
    stop("`covariates` must be provided (e.g., c('sex', 'median_cv_coverage'))")
  }
  if (is.null(ref_level)) {
    stop("`ref_level` must be provided (e.g., 'No_Treatment_RV16')")
  }
  if (is.null(treat_levels)) {
    stop(
      "`treat_levels` must be provided (e.g., c('IL-13_Low_RV16', 'IL-13_High_RV16'))"
    )
  }

  required_cols <- c(
    treat_col,
    mediator_col,
    outcome_col,
    covariates,
    random_effect,
    gene_col
  )
  missing_cols <- setdiff(required_cols, colnames(dat))
  if (length(missing_cols) > 0) {
    stop(
      "Column(s) not found in `dat`: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  observed_levels <- unique(dat[[treat_col]])
  missing_ref <- setdiff(ref_level, observed_levels)
  if (length(missing_ref) > 0) {
    stop(
      "`ref_level` '",
      ref_level,
      "' not found in column '",
      treat_col,
      "'. Available levels: ",
      paste(observed_levels, collapse = ", ")
    )
  }
  missing_treat <- setdiff(treat_levels, observed_levels)
  if (length(missing_treat) > 0) {
    stop(
      "`treat_levels` not found in column '",
      treat_col,
      "': ",
      paste(missing_treat, collapse = ", "),
      ". Available levels: ",
      paste(observed_levels, collapse = ", ")
    )
  }

  # --- Setup ------------------------------------------------------------------
  library(lmerTest)
  library(tidyverse)

  if (is.null(genes)) {
    genes <- unique(dat[[gene_col]])
    message(sprintf("No genes specified — running all %d genes", length(genes)))
  }

  covar_formula <- paste(covariates, collapse = " + ")
  re_formula <- paste0("(1 | ", random_effect, ")")

  # Mediator model: geneExp ~ treatment + covariates + (1 | donor)
  med_formula <- as.formula(
    paste(mediator_col, "~", treat_col, "+", covar_formula, "+", re_formula)
  )

  # Outcome model: viral_load ~ treatment + geneExp + covariates + (1 | donor)
  out_formula <- as.formula(
    paste(
      outcome_col,
      "~",
      treat_col,
      "+",
      mediator_col,
      "+",
      covar_formula,
      "+",
      re_formula
    )
  )

  med_formula_str <- deparse(med_formula, width.cutoff = 500)
  out_formula_str <- deparse(out_formula, width.cutoff = 500)

  results_list <- list()
  n_skipped <- 0
  n_failed <- 0

  # --- Helper to extract lmerTest summary into a tibble -----------------------
  extract_coefs <- function(fit, model_name, formula_str, gene, n_obs) {
    s <- summary(fit)
    coef_df <- as.data.frame(s$coefficients)
    tibble(
      gene = gene,
      model = model_name,
      model_formula = formula_str,
      term = rownames(coef_df),
      estimate = coef_df[["Estimate"]],
      std_error = coef_df[["Std. Error"]],
      df = coef_df[["df"]],
      t_value = coef_df[["t value"]],
      p_value = coef_df[["Pr(>|t|)"]],
      n_obs = n_obs
    )
  }

  # --- Main loop --------------------------------------------------------------
  for (i in seq_along(genes)) {
    gene <- genes[i]
    message(sprintf("[%d/%d] Fitting models for %s", i, length(genes), gene))

    dat_gene <- dat %>%
      filter(.data[[gene_col]] == gene) %>%
      mutate(
        !!treat_col := factor(
          .data[[treat_col]],
          levels = c(ref_level, treat_levels)
        )
      ) %>%
      drop_na(all_of(c(mediator_col, outcome_col, covariates, random_effect)))

    if (nrow(dat_gene) < min_obs) {
      message(sprintf(
        "  Skipping %s — only %d observations (min_obs = %d)",
        gene,
        nrow(dat_gene),
        min_obs
      ))
      n_skipped <- n_skipped + 1
      next
    }

    n <- nrow(dat_gene)

    # Mediator model
    tryCatch(
      {
        fit_m <- lmerTest::lmer(med_formula, data = dat_gene)
        results_list[[length(results_list) + 1]] <-
          extract_coefs(fit_m, "mediator", med_formula_str, gene, n)
      },
      error = function(e) {
        message(sprintf("  Mediator model failed for %s: %s", gene, e$message))
        n_failed <<- n_failed + 1
      }
    )

    # Outcome model
    tryCatch(
      {
        fit_y <- lmerTest::lmer(out_formula, data = dat_gene)
        results_list[[length(results_list) + 1]] <-
          extract_coefs(fit_y, "outcome", out_formula_str, gene, n)
      },
      error = function(e) {
        message(sprintf("  Outcome model failed for %s: %s", gene, e$message))
        n_failed <<- n_failed + 1
      }
    )
  }

  # --- Combine and FDR --------------------------------------------------------
  message(sprintf(
    "\nDone. %d genes succeeded, %d skipped (low obs), %d model failures",
    length(unique(sapply(results_list, function(x) x$gene[1]))),
    n_skipped,
    n_failed
  ))

  if (length(results_list) == 0) {
    warning("No successful model fits")
    return(tibble())
  }

  results_df <- bind_rows(results_list) %>%
    group_by(model, term) %>%
    mutate(fdr = p.adjust(p_value, method = "fdr")) %>%
    ungroup() %>%
    arrange(model, term, p_value)

  return(results_df)
}
