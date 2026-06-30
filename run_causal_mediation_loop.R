#' Run causal mediation analysis across multiple genes
#'
#' Loops through a set of genes fitting mediator and outcome mixed-effects
#' models, then runs quasi-Bayesian causal mediation analysis using the
#' \code{mediation} package. Returns a tidy data frame of results with
#' FDR-corrected p-values computed within each treatment contrast.
#'
#' @details
#' For each gene, two \code{lme4::lmer} models are fit:
#' \itemize{
#'   \item \strong{Mediator model}: \code{mediator ~ treatment + covariates + (1 | random_effect)}
#'   \item \strong{Outcome model}: \code{outcome ~ treatment + mediator + covariates + (1 | random_effect)}
#' }
#' \code{mediation::mediate()} then decomposes the total treatment effect into:
#' \itemize{
#'   \item \strong{ACME}: Average Causal Mediation Effect (indirect, through the gene)
#'   \item \strong{ADE}: Average Direct Effect (not through the gene)
#'   \item \strong{Total}: ACME + ADE
#'   \item \strong{Proportion mediated}: ACME / Total
#' }
#' Genes that fail model fitting or mediation (e.g., convergence issues) are
#' skipped with a warning message. Genes with fewer than \code{min_obs}
#' observations after dropping NAs are also skipped.
#'
#' FDR correction (Benjamini-Hochberg) is applied separately within each
#' treatment contrast.
#'
#' @param dat Data frame containing gene expression, outcome, treatment, and
#'   covariate columns. Must contain one row per gene-sample combination.
#' @param genes Character vector of gene names to analyse. If \code{NULL}
#'   (default), all unique values in \code{gene_col} are used.
#' @param treat_col String. Column name for the treatment variable.
#' @param mediator_col String. Column name for gene expression (the mediator).
#' @param outcome_col String. Column name for the outcome variable (e.g., viral load).
#' @param covariates Character vector of covariate column names to adjust for.
#' @param random_effect String. Column name for the random effect grouping
#'   variable (e.g., donor ID). Default \code{"donorId"}.
#' @param gene_col String. Column name identifying genes. Default \code{"hgnc_symbol"}.
#' @param ref_level String. Reference (control) level of the treatment factor.
#' @param treat_levels Character vector of treatment levels to contrast against
#'   \code{ref_level}. Each generates a separate mediation analysis.
#' @param sims Integer. Number of quasi-Bayesian simulations for
#'   \code{mediation::mediate()}. Default 5000.
#' @param seed Integer. Random seed for reproducibility. Default 123.
#' @param min_obs Integer. Minimum observations required per gene after
#'   dropping NAs. Default 10.
#'
#' @return A tibble with one row per gene-contrast combination, containing:
#' \describe{
#'   \item{gene}{Gene symbol}
#'   \item{contrast}{Treatment comparison (e.g., "IL-13_Low_RV16 vs No_Treatment_RV16")}
#'   \item{n_obs}{Number of complete observations used}
#'   \item{ACME_estimate, ACME_ci_lower, ACME_ci_upper, ACME_p}{Indirect effect}
#'   \item{ADE_estimate, ADE_ci_lower, ADE_ci_upper, ADE_p}{Direct effect}
#'   \item{total_estimate, total_ci_lower, total_ci_upper, total_p}{Total effect}
#'   \item{prop_mediated, prop_med_ci_lower, prop_med_ci_upper, prop_med_p}{Proportion mediated}
#'   \item{ACME_fdr, ADE_fdr, total_fdr, prop_med_fdr}{FDR-corrected p-values (within contrast)}
#' }
#'
#' @examples
#' # Single gene
#' results_irf1 <- run_causal_mediation_loop(
#'   dat         = dat_Meta_w_clinical_ALI_Batch_5_RV16_Baseline_postinfection_viral_load,
#'   genes       = "IRF1",
#'   treat_col   = "Treatment_Infection_alt",
#'   mediator_col = "geneExp",
#'   outcome_col = "log10_viral_load_24hr",
#'   covariates  = c("sex", "median_cv_coverage"),
#'   ref_level   = "No_Treatment_RV16",
#'   treat_levels = c("IL-13_Low_RV16", "IL-13_High_RV16"),
#'   sims        = 1000
#' )
#'
#' # Multiple genes
#' my_genes <- c("IRF1", "OAS1", "MX1", "IFIT1", "IFITM1")
#' results <- run_causal_mediation_loop(
#'   dat         = dat_Meta_w_clinical_ALI_Batch_5_RV16_Baseline_postinfection_viral_load,
#'   genes       = my_genes,
#'   treat_col   = "Treatment_Infection_alt",
#'   mediator_col = "geneExp",
#'   outcome_col = "log10_viral_load_24hr",
#'   covariates  = c("sex", "median_cv_coverage"),
#'   ref_level   = "No_Treatment_RV16",
#'   treat_levels = c("IL-13_Low_RV16", "IL-13_High_RV16")
#' )
#'
#' # All genes in the dataframe (will take a while)
#' results_all <- run_causal_mediation_loop(
#'   dat         = dat_Meta_w_clinical_ALI_Batch_5_RV16_Baseline_postinfection_viral_load,
#'   treat_col   = "Treatment_Infection_alt",
#'   mediator_col = "geneExp",
#'   outcome_col = "log10_viral_load_24hr",
#'   covariates  = c("sex", "median_cv_coverage"),
#'   ref_level   = "No_Treatment_RV16",
#'   treat_levels = c("IL-13_Low_RV16", "IL-13_High_RV16"),
#'   sims        = 1000
#' )
#'
#' # Top ACME hits for IL-13 Low
#' results_all %>%
#'   filter(contrast == "IL-13_Low_RV16 vs No_Treatment_RV16") %>%
#'   select(gene, ACME_estimate, ACME_p, ACME_fdr, prop_mediated) %>%
#'   head(20)
#'
#' @export
run_causal_mediation_loop <- function(
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
  sims = 1000,
  seed = 123,
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

  # Check all required columns exist in dat
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

  # Check treatment levels exist in the data
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
  library(lme4)
  library(mediation)
  library(tidyverse)

  if (is.null(genes)) {
    genes <- unique(dat[[gene_col]])
    message(sprintf("No genes specified — running all %d genes", length(genes)))
  }

  covar_formula <- paste(covariates, collapse = " + ")
  re_formula <- paste0("(1 | ", random_effect, ")")

  med_formula <- as.formula(
    paste(mediator_col, "~", treat_col, "+", covar_formula, "+", re_formula)
  )

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

  results_list <- list()
  n_skipped <- 0
  n_failed <- 0

  # --- Main loop --------------------------------------------------------------
  for (i in seq_along(genes)) {
    gene <- genes[i]
    message(sprintf("[%d/%d] Running mediation for %s", i, length(genes), gene))

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

    fits <- tryCatch(
      {
        fit_m <- lme4::lmer(med_formula, data = dat_gene)
        fit_y <- lme4::lmer(out_formula, data = dat_gene)
        list(m = fit_m, y = fit_y)
      },
      error = function(e) {
        message(sprintf("  Model fitting failed for %s: %s", gene, e$message))
        return(NULL)
      }
    )

    if (is.null(fits)) {
      n_failed <- n_failed + 1
      next
    }

    for (treat_val in treat_levels) {
      med_result <- tryCatch(
        {
          set.seed(seed)
          mediate(
            model.m = fits$m,
            model.y = fits$y,
            mediator = mediator_col,
            treat = treat_col,
            treat.value = treat_val,
            control.value = ref_level,
            sims = sims
          )
        },
        error = function(e) {
          message(sprintf(
            "  Mediation failed for %s [%s]: %s",
            gene,
            treat_val,
            e$message
          ))
          return(NULL)
        }
      )

      if (is.null(med_result)) {
        n_failed <- n_failed + 1
        next
      }

      s <- summary(med_result)

      results_list[[length(results_list) + 1]] <- tibble(
        gene = gene,
        contrast = paste0(treat_val, " vs ", ref_level),
        n_obs = nrow(dat_gene),
        ACME_estimate = s$d0,
        ACME_ci_lower = s$d0.ci[1],
        ACME_ci_upper = s$d0.ci[2],
        ACME_p = s$d0.p,
        ADE_estimate = s$z0,
        ADE_ci_lower = s$z0.ci[1],
        ADE_ci_upper = s$z0.ci[2],
        ADE_p = s$z0.p,
        total_estimate = s$tau.coef,
        total_ci_lower = s$tau.ci[1],
        total_ci_upper = s$tau.ci[2],
        total_p = s$tau.p,
        prop_mediated = s$n0,
        prop_med_ci_lower = s$n0.ci[1],
        prop_med_ci_upper = s$n0.ci[2],
        prop_med_p = s$n0.p
      )
    }
  }

  # --- Summary ----------------------------------------------------------------
  message(sprintf(
    "\nDone. %d genes succeeded, %d skipped (low obs), %d failed (model/mediation error)",
    length(unique(sapply(results_list, function(x) x$gene))),
    n_skipped,
    n_failed
  ))

  if (length(results_list) == 0) {
    warning("No successful mediation results")
    return(tibble())
  }

  results_df <- bind_rows(results_list) %>%
    group_by(contrast) %>%
    mutate(
      ACME_fdr = p.adjust(ACME_p, method = "fdr"),
      ADE_fdr = p.adjust(ADE_p, method = "fdr"),
      total_fdr = p.adjust(total_p, method = "fdr"),
      prop_med_fdr = p.adjust(prop_med_p, method = "fdr")
    ) %>%
    ungroup() %>%
    arrange(contrast, ACME_p)

  return(results_df)
}
