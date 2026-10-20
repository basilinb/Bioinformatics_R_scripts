#' Run linear mixed-effects or linear models across genes/modules
#'
#' Fits an lmer (or lm) model independently for each unique gene or module in
#' the input data, extracts coefficient tables, and applies BH FDR correction
#' within each model term.
#'
#' @param expressionDatawMeta A data frame containing gene/module expression
#'   values, metadata covariates, and a grouping column (e.g. gene name). Not
#'   needed when \code{voomObj} is supplied.
#' @param lmerformula A character string for the right-hand side of the model
#'   formula, e.g. \code{"~ viral_load + sex + (1|donorId)"}. The expression
#'   column is prepended automatically.
#' @param nameCol Character. Column name in \code{expressionDatawMeta} that
#'   identifies each gene or module to iterate over. Not needed with
#'   \code{voomObj}.
#' @param expCol Character. Column name for the expression/response variable.
#'   Not needed with \code{voomObj}.
#' @param use_lmer Logical. If \code{TRUE} (default), fits
#'   \code{lmerTest::lmer}; if \code{FALSE}, fits \code{stats::lm}.
#' @param keep_models Logical. If \code{TRUE}, returns a list containing both
#'   the results data frame and a named list of fitted model objects. If
#'   \code{FALSE} (default), returns only the results data frame. The matrix
#'   lm path (see Details) does not create model objects.
#' @param voomObj Optional limma voom \code{EList}. When supplied, expression is
#'   taken from \code{$E} (genes x samples) and metadata from \code{$targets}
#'   (one row per sample, in the same order as the columns of \code{$E}), so no
#'   long data frame is needed. With \code{use_lmer = FALSE} this uses the
#'   matrix lm path. Voom precision weights (\code{$weights}) are not used.
#' @param genes Optional character vector of genes/modules to fit. Defaults to
#'   all rownames of \code{voomObj$E}, or all unique values of \code{nameCol}.
#' @param sampleCol Optional. Column name in \code{expressionDatawMeta} that
#'   identifies samples. Supplying it with \code{use_lmer = FALSE} turns on the
#'   matrix lm path for long-format input.
#' @param n_cores Integer. Number of cores for the per-gene loop (default 2).
#'   Values above 1 fork workers with \code{pbmcapply::pbmclapply}
#'   (falling back to \code{parallel::mclapply} if pbmcapply is not
#'   installed); 1 runs serially. Results are identical either way. Forking
#'   is unavailable on Windows, where this is set to 1. Not used by matrix lm.
#' @param on_error Either \code{"stop"} (default) to stop at the first gene
#'   whose model fails, or \code{"warn"} to skip failed genes, warn, and
#'   return results for the rest.
#'
#' @return If \code{keep_models = FALSE}, a tibble with columns \code{variable},
#'   \code{estimate}, \code{Std. Error}, \code{df} (lmer only), \code{t value},
#'   \code{pval}, \code{gene}, and \code{FDR}. If \code{keep_models = TRUE}, a
#'   named list with elements \code{resultsDf} and \code{modelObjList}.
#'
#' @details
#' \strong{Per-gene loop.} By default each gene/module is fit independently.
#' The data are split by gene once up front, instead of being scanned on every
#' iteration, which is where the old version spent most of its time. Results
#' are identical to fitting each gene separately.
#'
#' \strong{Matrix lm.} With \code{use_lmer = FALSE} and either \code{voomObj}
#' or \code{sampleCol}, all genes share one design matrix, so a single QR
#' decomposition gives every gene's coefficients, standard errors and t tests
#' (the same numbers \code{lm()} gives gene by gene). Libraries with an NA in
#' any model variable are dropped, and unused factor levels are removed, as
#' \code{lm()} does. Genes with missing expression in some libraries are fit
#' with the per-gene loop instead. If the design is rank deficient, all genes
#' fall back to the loop.
#'
#' FDR correction is applied per model term (grouped by \code{variable}) using
#' \code{stats::p.adjust(method = "BH")}. Note that storing model objects
#' (\code{keep_models = TRUE}) can be very memory-intensive, as each lmer
#' object retains a copy of its input data; with \code{n_cores > 1} they are
#' also copied back from each worker.
#'
#' @examples
#' # dat_expression: long data frame, one row per gene per library, e.g.
#' #   | geneName | geneExp | libId | log10_viral_load | sex | donorId |
#' # vwts: limma voom EList (vwts$E genes x libraries, vwts$targets one row
#' #   per library, in the same order as the columns of vwts$E)
#'
#' # ── lmer ──────────────────────────────────────────────────────────────────
#' # Long data frame, default 2 cores (same call and results as the old
#' # version)
#' results <- runLmerLoop(
#'   expressionDatawMeta = dat_expression,
#'   lmerformula = "~ log10_viral_load + sex + (1|donorId)",
#'   nameCol = "geneName",
#'   expCol = "geneExp"
#' )
#'
#' # Serial, e.g. if forking misbehaves in RStudio
#' results <- runLmerLoop(
#'   expressionDatawMeta = dat_expression,
#'   lmerformula = "~ log10_viral_load + sex + (1|donorId)",
#'   nameCol = "geneName",
#'   expCol = "geneExp",
#'   n_cores = 1
#' )
#'
#' # More cores; results are identical to the serial run
#' results <- runLmerLoop(
#'   expressionDatawMeta = dat_expression,
#'   lmerformula = "~ log10_viral_load + sex + (1|donorId)",
#'   nameCol = "geneName",
#'   expCol = "geneExp",
#'   n_cores = parallel::detectCores() - 1
#' )
#'
#' # Straight from a voom object: no long data frame, nameCol or expCol needed
#' results <- runLmerLoop(
#'   voomObj = vwts,
#'   lmerformula = "~ log10_viral_load + sex + (1|donorId)",
#'   n_cores = parallel::detectCores() - 1
#' )
#'
#' # ── lm ────────────────────────────────────────────────────────────────────
#' # All genes at once (matrix lm) from a voom object - the fastest route
#' results_lm <- runLmerLoop(
#'   voomObj = vwts,
#'   lmerformula = "~ log10_viral_load + sex",
#'   use_lmer = FALSE
#' )
#'
#' # Matrix lm from the long data frame: sampleCol says which rows belong to
#' # the same library, so genes can share one design matrix
#' results_lm <- runLmerLoop(
#'   expressionDatawMeta = dat_expression,
#'   lmerformula = "~ log10_viral_load + sex",
#'   nameCol = "geneName",
#'   expCol = "geneExp",
#'   sampleCol = "libId",
#'   use_lmer = FALSE
#' )
#'
#' # lm one gene at a time (no sampleCol/voomObj); needed to keep lm models
#' results_lm <- runLmerLoop(
#'   expressionDatawMeta = dat_expression,
#'   lmerformula = "~ log10_viral_load + sex",
#'   nameCol = "geneName",
#'   expCol = "geneExp",
#'   use_lmer = FALSE,
#'   n_cores = 4
#' )
#'
#' # ── Options ───────────────────────────────────────────────────────────────
#' # Fit only some genes
#' results_isg <- runLmerLoop(
#'   voomObj = vwts,
#'   lmerformula = "~ log10_viral_load + sex + (1|donorId)",
#'   genes = c("IFNL1", "ISG15", "MX1", "OAS1")
#' )
#'
#' # Skip genes whose model fails instead of stopping; failures are listed
#' # in a warning
#' results <- runLmerLoop(
#'   voomObj = vwts,
#'   lmerformula = "~ log10_viral_load + sex + (1|donorId)",
#'   n_cores = 4,
#'   on_error = "warn"
#' )
#'
#' # Keep model objects for diagnostics
#' results_full <- runLmerLoop(
#'   expressionDatawMeta = dat_expression,
#'   lmerformula = "~ log10_viral_load + sex + (1|donorId)",
#'   nameCol = "geneName",
#'   expCol = "geneExp",
#'   keep_models = TRUE
#' )
#' results_full$resultsDf
#' results_full$modelObjList[["IFNL1"]]
#'
#' # Matrix lm makes no model objects: this warns and modelObjList is NULL
#' results_lm_full <- runLmerLoop(
#'   voomObj = vwts,
#'   lmerformula = "~ log10_viral_load + sex",
#'   use_lmer = FALSE,
#'   keep_models = TRUE
#' )
#' results_lm_full$resultsDf
#'
#' @importFrom lmerTest lmer
#' @importFrom data.table rbindlist uniqueN as.data.table
#' @importFrom dplyr rename group_by mutate arrange any_of
#' @importFrom tibble rownames_to_column as_tibble
#' @importFrom stats as.formula lm p.adjust model.matrix complete.cases pt
#' @export
runLmerLoop <- function(
  expressionDatawMeta = NULL,
  lmerformula = NULL,
  nameCol = NULL,
  expCol = NULL,
  use_lmer = TRUE,
  keep_models = FALSE,
  voomObj = NULL,
  genes = NULL,
  sampleCol = NULL,
  n_cores = 2L,
  on_error = c("stop", "warn")
) {
  require(lmerTest, quietly = TRUE)
  require(dplyr, quietly = TRUE)
  require(tibble, quietly = TRUE)
  require(data.table, quietly = TRUE)

  on_error <- match.arg(on_error)
  useVoom <- !is.null(voomObj)

  if (useVoom && !is.null(expressionDatawMeta)) {
    stop(
      "Error: Please provide either expressionDatawMeta or voomObj, not both."
    )
  }
  if (useVoom) {
    if (is.null(lmerformula)) {
      stop("Error: Please provide an lmer formula.")
    }
  } else {
    if (is.null(lmerformula) && is.null(nameCol) && is.null(expCol)) {
      stop(
        "Error: Please provide both an lmer formula and a gene/module column name."
      )
    }
    if (is.null(lmerformula)) {
      stop("Error: Please provide an lmer formula.")
    }
    if (is.null(nameCol)) {
      stop("Error: Please provide a gene/module column name.")
    }
    if (is.null(expCol)) {
      stop("Error: Please provide a gene/module expression column name.")
    }
  }

  if (.Platform$OS.type == "windows" && n_cores > 1) {
    message("Forking is not available on Windows; running on 1 core.")
    n_cores <- 1L
  }

  # ── Inputs: gene list, and a way to get one gene's data ─────────────────────
  if (useVoom) {
    exprAll <- voomObj$E
    targets <- as.data.frame(voomObj$targets)
    if (nrow(targets) != ncol(exprAll)) {
      stop("Error: voomObj$targets must have one row per column of voomObj$E.")
    }
    if (is.null(rownames(exprAll))) {
      stop("Error: voomObj$E must have gene rownames.")
    }
    varNames <- if (is.null(genes)) rownames(exprAll) else genes
    missingGenes <- setdiff(varNames, rownames(exprAll))
    if (length(missingGenes) > 0) {
      stop(
        "Error: ",
        length(missingGenes),
        " gene(s) not in voomObj$E, e.g. ",
        paste(head(missingGenes, 5), collapse = ", ")
      )
    }
    exprMat <- exprAll[varNames, , drop = FALSE]
    respName <- if (is.null(expCol)) ".exp" else expCol
    formula <- as.formula(paste0(respName, lmerformula))
    getData <- function(i) {
      dat <- targets
      dat[[respName]] <- exprMat[i, ]
      dat
    }
  } else {
    varNames <- unique(expressionDatawMeta[[nameCol]])
    if (!is.null(genes)) {
      varNames <- varNames[varNames %in% genes]
    }
    formula <- as.formula(paste0(expCol, lmerformula))
    # Row indices per gene, found once instead of scanning the whole frame
    # for every gene. Built on first use (the matrix lm path only needs it if
    # some genes fall back to the loop), and always before forking.
    rowsByGene <- NULL
    splitRows <- function() {
      if (is.null(rowsByGene)) {
        rowsByGene <<- split(
          seq_len(nrow(expressionDatawMeta)),
          expressionDatawMeta[[nameCol]]
        )
      }
    }
    getData <- function(i) {
      expressionDatawMeta[rowsByGene[[as.character(varNames[i])]], ]
    }
  }

  # ── Matrix lm: every gene at once ───────────────────────────────────────────
  useMatrix <- !use_lmer && (useVoom || !is.null(sampleCol))

  if (useMatrix) {
    if (!useVoom) {
      longInput <- .runLmerLoop_long_to_matrix(
        expressionDatawMeta,
        lmerformula,
        nameCol,
        expCol,
        sampleCol,
        varNames
      )
      exprMat <- longInput$exprMat
      targets <- longInput$targets
    }
    matrixFit <- .runLmerLoop_lm_matrix(exprMat, targets, lmerformula, varNames)

    if (is.null(matrixFit)) {
      message(
        "Design matrix is rank deficient; fitting lm gene by gene instead."
      )
    } else {
      if (keep_models) {
        warning(
          "Using matrix lm, which does not save model objects; ",
          "modelObjList will be NULL."
        )
      }
      loopDfs <- list()
      if (length(matrixFit$fallbackIdx) > 0) {
        message(
          length(matrixFit$fallbackIdx),
          " gene(s) have missing expression values; fitting those gene by gene."
        )
        if (!useVoom) {
          splitRows()
        }
        loopDfs <- .runLmerLoop_fit_genes(
          matrixFit$fallbackIdx,
          varNames,
          getData,
          formula,
          use_lmer = FALSE,
          keep_models = FALSE,
          n_cores,
          on_error
        )$resultsList
      }
      combined <- data.table::rbindlist(c(list(matrixFit$results), loopDfs))
      # Gene-major order, as the per-gene loop gives, so FDR ties sort the same
      geneOrder <- match(as.character(combined$gene), as.character(varNames))
      combined <- combined[order(geneOrder, method = "radix"), ]
      resultsDf <- .runLmerLoop_finalize(combined)

      if (keep_models) {
        return(list("resultsDf" = resultsDf, "modelObjList" = NULL))
      } else {
        return(resultsDf)
      }
    }
  } else if (!use_lmer && !useVoom && length(varNames) > 1000) {
    message(
      "Tip: supply sampleCol (or voomObj) to fit all genes at once with ",
      "matrix lm, which is much faster."
    )
  }

  # ── Per-gene loop ───────────────────────────────────────────────────────────
  if (!useVoom) {
    splitRows()
  }
  fits <- .runLmerLoop_fit_genes(
    seq_along(varNames),
    varNames,
    getData,
    formula,
    use_lmer,
    keep_models,
    n_cores,
    on_error
  )

  resultsDf <- .runLmerLoop_finalize(
    data.table::rbindlist(fits$resultsList)
  )

  if (keep_models) {
    return(list("resultsDf" = resultsDf, "modelObjList" = fits$modelObjList))
  } else {
    return(resultsDf)
  }
}


# Fit one model per gene for the gene indices in idx, serially or forked
# across n_cores. Returns the per-gene coefficient tables (and models) in idx
# order, with failed genes handled according to on_error.
.runLmerLoop_fit_genes <- function(
  idx,
  varNames,
  getData,
  formula,
  use_lmer,
  keep_models,
  n_cores,
  on_error
) {
  fitOne <- function(i) {
    tryCatch(
      {
        dat_sub <- getData(i)
        if (use_lmer) {
          model <- lmerTest::lmer(formula, data = dat_sub)
        } else {
          model <- lm(formula, data = dat_sub)
        }
        coefs <- summary(model)$coefficients %>%
          as.data.frame() %>%
          rownames_to_column("variable") %>%
          mutate(gene = varNames[i])
        list(coefs = coefs, model = if (keep_models) model)
      },
      error = function(e) list(error = conditionMessage(e))
    )
  }

  failMessage <- function(i, msg) {
    sprintf("Model failed for '%s': %s", as.character(varNames[i]), msg)
  }

  if (n_cores > 1) {
    if (requireNamespace("pbmcapply", quietly = TRUE)) {
      fits <- pbmcapply::pbmclapply(idx, fitOne, mc.cores = n_cores)
    } else {
      message("Install pbmcapply for a progress bar; using parallel::mclapply.")
      fits <- parallel::mclapply(idx, fitOne, mc.cores = n_cores)
    }
  } else {
    fits <- vector("list", length(idx))
    pb <- txtProgressBar(
      min = 0,
      max = length(idx),
      style = 3,
      width = 50,
      char = "="
    )
    for (j in seq_along(idx)) {
      fits[[j]] <- fitOne(idx[j])
      if (on_error == "stop" && !is.null(fits[[j]]$error)) {
        close(pb)
        stop(failMessage(idx[j], fits[[j]]$error), call. = FALSE)
      }
      setTxtProgressBar(pb, j)
    }
    close(pb)
  }

  # A worker that dies returns a try-error rather than our list
  errors <- vapply(
    fits,
    function(f) {
      if (inherits(f, "try-error")) {
        as.character(f)
      } else if (!is.null(f$error)) {
        f$error
      } else {
        NA_character_
      }
    },
    character(1)
  )
  failed <- !is.na(errors)

  if (any(failed)) {
    if (on_error == "stop") {
      first <- which(failed)[1]
      stop(failMessage(idx[first], errors[first]), call. = FALSE)
    }
    warning(
      sum(failed),
      " model(s) failed and were skipped:\n",
      paste(
        head(failMessage(idx[failed], errors[failed]), 10),
        collapse = "\n"
      ),
      call. = FALSE
    )
  }

  ok <- which(!failed)
  modelObjList <- NULL
  if (keep_models) {
    modelObjList <- lapply(fits[ok], `[[`, "model")
    names(modelObjList) <- as.character(varNames[idx[ok]])
  }

  list(
    resultsList = lapply(fits[ok], `[[`, "coefs"),
    modelObjList = modelObjList
  )
}


# lm for every gene at once. Every gene shares one design matrix, so a single
# QR decomposition gives every gene's coefficients, standard errors and t
# tests - exactly what lm() gives gene by gene, without the loop.
#   exprMat     - genes x samples matrix, rows in varNames order
#   targets     - one row per column of exprMat
# Returns NULL if the design is rank deficient (the caller then loops), else
# list(results = coefficient table in gene-major order, fallbackIdx = genes
# with missing expression, which need their own design and so are not fit).
.runLmerLoop_lm_matrix <- function(exprMat, targets, lmerformula, varNames) {
  modelFormula <- as.formula(lmerformula)
  modelVars <- all.vars(modelFormula)
  # lm() drops any library with an NA in a model variable (na.omit); do the same
  keepLibs <- complete.cases(targets[, modelVars, drop = FALSE])
  # and, like lm(), ignore factor levels that no longer have any rows
  modelTargets <- droplevels(targets[keepLibs, , drop = FALSE])
  designMatrix <- model.matrix(modelFormula, modelTargets)
  qrDesign <- qr(designMatrix)
  if (
    qrDesign$rank < ncol(designMatrix) ||
      any(qrDesign$pivot != seq_len(ncol(designMatrix)))
  ) {
    return(NULL)
  }

  exprSub <- exprMat[, keepLibs, drop = FALSE]
  hasNA <- rowSums(is.na(exprSub)) > 0
  fitIdx <- which(!hasNA)

  exprByLibrary <- t(exprSub[fitIdx, , drop = FALSE])
  coefMatrix <- qr.coef(qrDesign, exprByLibrary)
  residualMatrix <- qr.resid(qrDesign, exprByLibrary)
  residualDf <- nrow(designMatrix) - ncol(designMatrix)
  residualSigma <- sqrt(colSums(residualMatrix^2) / residualDf)
  unscaledSE <- sqrt(diag(chol2inv(qr.R(qrDesign))))

  # coefficients x genes; as.vector() walks coefficients within each gene,
  # the same row order as binding the per-gene lm() tables
  standardError <- outer(unscaledSE, residualSigma)
  tStatistic <- coefMatrix / standardError
  results <- data.frame(
    variable = rep(colnames(designMatrix), times = length(fitIdx)),
    Estimate = as.vector(coefMatrix),
    `Std. Error` = as.vector(standardError),
    `t value` = as.vector(tStatistic),
    `Pr(>|t|)` = as.vector(
      2 * pt(abs(tStatistic), residualDf, lower.tail = FALSE)
    ),
    check.names = FALSE,
    stringsAsFactors = FALSE
  )
  results$gene <- rep(varNames[fitIdx], each = ncol(designMatrix))

  list(results = results, fallbackIdx = which(hasNA))
}


# Long data frame -> genes x samples matrix (rows in varNames order) plus one
# metadata row per sample, for the matrix lm path.
.runLmerLoop_long_to_matrix <- function(
  expressionDatawMeta,
  lmerformula,
  nameCol,
  expCol,
  sampleCol,
  varNames
) {
  modelVars <- all.vars(as.formula(lmerformula))
  sampleIds <- expressionDatawMeta[[sampleCol]]
  samples <- unique(sampleIds)

  # Metadata must be the same on every row for a sample
  nSampleMeta <- data.table::uniqueN(
    data.table::as.data.table(
      expressionDatawMeta[, c(sampleCol, modelVars), drop = FALSE]
    )
  )
  if (nSampleMeta != length(samples)) {
    stop(
      "Error: model variables (",
      paste(modelVars, collapse = ", "),
      ") are not constant within each ",
      sampleCol,
      ", so genes do not share one design. Drop sampleCol to fit gene by gene."
    )
  }
  targets <- as.data.frame(expressionDatawMeta)[
    match(samples, sampleIds),
    modelVars,
    drop = FALSE
  ]

  rowIdx <- match(expressionDatawMeta[[nameCol]], varNames)
  colIdx <- match(sampleIds, samples)
  keep <- !is.na(rowIdx)
  cellIdx <- rowIdx[keep] + (colIdx[keep] - 1) * length(varNames)
  if (anyDuplicated(cellIdx)) {
    stop(
      "Error: more than one row for the same ",
      nameCol,
      " and ",
      sampleCol,
      "; matrix lm needs one expression value per gene per sample."
    )
  }
  # Genes absent for a sample stay NA and are fit gene by gene
  exprMat <- matrix(NA_real_, nrow = length(varNames), ncol = length(samples))
  exprMat[cellIdx] <- expressionDatawMeta[[expCol]][keep]

  list(exprMat = exprMat, targets = targets)
}


# Coefficient table -> tibble with pval/estimate names and BH FDR per term
.runLmerLoop_finalize <- function(resultsDt) {
  resultsDt %>%
    as_tibble() %>%
    rename(any_of(c("pval" = "Pr(>|t|)", "estimate" = "Estimate"))) %>%
    group_by(variable) %>%
    mutate(FDR = stats::p.adjust(pval, method = "BH")) %>%
    arrange(FDR)
}
