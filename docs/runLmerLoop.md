# runLmerLoop

**File:** `runLmerLoop.R`

## Description

Fits a linear mixed-effects model (`lmerTest::lmer`) or a standard linear model (`lm`) to every gene or module and returns one combined coefficient table, with BH FDR computed within each model term.

Input can be a long-format expression + metadata data frame (one row per gene per sample), or a limma voom `EList` (`$E` genes × samples, `$targets` one row per sample).

- **Per-gene loop** (lmer, or lm without `sampleCol`/`voomObj`): the data are split by gene once, and each gene is fit separately. By default genes run on 2 cores (forked with `pbmcapply`). Set `n_cores` higher for more speed, or to 1 for a serial run. Results are identical either way.
- **Matrix lm** (`use_lmer = FALSE` with `voomObj` or `sampleCol`): every gene shares one design matrix, so a single QR decomposition gives every gene's coefficients, standard errors and t tests. The numbers are the same as `lm()` gene by gene, and it runs orders of magnitude faster. The following fall back to the per-gene loop automatically:
  - genes with missing expression in some libraries;
  - rank-deficient designs.

Existing calls return the same results as before. Matrix lm is opt-in, and the per-gene loop runs on 2 cores by default.

## Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `expressionDatawMeta` | data.frame | `NULL` | Long-format data frame with expression values and metadata. Each row is one observation for one gene/module. Omit when using `voomObj`. |
| `lmerformula` | character | `NULL` | Right-hand side of the model formula as a string (starting with `~`), e.g. `"~ Treatment + (1\|donorId)"`. The expression column is prepended automatically. |
| `nameCol` | character | `NULL` | Column that identifies genes/modules. Long input only. |
| `expCol` | character | `NULL` | Column containing expression values (the response). Long input only. |
| `use_lmer` | logical | `TRUE` | `TRUE` fits `lmerTest::lmer`; `FALSE` fits `lm`. |
| `keep_models` | logical | `FALSE` | Also return the fitted model objects. Matrix lm creates none: it warns and returns `modelObjList = NULL`. |
| `voomObj` | EList | `NULL` | voom object. Expression comes from `$E` and metadata from `$targets`, in the same sample order. Voom `$weights` are not used. |
| `genes` | character | `NULL` | Subset of genes/modules to fit. Defaults to all. |
| `sampleCol` | character | `NULL` | Sample ID column in the long data frame. With `use_lmer = FALSE`, turns on matrix lm. Model variables must be constant within each sample. |
| `n_cores` | integer | `2` | Cores for the per-gene loop (1 = serial). Above 1, genes run in forked workers (`pbmcapply::pbmclapply`, or `parallel::mclapply` if pbmcapply is missing). Always 1 on Windows. |
| `on_error` | character | `"stop"` | `"stop"`: stop at the first gene whose model fails (the old behavior). `"warn"`: skip failed genes, warn, and return the rest. |

## Returns

When `keep_models = FALSE`, a tibble grouped by `variable` with columns `variable`, `estimate`, `Std. Error`, `df` (lmer only), `t value`, `pval`, `gene` and `FDR`, sorted by `FDR`.

When `keep_models = TRUE`, a named list:

| Element | Description |
|---------|-------------|
| `resultsDf` | the tibble above |
| `modelObjList` | a list of fitted model objects named by gene/module (`NULL` for matrix lm) |

## Dependencies

`lmerTest`, `dplyr`, `tibble`, `data.table`. `pbmcapply` is optional, used for a progress bar when `n_cores > 1`.

## Example

```r
source("runLmerLoop.R")

# Long data frame, e.g.:
# | geneName | expression | libId | Treatment | donorId | ... |
res <- runLmerLoop(
  expressionDatawMeta = long_dat,
  lmerformula = "~ Treatment + TimePoint + (1|donorId)",
  nameCol     = "geneName",
  expCol      = "expression",
  n_cores     = 8
)

# Straight from voom, no long data frame needed
res <- runLmerLoop(
  voomObj     = vwts,
  lmerformula = "~ Treatment + TimePoint + (1|donorId)",
  n_cores     = 8
)

# lm for every gene at once (matrix lm)
res_lm <- runLmerLoop(
  voomObj     = vwts,
  lmerformula = "~ Treatment + median_cv_coverage",
  use_lmer    = FALSE
)

# Keep model objects for diagnostics
res_full <- runLmerLoop(
  expressionDatawMeta = long_dat,
  lmerformula = "~ Treatment + (1|donorId)",
  nameCol     = "geneName",
  expCol      = "expression",
  keep_models = TRUE
)
summary(res_full$modelObjList[["ENSG00000001"]])
```

Tests: `Rscript Testing/test_runLmerLoop.R` checks default calls against the previous version, and checks matrix lm and parallel lmer against lm/lmer fit gene by gene.
