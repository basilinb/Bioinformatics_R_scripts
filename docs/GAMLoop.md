# GAMLoop

**File:** `GAMLoop.R`
**Author:** Basilin Benson
**License:** GPL v3

## Description

Loops over every gene (from a voom object) or every module (from a data frame) and fits a Generalized Additive Model (GAM) using `mgcv::gam`. Useful for longitudinal / time-course RNA-seq designs with smooth terms and random effects. Returns fixed effects, smooth effects, model summaries, fitted model objects, and a list of any genes that failed to converge.

## Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `voomObject` | EList | `NULL` | `limma::voom` output. Must have `$E` (expression matrix) and `$targets` (metadata with `libid` and `donorId` columns). Used when `moduleDatwMeta` is not provided. |
| `moduleDatwMeta` | data.frame | `NULL` | Module-level expression + metadata in a single data frame. Column names for modules must start with `"module"`. Overrides `voomObject` if provided. |
| `gamFormula` | character | `NULL` | Right-hand side of the GAM formula as a string, e.g. `"~ s(Timepoint, bs = 'cs') + Treatment + s(donorId, bs = 're')"`. |
| `returnFixedEffects` | logical | `TRUE` | Return FDR-corrected table of parametric (fixed) effects. |
| `returnSmoothEffects` | logical | `TRUE` | Return FDR-corrected table of smooth terms. |
| `returnGeneModelSummary` | logical | `TRUE` | Return a named list of `summary(gam)` objects per gene/module. |
| `returnModelObject` | logical | `TRUE` | Return a named list of fitted `gam` objects per gene/module. |

## Returns

A named list:

| Element | Description |
|---------|-------------|
| `FixedEffects` | data.frame — parametric terms with `variable`, `geneName`, p-value columns, and BH `FDR` |
| `SmoothEffects` | data.frame — smooth terms with `variable`, `geneName`, p-value columns, and BH `FDR` |
| `GeneModelSummary` | named list of `summary.gam` objects |
| `ModelObject` | named list of fitted `gam` objects |
| `FailedGeneList` | character vector of genes/modules where fitting failed |
| `ExecutionTime` | `difftime` object — total wall-clock time |

## Dependencies

`tidyverse`, `mgcv`, `data.table`

## Example

```r
source("GAMLoop.R")

res <- GAMLoop(
  voomObject = v,                                          # limma voom object
  gamFormula  = "~ s(Timepoint, bs = 'cs') + Treatment + s(donorId, bs = 're')",
  returnFixedEffects    = TRUE,
  returnSmoothEffects   = TRUE,
  returnGeneModelSummary = FALSE,
  returnModelObject      = FALSE
)

head(res$FixedEffects)
head(res$SmoothEffects)
res$FailedGeneList
```

> **Note:** `donorId` in `voomObject$targets` is automatically coerced to a factor. Ensure your metadata contains a `donorId` column if you include it in the formula.
