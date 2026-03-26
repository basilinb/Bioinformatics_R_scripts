# predictGamValues

**File:** `predictGamValues.R`
**Author:** Basilin Benson
**License:** GPL v3

## Description

Extracts fitted/predicted values from a list of GAM model objects (typically the `ModelObject` element returned by `GAMLoop`). Uses `tidymv::predict_gam` to generate predictions while excluding specified terms (e.g. donor random effects), then returns the result as a wide-format gene × condition matrix.

## Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `moduleGAM` | list | — | Output from `GAMLoop` containing a `$ModelObject` named list of fitted `gam` objects. |
| `excludeTerms` | character vector | `c("s(donorId)")` | Smooth terms to exclude from prediction (typically random-effect smooths). |
| `treat` | character | `NULL` | Name of the treatment variable column used to form wide-format column names. Must match a predictor in the model. |
| `time` | character | `NULL` | Name of the time variable column used to form wide-format column names. Must match a predictor in the model. |

## Returns

A numeric matrix:
- **Rows:** gene/module names
- **Columns:** combined `treat_time` labels (formed by `pivot_wider` over `treat` and `time`)
- **Values:** predicted GAM fit (on the log-CPM scale if input was voom-transformed)

Genes where prediction fails are silently skipped (a message is printed).

## Dependencies

`tidymv`, `dplyr`, `tidyr`, `data.table`

## Example

```r
source("GAMLoop.R")
source("predictGamValues.R")

# 1. Fit GAM models (returnModelObject must be TRUE)
gamRes <- GAMLoop(
  voomObject = v,
  gamFormula = "~ s(Timepoint, bs = 'cs') + Treatment + s(donorId, bs = 're')",
  returnModelObject = TRUE
)

# 2. Extract predicted values excluding donor random effect
predMat <- predictGamValues(
  moduleGAM    = gamRes,
  excludeTerms = c("s(donorId)"),
  treat        = "Treatment",
  time         = "Timepoint"
)

dim(predMat)   # genes × (treatment × timepoint combinations)
predMat[1:5, ]
```
