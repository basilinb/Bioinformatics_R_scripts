# runLmerLoop

**File:** `runLmerLoop.R`

## Description

Loops over every gene or module in a long-format expression + metadata data frame and fits a linear mixed-effects model (`lmerTest::lmer`) or a standard linear model (`lm`) to each feature. Returns a combined coefficient table with BH-corrected FDR values and a list of all fitted model objects.

## Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `expressionDatawMeta` | data.frame | — | Long-format data frame containing both expression values and metadata. Each row is one observation for one gene/module. |
| `lmerformula` | character | `NULL` | Right-hand side of the model formula as a string (starting with `~`), e.g. `"~ Treatment + (1|donorId)"`. The expression column name is prepended automatically. |
| `nameCol` | character | `NULL` | Column name that identifies gene/module names (used to split the data). |
| `expCol` | character | `NULL` | Column name containing the expression values (used as the response variable). |
| `use_lmer` | logical | `TRUE` | If `TRUE`, fits `lmerTest::lmer`; if `FALSE`, fits base `lm`. |

## Returns

A named list:

| Element | Description |
|---------|-------------|
| `resultsDf` | data.frame — all coefficients across all genes/modules with `variable`, `gene`, `estimate`, `pval`, and BH `FDR` |
| `modelObjList` | named list of fitted model objects, one per gene/module |

## Dependencies

`lmerTest`, `dplyr`, `tibble`

## Example

```r
source("runLmerLoop.R")

# expressionDatawMeta should be a long data frame, e.g.:
# | geneName | expression | Treatment | donorId | ... |

res <- runLmerLoop(
  expressionDatawMeta = long_dat,
  lmerformula = "~ Treatment + TimePoint + (1|donorId)",
  nameCol     = "geneName",
  expCol      = "expression",
  use_lmer    = TRUE
)

head(res$resultsDf)
# Access a specific model
summary(res$modelObjList[["ENSG00000001"]])
```
