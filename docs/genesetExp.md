# genesetExp

**File:** `genesetExp.R`

## Description

Calculates per-library median and mean expression for genes belonging to one or more specified gene sets. Supports named MSigDB gene sets (via `msigdbr`) and/or a custom gene list. Returns a list of data frames, one per gene set, joined with the sample metadata from the voom object.

## Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `data_voom` | EList | `NULL` | `limma::voom` output. Must have `$E` (expression matrix with `lib*` column names), `$genes` (with `geneName` and `hgnc_symbol` columns), and `$targets` (sample metadata with `libid`). |
| `custom_dat` | — | `NULL` | Currently unused; reserved for future use. |
| `genesetDb_cat` | character | `"H"` | MSigDB collection category (e.g. `"H"`, `"C2"`, `"C5"`). |
| `genesetDb_sub` | character | `NULL` | MSigDB subcollection (e.g. `"GO:BP"`). |
| `genesetName` | character vector | `NULL` | One or more exact MSigDB gene set names (e.g. `"HALLMARK_INTERFERON_ALPHA_RESPONSE"`). |
| `custom_geneset` | character vector | `NULL` | Custom gene set as a vector of HGNC symbols. |
| `custom_geneset_name` | character | `"custom_geneset"` | Name to assign to the custom gene set in the output list. |

At least one of `genesetName` or `custom_geneset` must be provided.

## Returns

A named list of data frames, one per gene set. Each data frame contains:

| Column | Description |
|--------|-------------|
| `libid` | Library / sample ID |
| `mean_expression` | Mean expression across genes in the set (log-CPM scale) |
| `median_expression` | Median expression across genes in the set (log-CPM scale) |
| *...metadata columns* | All columns from `data_voom$targets` joined by `libid` |

## Dependencies

`tidyverse`, `msigdbr`

## Example

```r
source("genesetExp.R")

# MSigDB gene set
ifn_exp <- genesetExp(
  data_voom   = v,
  genesetDb_cat = "H",
  genesetName = "HALLMARK_INTERFERON_ALPHA_RESPONSE"
)

# Custom gene set
custom_exp <- genesetExp(
  data_voom          = v,
  custom_geneset     = c("STAT1", "IRF7", "MX1"),
  custom_geneset_name = "IFN_signature"
)

# Both together
combined_exp <- genesetExp(
  data_voom          = v,
  genesetName        = "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  custom_geneset     = c("STAT1", "IRF7"),
  custom_geneset_name = "key_IFN_genes"
)
```
