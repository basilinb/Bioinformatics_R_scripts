# vizGO

**File:** `vizGo.R`

## Description

Reduces redundant Gene Ontology (GO) terms from enrichment results using semantic similarity via the `rrvgo` package. For each group of interest, significant GO terms are fetched, a semantic similarity matrix is computed, and terms are clustered and reduced to a representative subset. Useful for summarizing large C5 GO enrichment outputs into interpretable themes.

> **Note:** Requires GO IDs from MSigDB C5 (`gs_exact_source`), so input should come from `BIGprofiler` run with `collection = "C5"`.

## Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `gsea_result` | data.frame | — | Enrichment results containing columns: `group`, `FDR`, `qvalue`, and `pathway` (matching MSigDB C5 gs_name values). |
| `variable_of_interest` | character | `"all"` | Which groups to process. `"all"` runs all unique groups; otherwise provide a character vector of group names. |
| `FDR_cutoff` | numeric | `0.01` | FDR threshold for including GO terms in the reduction. |
| `rrvgo_threshold` | numeric | `0.7` | Similarity threshold for `rrvgo::reduceSimMatrix`. Higher values = more aggressive reduction (fewer representative terms). Range: 0–1. |

## Returns

A named list, one element per group (key: `"rrvgo_<group_name>"`). Each element contains:

| Sub-element | Description |
|-------------|-------------|
| `simMatrix` | Semantic similarity matrix from `rrvgo::calculateSimMatrix` |
| `reducedTerms` | data.frame of reduced/representative GO terms from `rrvgo::reduceSimMatrix` |

Groups with fewer than 2 significant terms are skipped; their names are appended to `$Failed_modules`.

## Dependencies

`tidyverse`, `rrvgo`, `msigdbr`, `org.Hs.eg.db`

## Example

```r
source("vizGo.R")

# enrich_res: output from BIGprofiler with collection = "C5"
go_reduced <- vizGO(
  gsea_result          = enrich_res,
  variable_of_interest = "all",
  FDR_cutoff           = 0.01,
  rrvgo_threshold      = 0.7
)

# Visualize reduced terms for one group
library(rrvgo)
module_res <- go_reduced[["rrvgo_module_01"]]
treemapPlot(module_res$reducedTerms)
scatterPlot(module_res$simMatrix, module_res$reducedTerms)
```
