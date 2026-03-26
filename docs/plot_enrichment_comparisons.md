# plot_enrichment_comparisons

**File:** `plot_enrichment_comparisons.R`

## Description

Creates a bubble plot visualizing enrichment results across multiple pairwise comparisons simultaneously. Each row is a pathway; each column is a comparison. Bubble size encodes the k/K ratio (fraction of pathway genes in the query) and shape/color indicates direction (up vs down in Group 1). Only pathways passing the FDR cutoff are shown.

## Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `df` | data.frame | — | Enrichment results (e.g. from `BIGprofiler`). Required columns: `group`, `FDR`, `pathway`, `k/K`. |
| `comparisons` | named list | — | Each element is a length-2 character vector of group names from `df$group`. Element 1 = "Group 1" (Up), element 2 = "Group 2" (Down). The list name becomes the column label. |
| `fdr_cutoff` | numeric | `0.05` | FDR threshold; only rows with FDR ≤ this value are plotted. |

## Returns

A `ggplot2` bubble plot object.

| Visual encoding | Meaning |
|-----------------|---------|
| Triangle up (dark red) | Up in Group 1 |
| Triangle down (dark blue) | Down in Group 1 (up in Group 2) |
| Size | k/K ratio |

## Dependencies

`dplyr`, `ggplot2`, `purrr`, `rlang`

## Example

```r
source("plot_enrichment_comparisons.R")

# Define which group pairs to compare
comparisons <- list(
  "T2 High vs T2 Low\n(Non-Asthma)" = c(
    "Up in T2 High in Non-Asthma",
    "Up in T2 Low in Non-Asthma"
  ),
  "T2 High vs T2 Low\n(Asthma)" = c(
    "Up in T2 High in Asthma",
    "Up in T2 Low in Asthma"
  ),
  "Asthma vs Non-Asthma\n(T2 High)" = c(
    "Up in Asthma in T2 High",
    "Up in Non-Asthma in T2 High"
  )
)

plt <- plot_enrichment_comparisons(
  df          = enrichment_results,
  comparisons = comparisons,
  fdr_cutoff  = 0.2
)
print(plt)
```
