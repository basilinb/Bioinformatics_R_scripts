# make_contrast_scatterplot

**File:** `make_contrast_scatterplot.R`

## Description

Creates a scatterplot comparing effect estimates between two differential expression contrasts. Each gene is classified by its significance pattern (significant in both, in x only, in y only, or neither), and the top N genes per significance class are labeled. Reference lines at zero and along y = ±x are included to aid interpretation.

## Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `dat` | data.frame | — | Long-format DE results. Required columns: `gene`, `hgnc_symbol`, `comparison`, `estimate`, `FDR`. |
| `comp_x` | character | — | Value in `comparison` to place on the x-axis. |
| `comp_y` | character | — | Value in `comparison` to place on the y-axis. |
| `xlab` | character | `comp_x` | X-axis label. |
| `ylab` | character | `comp_y` | Y-axis label. |
| `group_x` | character | — | Human-friendly name for `comp_x` (used in significance category labels). |
| `group_y` | character | — | Human-friendly name for `comp_y` (used in significance category labels). |
| `fdr_cutoff` | numeric | `0.05` | FDR threshold for significance. |
| `top_n` | integer | `10` | Number of top genes per significance class to label (ranked by maximum absolute effect across both contrasts). |

## Returns

A named list:

| Element | Description |
|---------|-------------|
| `scatter_data` | Wide-format data frame with estimates and FDRs for both contrasts plus significance classification |
| `top_genes` | Subset of `scatter_data` used for gene labels |
| `plot` | `ggplot2` scatterplot object |

## Dependencies

`ggplot2`, `ggrepel`, `dplyr`, `tidyr`, `rlang`

> **Note:** The plot uses `theme_classic_custom`. Ensure this theme object is defined or sourced before calling this function.

## Example

```r
source("make_contrast_scatterplot.R")

# de_results: long data frame with comparison, gene, hgnc_symbol, estimate, FDR
result <- make_contrast_scatterplot(
  dat       = de_results,
  comp_x    = "TreatmentA_vs_Control",
  comp_y    = "TreatmentB_vs_Control",
  xlab      = "Treatment A effect",
  ylab      = "Treatment B effect",
  group_x   = "Treatment A",
  group_y   = "Treatment B",
  fdr_cutoff = 0.05,
  top_n      = 10
)

print(result$plot)
head(result$scatter_data)
```
