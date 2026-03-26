# GSEA_plot

**File:** `GSEA_plot.R`
**Original author:** Kim Dill-McFarland (UW)
**License:** GPL v3

## Description

Creates a lollipop plot of GSEA results (NES on x-axis, pathways on y-axis) from the data frame returned by `GSEA_run`. Pathways are filtered by an FDR threshold. Points are colored by significance tier. An optional `plot.groups` argument restricts the display to pathways that are significant across multiple specified groups simultaneously.

## Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `gsea_result` | data.frame | — | Output from `GSEA_run`. Must contain `fgsea.FDR`, `fgsea.NES`, `pathway`, and `group` columns. |
| `plot.fdr` | numeric | `0.05` | FDR cutoff for pathway inclusion. |
| `plot.groups` | list of character vectors | `NULL` | If provided, only pathways significant in **all** specified group sets are plotted. Each element of the list is a character vector of group names; the intersection across list elements is kept. |

## Returns

A `ggplot2` object — a faceted lollipop plot (faceted by `group`).

Color scale:

| Fill | Threshold |
|------|-----------|
| Dark red | FDR < 0.01 |
| Orange-red | FDR < 0.05 |
| Light orange | FDR < 0.1 |
| Yellow | FDR < 0.2 |
| Grey | NS |

## Dependencies

`tidyverse`, `ggplot2`

## Example

```r
source("GSEA_run.R")
source("GSEA_plot.R")

gsea_res <- GSEA_run(gene_list = gene_list, group = "H")

# Plot all pathways significant at FDR < 0.05
plt <- GSEA_plot(gsea_res, plot.fdr = 0.05)
print(plt)

# Only show pathways significant in both "TreatmentA" and "TreatmentB"
plt_shared <- GSEA_plot(
  gsea_res,
  plot.fdr    = 0.05,
  plot.groups = list(c("TreatmentA"), c("TreatmentB"))
)
print(plt_shared)
```
