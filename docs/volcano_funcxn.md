# volcano_funcxn

**File:** `Volcano_fxn.R`

## Description

Creates volcano plots from a differential expression result data frame (e.g. output of a limma pipeline). For each coefficient/comparison, plots effect estimate on the x-axis and −log10(FDR) on the y-axis. Labels the top N up- and down-regulated significant genes by FDR, and optionally highlights additional genes of interest regardless of significance.

## Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `limma_output` | data.frame | — | Must contain columns: `volcano_comp` (comparison label), `gene` (gene ID), `hgnc_symbol`, `estimate` (effect size), `FDR`. |
| `coeff` | character vector | `NULL` | Comparison(s) to plot. Must match values in `volcano_comp`. If `NULL`, all unique comparisons are plotted. |
| `fdr` | numeric | `0.01` | FDR threshold for calling significance. |
| `n_top_gene` | integer | `10` | Number of top up- and down-regulated genes to label per plot. |
| `add_genes_int` | character vector | `NULL` | Additional gene symbols to label regardless of FDR (shown as "Gene of Interest (NS)" if not significant). |

## Returns

A named list of `ggplot2` objects, one per coefficient/comparison. Names match the values in `volcano_comp`.

Color coding:

| Color | Category |
|-------|----------|
| Red | Upregulated (significant) |
| Blue | Downregulated (significant) |
| Steel blue | Gene of Interest (not significant) |
| Grey | Not significant |

## Dependencies

`tidyverse`, `ggplot2`, `ggrepel`

## Example

```r
source("Volcano_fxn.R")

# limma_output must have: volcano_comp, gene, hgnc_symbol, estimate, FDR
plt_list <- volcano_funcxn(
  limma_output  = de_results,
  coeff         = c("TreatmentA_vs_Control"),
  fdr           = 0.05,
  n_top_gene    = 15,
  add_genes_int = c("STAT1", "IRF7", "CXCL10")
)

# Print a single plot
print(plt_list[["TreatmentA_vs_Control"]])

# Save
ggsave("volcano_TreatmentA.pdf", plt_list[["TreatmentA_vs_Control"]],
       width = 8, height = 6)
```
