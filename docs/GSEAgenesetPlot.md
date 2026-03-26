# GSEAgenesetPlot

**File:** `GSEAgenesetPlot.R`
**Author:** Basilin Benson
**License:** GPL v3

## Description

Creates a detailed three-panel GSEA figure for a **single** gene set:
1. **Enrichment curve** — the classic `fgsea::plotEnrichment` running-sum plot.
2. **Direction heatmap** — a 1-row color tile showing high → low statistic across the ranked list, labeled with custom group names.
3. **Ranked statistic bar plot** — bar chart of per-gene ranking statistics with a median line.

Panels are assembled with `patchwork`. Supports both MSigDB collections and a custom database.

## Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `FC_list` | named numeric vector | — | Ranking statistics (e.g. logFC or t-stat) with gene symbols as names. |
| `geneset_db` | character | `"H"` | MSigDB collection code (e.g. `"H"`, `"C2"`, `"C5"`). Ignored if `custom_db` is provided. |
| `geneset_name` | character | — | Exact name of the gene set to plot (e.g. `"HALLMARK_INTERFERON_ALPHA_RESPONSE"`). |
| `custom_db` | data.frame | `NULL` | Custom gene set database with columns `gs_name` and `gene_symbol`. Overrides `geneset_db`. |
| `right_label` | character | `"Upregulated"` | Label on the right side of the heatmap bar (high statistic side). |
| `left_label` | character | `"Downregulated"` | Label on the left side of the heatmap bar (low statistic side). |
| `highlight_gene` | character | `NULL` | A single gene symbol to annotate on the enrichment curve. |
| `right_label_col` | character | `"red"` | Color for `right_label` text and heatmap high end. |
| `left_label_col` | character | `"blue"` | Color for `left_label` text and heatmap low end. |

## Returns

A `patchwork` ggplot object (three panels combined). Print or save with `ggsave`.

## Dependencies

`ggplot2`, `msigdbr`, `fgsea`, `dplyr`, `stringr`, `patchwork`

## Example

```r
source("GSEAgenesetPlot.R")

# Using MSigDB Hallmark
plt <- GSEAgenesetPlot(
  FC_list      = setNames(logFC_vec, gene_symbols),
  geneset_db   = "H",
  geneset_name = "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  right_label  = "Severe",
  left_label   = "Non-Severe",
  highlight_gene   = "CXCL10",
  right_label_col  = "darkred",
  left_label_col   = "black"
)
print(plt)

# Using a custom database
# custom_db must have columns: gs_name, gene_symbol
plt_custom <- GSEAgenesetPlot(
  FC_list      = setNames(logFC_vec, gene_symbols),
  custom_db    = my_custom_db,
  geneset_name = "MY_CUSTOM_GENESET"
)
print(plt_custom)
```
