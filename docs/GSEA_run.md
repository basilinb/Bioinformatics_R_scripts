# GSEA_run

**File:** `GSEA_run.R`
**Original author:** Kim Dill-McFarland (UW)
**License:** GPL v3

## Description

Runs gene set enrichment analysis using `fgsea::fgseaSimple` across one or more ranked gene lists. Gene sets can be loaded from MSigDB via `msigdbr` (by specifying `group`/`subgroup`) or supplied as a custom named list. Automatically detects score type (standard, positive, or negative) from the range of the ranking statistic.

## Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `gene_list` | named list | — | Each element is a **named numeric vector**: gene symbols as names, ranking statistic (e.g. logFC or t-statistic) as values. One GSEA run is performed per list element. |
| `gmt_ls` | named list | `NULL` | Custom gene set list (list of character vectors of gene symbols, named by pathway). Used if `group` is not provided. |
| `nperm` | integer | `1000` | Number of permutations for p-value estimation. |
| `group` | character | `NULL` | MSigDB collection code, e.g. `"H"`, `"C2"`, `"C5"`. |
| `subgroup` | character | `NULL` | MSigDB subcollection, e.g. `"CP:REACTOME"`, `"GO:BP"`. |
| `S` | character | `"human"` | Species (currently only `"human"` is used internally). |

### Available MSigDB collections

| `group` | `subgroup` |
|---------|-----------|
| `H` | *(none)* |
| `C2` | `CGP`, `CP`, `CP:BIOCARTA`, `CP:KEGG`, `CP:PID`, `CP:REACTOME`, `CP:WIKIPATHWAYS` |
| `C5` | `GO:BP`, `GO:CC`, `GO:MF`, `HPO` |
| `C6` | *(none)* |
| `C1` | *(none)* |

## Returns

A data.frame with one row per pathway × group combination:

| Column | Description |
|--------|-------------|
| `group` | Name from `gene_list` |
| `pathway` | Pathway name (prefix and underscores stripped) |
| `fgsea.pval` | Nominal p-value |
| `fgsea.FDR` | BH-adjusted p-value |
| `fgsea.ES` | Enrichment score |
| `fgsea.NES` | Normalized enrichment score |
| `fgsea.size` | Number of genes in pathway |
| `fgsea.FC` | Direction: `"up"` or `"down"` |
| `fgsea.leadingEdge` | Semicolon-separated leading-edge genes |

## Dependencies

`tidyverse`, `fgsea`, `msigdbr`

## Example

```r
source("GSEA_run.R")
source("GSEA_plot.R")

# gene_list: named list of named numeric vectors
# e.g. gene_list <- list("TreatmentA" = setNames(logFC_vec, gene_symbols))

gsea_res <- GSEA_run(
  gene_list = gene_list,
  group     = "H"          # Hallmark gene sets
)

# Plot results
plt <- GSEA_plot(gsea_res, plot.fdr = 0.05)
print(plt)

# With a subcollection
gsea_res_reactome <- GSEA_run(
  gene_list = gene_list,
  group     = "C2",
  subgroup  = "CP:REACTOME"
)
```
