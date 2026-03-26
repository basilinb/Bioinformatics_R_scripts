# BIGprofiler

**File:** `BIGprofiler.R`

## Description

Over-representation analysis (ORA) wrapper around `clusterProfiler::enricher`. Runs enrichment for one or more gene lists against MSigDB collections (via `msigdbr`) or a custom TERM2GENE-style data frame. Supports SYMBOL, ENSEMBL, and ENTREZ gene IDs. Adds k/K ratios, GO IDs (for C5), and per-group query gene counts to the output.

## Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `gene_list` | named list of character vectors | `NULL` | Each element is a character vector of gene IDs for one group. |
| `gene_df` | data.frame | `NULL` | Two-column data frame: column 1 = group label, column 2 = gene ID. Alternative to `gene_list`. |
| `ID` | character | `"SYMBOL"` | Gene identifier type: `"SYMBOL"`, `"ENSEMBL"`, or `"ENTREZ"`. |
| `species` | character | `"human"` | Species: `"human"` or `"mouse"`. |
| `collection` | character | `NULL` | MSigDB collection, e.g. `"H"`, `"C2"`, `"C5"`. Use `msigdbr::msigdbr_collections()` to see all options. |
| `subcollection` | character | `NULL` | MSigDB subcollection, e.g. `"CP:REACTOME"`, `"GO:BP"`. |
| `db` | data.frame | `NULL` | Custom database — two columns: pathway name, gene ID. Overrides `collection`. |
| `minGSSize` | integer | `10` | Minimum gene set size. |
| `maxGSSize` | integer | `500` | Maximum gene set size. |

> **Note:** `category` and `subcategory` are deprecated aliases for `collection` and `subcollection`.

## Returns

A data.frame with one row per enriched pathway × group:

| Column | Description |
|--------|-------------|
| `group` | Group name from input |
| `n_query_genes` | Total genes in the query list |
| `gs_collection` / `gs_subcollection` | Database labels |
| `pathway` | Pathway name |
| `pathway_GOID` | GO term ID (C5 only) |
| `n_pathway_genes` | Total genes in pathway |
| `n_query_genes_in_pathway` | Query genes that overlap the pathway |
| `k/K` | Overlap ratio |
| `pval` | Hypergeometric p-value |
| `FDR` | BH-adjusted p-value |
| `qvalue` | q-value |
| `genes` | List of overlapping gene IDs |

## Dependencies

`clusterProfiler`, `msigdbr`, `dplyr`, `tidyr`, `tibble`

## Example

```r
source("BIGprofiler.R")

# Named list of gene symbol vectors
gene_list <- list(
  "TreatmentA_Up"   = up_genes_A,
  "TreatmentB_Up"   = up_genes_B
)

res <- BIGprofiler(
  gene_list   = gene_list,
  ID          = "SYMBOL",
  species     = "human",
  collection  = "H"        # Hallmark
)

head(res)

# With subcollection
res_go <- BIGprofiler(
  gene_list      = gene_list,
  collection     = "C5",
  subcollection  = "GO:BP"
)

# With a custom database (two-column data frame: pathway, gene_symbol)
res_custom <- BIGprofiler(
  gene_list = gene_list,
  db        = my_custom_db
)
```
