# overlapGroups

**File:** `overlapGroups.R`

## Description

Helper function for UpSet-style analysis. Takes a named list of element sets and returns a grouped list where each group corresponds to a unique set-membership pattern (e.g. "only in A", "in A and B but not C", etc.). Designed to complement `UpSetR::fromList` / `ComplexUpset` workflows to allow downstream access to the actual elements in each intersection.

## Parameters

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `listInput` | named list of character vectors | — | Each element is a set of items (e.g. gene symbols). Names become the set labels. |
| `sort` | logical | `TRUE` | If `TRUE`, groups are sorted by size (largest first). |

## Returns

A named list where:
- **Names** are colon-separated set labels reflecting membership (e.g. `"setA:setB"` for elements in A and B only).
- **Values** are named integer vectors of element indices (names are the element values; indices reference the unique element list stored in `attr(result, "elements")`).
- `attr(result, "elements")` — character vector of all unique elements across all input sets.

## Dependencies

`UpSetR` (for `fromList`)

## Example

```r
library(UpSetR)
source("overlapGroups.R")

gene_sets <- list(
  TreatmentA = c("GeneA", "GeneB", "GeneC", "GeneD"),
  TreatmentB = c("GeneA", "GeneB", "GeneE"),
  TreatmentC = c("GeneA", "GeneC", "GeneF")
)

groups <- overlapGroups(gene_sets, sort = TRUE)

# Elements shared across all three sets
groups[["TreatmentA:TreatmentB:TreatmentC"]]

# Elements unique to TreatmentA
groups[["TreatmentA"]]

# All unique elements
attr(groups, "elements")
```
