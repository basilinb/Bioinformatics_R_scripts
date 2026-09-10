# `ensembl_gene_key.R` — version-pinned Ensembl gene keys

Reproducible, trackable replacement for the hand-rolled
`if (file.exists("key_pc_mouse.csv")) ... else <biomaRt download>` blocks that were
duplicated across projects.

## Why

The old pattern had four failure modes, all observed in practice:

| Problem | Consequence |
|---|---|
| `useEnsembl()` with no `version=` | Connects to whatever release is current *at run time*. Ensembl ships ~quarterly, so the annotation drifts silently and no record of it survives. |
| `file.exists()` on an untagged filename | A stale download is preferred forever. A "refresh every 3 months" comment cannot enforce itself. |
| Relative cache path | Which copy you get depends on the working directory. One project ended up with an Ensembl 115 key at its root and an Ensembl 113 key one directory down, differing by 58 mouse / 2 human genes. |
| `filter(!is.na(symbol))` | A no-op: biomaRt returns `""`, never `NA`. And `""` is not round-trip stable — `write_csv` writes `""` as an empty field but `NA` as literal `"NA"`, while `read_csv` reads both back as `NA`. The in-memory key therefore differed on a cache hit vs a fresh download. |

## Usage

```r
options(ensembl_key_dir = path.expand("~/path/to/gene_annotation"))
source("ensembl_gene_key.R")

ensembl_release <- 115L

key.pc_mouse <- get_ensembl_pc_key("mouse", release = ensembl_release)
key.pc_human <- get_ensembl_pc_key("human", release = ensembl_release)

ensembl_release_report(ensembl_release)
#>   pinned_release live_release                status checked_on
#> 1            115          116 stale by 1 release(s) 2026-09-04
```

Joining onto a counts frame — always narrow first:

```r
species <- ensembl_species_from_ids(rownames(counts))
key.pc  <- ensembl_gene_symbol_map(switch(species,
                                          mouse = key.pc_mouse,
                                          human = key.pc_human))
stopifnot(identical(names(key.pc), c("geneName", "gene_symbol")))

counts %>% as.data.frame() %>% rownames_to_column("geneName") %>% left_join(key.pc)
```

## Functions

| Function | Returns |
|---|---|
| `get_ensembl_pc_key(species, release, dir, drop_missing_symbol = FALSE, check_release = TRUE, overwrite = FALSE)` | `geneName`, `<symbol>`, `gene_biotype`, `ensembl_release`; protein-coding only, deduplicated, sorted by `geneName` |
| `get_ensembl_full_key(species, release, dir, overwrite = FALSE)` | The unfiltered BioMart key. Cached layer; you rarely call this directly |
| `ensembl_gene_symbol_map(key)` | Exactly `geneName` + `gene_symbol` |
| `ensembl_species_from_ids(gene_ids)` | `"mouse"` or `"human"` |
| `ensembl_live_release(timeout = 10, cache_hours = 24, dir, sources)` | Current release as integer, or `NA_integer_` |
| `ensembl_release_report(release)` | One-row data.frame: pinned, live, status, date |

`species` is `"mouse"` or `"human"`. `release` is an integer Ensembl release.

## Design notes

**The release is the cache key.** `get_ensembl_full_key` caches to
`<dir>/ensembl<release>_<species>_full_gene_key.rds` via `HoardR::make_path()`. Bumping
`release` therefore misses the cache and genuinely re-downloads, while the previous
release's file stays on disk so earlier results remain reproducible.
`list.files(ensembl_key_dir(), pattern = "\\.rds$")` is a full history of every release
used.

**One variable, three uses.** `release` is passed to `version =`, interpolated into the
cache filename, and is the sole source of the `ensembl_release` column. This is a
deliberate fix for a real bug where a script pinned `version = 112`, stamped
`ensembl_release = 113`, and wrote a file named `...Ensembl112.csv`.

**Cached as `.rds`, not `.csv`.** `.rds` preserves `NA` exactly, so a cache hit and a
fresh download return identical objects. The `""` → `NA_character_` normalization happens
at download time, inside the cached expression.

**Two tiers.** The expensive `getBM` download is cached; the cheap protein-coding filter
is not. Filter policy (including `drop_missing_symbol`) can therefore change without
re-downloading.

**`drop_missing_symbol` defaults to `FALSE`.** ~592 human / ~94 mouse protein-coding genes
have no assigned symbol. For `left_join` consumers dropping them is byte-for-byte neutral
(they yield `NA` either way), but it would silently remove real genes from any `%in%`-style
gene universe. Opt in explicitly if you want the stricter key.

**Absolute cache dir, deliberately not `here::here()`.** In a project with two `.Rproj`
files, `here()` resolves differently depending on which one is open — reproducing the exact
working-directory ambiguity this is meant to prevent. Set `options(ensembl_key_dir = )` or
the `ENSEMBL_KEY_DIR` environment variable to an absolute path.

**Release detection is offline-tolerant.** `ensembl_live_release()` tries three independent
sources in order — `rest.ensembl.org/info/data`, `rest.ensembl.org/info/software`, and the
`ftp.ensembl.org/pub/` directory listing — because the REST endpoints intermittently return
HTTP 500. The answer is memoised for 24h, so a typical run makes no network call. If every
source fails you get a *message*, not a warning: **absence of a warning is not proof you
are current**, which is why `ensembl_release_report()` reports `status = "unknown"`
separately from `"current"`. `biomaRt::listEnsemblArchives()` is deliberately unused — it
probes every archive host and times out.

## Updating the pinned release

1. A run warns `Pinned Ensembl release 115 is behind the live release 116.`
2. Change one number: `ensembl_release <- 116L`.
3. Re-run. A new cache file is downloaded; the old one is kept.
4. Diff before it reaches anything downstream:
   ```r
   old <- get_ensembl_pc_key("mouse", 115L); new <- get_ensembl_pc_key("mouse", 116L)
   length(setdiff(new$geneName, old$geneName))   # gained
   length(setdiff(old$geneName, new$geneName))   # lost
   ```
5. Re-process projects only if you want the new annotation. Bumping the pin does not
   retroactively change anything already delivered.

## Requirements

`HoardR`, `biomaRt`, `dplyr`, `curl`, `jsonlite`. Assumes `dplyr`/`magrittr` are attached
(the script uses `%>%`).

Release calendar for reference: 113 = Oct 2024, 114 = May 2025, 115 = Sep 2025,
116 = Jun 2026.
