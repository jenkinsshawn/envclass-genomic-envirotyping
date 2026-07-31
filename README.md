# Genomic envirotyping in the Northern Uniform Soybean Tests

Analysis code for Jenkins and Graef, *Genomic envirotyping: classifying
soybean production environments by the genome's marker-yield response*,
Theoretical and Applied Genetics.

## Running

Open the project at its top level and run the scripts in `R/` in order.
`00_setup.R` defines every path from the working directory; no absolute
paths appear anywhere else.

```r

To run everything and check the output against the published numbers:

```r
setwd("<repository root>")
source("run_all.R")
```

`run_all.R` writes `out/Table_Reproducibility_Check.csv`, which compares
the regenerated clusters, silhouettes, Table 1 correlations and
disjoint-genotype means against the values reported in the paper.
for (f in sort(list.files('R', full.names = TRUE))) source(f, echo = TRUE)
```

## What each script regenerates

| Script | Regenerates |
|---|---|
| `01_marker_profiles.R` | marker response profiles, clusters, subsampling consensus, temporal windows |
| `02_phenotypic_mantel.R` | Table 1, Mantel and partial Mantel |
| `03_disjoint_genotype_control.R` | Fig. 3, disjoint-genotype control |
| `04_method_grid.R` | 15-method grid, adjusted Rand index against the reported partition |
| `05_permutation_null.R` | permutation null for cluster separation |
| `06_county_projection.R` | county projection, Fig. 4 |
| `07_publication_figures.R` | Figs. 1, 2 |

## Data

`data/nust_ecw_subset` is the compiled NUST phenotype and environmental
record (Wartha et al. 2025) reduced to the columns used here. Full data:
SoyBase (https://www.soybase.org/tools/nust/) and European Variation
Archive study PRJEB86055. Cached intermediate objects are included in
`data/` so the scripts run without regenerating upstream steps.

## Environment

R version and package versions in `sessionInfo.txt`.

