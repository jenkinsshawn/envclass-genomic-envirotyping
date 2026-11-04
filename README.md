# Genomic envirotyping in the Uniform Soybean Tests, Northern Region

Analysis code for Jenkins and Graef, *Genomic envirotyping: classifying
soybean production environments by marker-yield response*, Theoretical and
Applied Genetics.

## Running

Open the project at its top level. `R/00_setup.R` defines every path from the
working directory; no absolute paths appear anywhere else. To run the scripts
and check the output against the numbers in the paper:

```r
setwd("<repository root>")
source("run_all.R")
```

`run_all.R` writes `out/run_log.txt` and `out/Table_Reproducibility_Check.csv`,
which compares regenerated values (clusters, silhouettes, Table 1, the
split-half control, the permutation test, county cross-validation, yield
divergence, latitude models, the marker-count permutation test, the county
covariates by type and the mixed-model comparison) with the values
reported in the paper. Slow scripts are switched off at the top of
`run_all.R`; switch them on for the final check. Scripts that take arguments
can also be run from PowerShell, for example

```powershell
Rscript R/12_county_projection_v2.R 25
Rscript R/13_county_projection_permutation.R 1000
```

## What each script regenerates

| Script | Regenerates |
|---|---|
| `01_marker_profiles.R` | marker response profiles, clusters, subsampling consensus |
| `02_phenotypic_mantel.R` | Table 1: Mantel, partial Mantel and multiple regression on distance matrices |
| `03_disjoint_genotype_control.R` | split-half control (Fig. 2e) |
| `04_method_grid.R` | 15 combinations of distance and clustering method |
| `05_permutation_null.R` | earlier permutation test on chromosome-level profiles; kept for reference, not used in the paper |
| `06_county_projection.R` | earlier 92-feature county model including NASS yield-history variables; kept for reference, not used in the paper |
| `07_publication_figures.R` | earlier figure versions; kept for reference (the figures in the paper come from `figures_src/make_figures.py`) |
| `08_shared_lines_comparison.R` | groupings from shared lines (Fig. 3) |
| `09_fa_loadings_local.R` | factor-analytic groupings, Online Resource Table S3 (needs sommer) |
| `10_permutation_null_v2.R` | permutation test of cluster separation, Online Resource Fig. S1 |
| `11_partial_mantel_power.R` | power of the partial Mantel test |
| `12_county_projection_v2.R` | 87-feature county model without yield history: cross-validation, calibration, county predictions and variable importance (Online Resource Table S4; data for Fig. 5 and Online Resource Fig. S3) |
| `13_county_projection_permutation.R` | permutation test of county cross-validation accuracy (run after 12) |
| `14_raw_yield_response.R` | phenotypic response distance from raw yields, Online Resource 1 Supplementary Methods |
| `15_type_profile_null.R` | permutation test of the number of markers that differ between types (Table 2), and of the correlation between type mean profiles (MG I, II, IV) |
| `16_figure_data.R` | location pairs, two-dimensional scaling, location coordinates and county growing degree days and vapor pressure deficit for the figures (`out/figure_data/`) |
| `17_covariates_by_type.R` | share of the variation in each county covariate that lies between types, and correlation of growing degree days and vapor pressure deficit with latitude (Fig. 1, Results) |
| `20_yield_divergence.R` | yield divergence between types, random-grouping comparison, yield penalty |
| `21_latitude_test.R` | latitude models, latitude overlap, same-latitude pairs |
| `22_marker_density.R` | marker subsets, Online Resource Table S2 |
| `23_blup_comparison.R` | mixed-model comparison, Online Resource Table S1 |
| `24_era_stability.R` | retention between eras, Online Resource Table S6 |
| `25_table_S5_check.R` | Online Resource Table S5 and the Gene Ontology test with the 10-gene rule |

The county maps and county cross-validation in the paper come from scripts 12
and 13. Script 06 fits the earlier 92-feature model and is kept for
reference.

## Data

`data/nust_ecw_subset.csv` is the compiled phenotype and environmental record
for the Uniform Soybean Tests, Northern Region (Wartha et al. 2025), reduced to
the columns used here. Full data: SoyBase (https://www.soybase.org/tools/nust/)
and European Variation Archive study PRJEB86055. Cached intermediate objects
are in `data/` so the scripts run without regenerating upstream steps.

Added in v1.1.0:

- `data/county_training_locations.csv`: testing locations with the county
  (FIPS) and state they fall in, so scripts 12 and 13 run without downloading
  census shapefiles.
- `data/county_centroids.csv`: centroids of the 2,305 counties used for the
  NASA POWER weather queries.
- `data/loc_fips_offline.csv`: location to county FIPS table used by
  script 11 when `sf` is unavailable.
- `data/annotation/Table_S_Discriminating_Genes_a2.csv` and
  `data/annotation/Table_GO_enrichment_a2.csv`: nearest-gene annotation and
  Gene Ontology test inputs for script 25.

`out/` and `figures/` hold the outputs reported in the paper, for comparison
with a fresh run.

## Figures

The figures in the paper (Figs. 1 to 5 and Online Resource Figs. S1 to S3)
are drawn in Python from the tables in `out/` after `run_all.R` has run:

```powershell
python figures_src/make_figures.py <repository root>
```

This needs pandas, numpy, matplotlib, cartopy and pyshp. Cartopy downloads
the Natural Earth county and state outlines on first use. Output goes to
`figures/` as TIFF (600 dpi), PDF and PNG.

## Environment

R version and package versions in `sessionInfo.txt`. Scripts 20 to 25 also
use lme4, nnet, MASS and mclust; script 09 uses sommer.
