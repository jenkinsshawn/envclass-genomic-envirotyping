# Changes in v1.1.2

- Fig. 1 now shows the environment types on two county backgrounds side by
  side, growing degree days and vapor pressure deficit (May to September), each
  in seven classes holding equal numbers of counties of the mapped region.
- Added R/17_covariates_by_type.R: for each county covariate of the county
  model, the share of its variation among the training counties that lies
  between types (`out/Table_Covariates_by_Type.csv`), and the correlation of
  growing degree days and vapor pressure deficit with latitude
  (`out/Table_Covariates_Latitude.csv`).
- R/16 writes `out/figure_data/county_climate.csv` (growing degree days and
  vapor pressure deficit) in place of `county_gdd.csv`.
- run_all.R runs script 17 and checks two of its values.

# Changes in v1.1.1

- Added `figures_src/make_figures.py`, which draws Figs. 1 to 5 and Online
  Resource Figs. S1 to S3 of the submitted paper from the tables in `out/`.
- Added R/16_figure_data.R, which writes the location pairs, two-dimensional
  scaling, location coordinates and county growing degree days used by the
  figures to `out/figure_data/`.
- R/15_type_profile_null.R now covers every maturity group, with the
  Kruskal-Wallis test for any number of types, and writes
  `out/Table_DiscriminatingMarkers_Null.csv` (Table 2).
- run_all.R runs scripts 15 and 16 and checks the marker counts in Table 2.
- README.md lists the figure numbers of the submitted paper.

# Changes in v1.1.0

- Added R/08 and R/09 (groupings from shared lines; factor-analytic model).
- Added R/10_permutation_null_v2.R, which permutes line yields within each
  location of each maturity group and recomputes the marker profiles. It
  replaces the permutation figure in R/07, which pooled profiles across
  maturity groups; that call in R/07 is commented out.
- Added R/11_partial_mantel_power.R (power of the partial Mantel test).
- Added R/12_county_projection_v2.R and R/13_county_projection_permutation.R,
  the 87-feature county model without yield-history variables, with
  leave-one-out and leave-one-state-out cross-validation, calibration,
  neighbor agreement and a relabeling permutation test. R/06 is kept for
  reference.
- Added R/14_raw_yield_response.R (response distance from raw yields).
- Added R/15_type_profile_null.R (permutation test of the number of markers that differ
  between types, with profiles reclustered after line yields are shuffled within
  each location).
- Added R/20 to R/25 (yield divergence, latitude, marker density,
  mixed-model comparison, era retention, Table S5).
- Added the small data files listed in README.md and the output tables of the
  checked run in out/.
- run_all.R runs the new scripts and checks the new numbers, counting a table
  as checked only if it was written during the current run; R/05 and R/06 are
  no longer run. All 42 checks passed in a fresh R 4.6.1 session (sessionInfo.txt).
