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
