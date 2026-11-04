## ============================================================================
## run_all.R
##
## Place in the repository root and run from there. Executes the staged scripts
## in order, then checks that the regenerated outputs carry the numbers the
## manuscript reports.
##
## This is the last gate before pushing. Everything up to now has verified that
## the deposit is well-formed. This verifies that it works.
##
##   setwd("<repository root>")
##   source("run_all.R")
##
## Writes run_log.txt and Table_Reproducibility_Check.csv to out/.
## ============================================================================

library(data.table)
library(mclust)

## 05 and 06 are earlier versions kept for reference and are not run: 05 tests
## the older chromosome-level profiles and 06 fits the 92-feature county model
## that includes NASS yield-history variables. The manuscript uses 10 and 12-13.
## 08, 09, 10, 11 and 13 are slow (see README for run times); set them FALSE
## for a first pass, then TRUE for the final check. 09 needs sommer.
RUN <- c("01_marker_profiles.R"               = TRUE,
         "02_phenotypic_mantel.R"             = TRUE,
         "03_disjoint_genotype_control.R"     = TRUE,
         "04_method_grid.R"                   = TRUE,
         "07_publication_figures.R"           = TRUE,
         "08_shared_lines_comparison.R"       = FALSE,
         "09_fa_loadings_local.R"             = FALSE,
         "10_permutation_null_v2.R"           = TRUE,
         "11_partial_mantel_power.R"          = TRUE,
         "12_county_projection_v2.R"          = TRUE,
         "13_county_projection_permutation.R" = TRUE,
         "14_raw_yield_response.R"            = TRUE,
         "15_type_profile_null.R"             = TRUE,
         "20_yield_divergence.R"              = TRUE,
         "21_latitude_test.R"                 = TRUE,
         "22_marker_density.R"                = TRUE,
         "23_blup_comparison.R"               = TRUE,
         "24_era_stability.R"                 = TRUE,
         "25_table_S5_check.R"                = TRUE,
         "16_figure_data.R"                   = TRUE,
         "17_covariates_by_type.R"            = TRUE,
         "18_type_prediction.R"               = TRUE)

source("R/00_setup.R")

## ---- keep a reference copy of the object 01 will overwrite -----------------
ref_path <- file.path(DATA, "marker_v2_results_reference.rds")
if (!file.exists(ref_path)) {
  file.copy(file.path(DATA, "marker_v2_results.rds"), ref_path)
  cat("reference copy saved:", basename(ref_path), "\n")
}
reference <- readRDS(ref_path)

## ---- run ------------------------------------------------------------------
T_START <- Sys.time()
log <- list()
for (nm in names(RUN)) {
  f <- file.path("R", nm)
  if (!file.exists(f)) {
    log[[nm]] <- data.table(script = nm, status = "missing", secs = NA_real_, msg = "")
    cat(sprintf("%-32s MISSING\n", nm)); next
  }
  if (!RUN[[nm]]) {
    log[[nm]] <- data.table(script = nm, status = "skipped", secs = NA_real_, msg = "")
    cat(sprintf("%-32s skipped\n", nm)); next
  }
  cat(sprintf("%-32s running...\n", nm))
  t0 <- Sys.time()
  err <- tryCatch({ source(f, local = new.env()); "" },
                  error = function(e) conditionMessage(e))
  secs <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  log[[nm]] <- data.table(script = nm,
                          status = if (nzchar(err)) "ERROR" else "ok",
                          secs = round(secs, 1), msg = err)
  cat(sprintf("%-32s %s  (%.0f s)\n", nm,
              if (nzchar(err)) paste("ERROR:", substr(err, 1, 90)) else "ok", secs))
}
log <- rbindlist(log)
writeLines(capture.output(print(log)), file.path(OUT, "run_log.txt"))

## ---- reproducibility checks -----------------------------------------------

chk <- list()
## Only outputs written during this run are checked; a table left over from an
## earlier run (or shipped in out/ for reference) counts as "not run".
fresh <- function(path) file.exists(path) && file.mtime(path) >= T_START
add <- function(item, expected, observed, tol = 0.005) {
  ok <- if (is.null(observed) || is.na(observed)) NA else abs(observed - expected) <= tol
  chk[[length(chk) + 1]] <<- data.table(item, expected, observed = if (is.null(observed)) NA_real_ else observed,
                                        pass = ok)
}
grab <- function(path, row, col) {
  if (!fresh(path)) return(NA_real_)
  d <- fread(path)
  v <- d[eval(row), get(col)]
  if (!length(v)) NA_real_ else as.numeric(v[1])
}

## clusters and silhouettes must reproduce exactly
if (fresh(file.path(DATA, "marker_v2_results.rds"))) {
  new <- readRDS(file.path(DATA, "marker_v2_results.rds"))
  for (mg in c("0", "I", "II", "III", "IV")) {
    a <- reference[[mg]]$F$clusters; b <- new[[mg]]$F$clusters
    if (!is.null(a) && !is.null(b)) {
      shared <- intersect(names(a), names(b))
      add(sprintf("MG %s clusters reproduce (ARI)", mg), 1,
          adjustedRandIndex(a[shared], b[shared]), tol = 1e-8)
    }
    add(sprintf("MG %s silhouette", mg),
        c("0"=0.344,"I"=0.336,"II"=0.377,"III"=0.203,"IV"=0.395)[[mg]],
        new[[mg]]$F$sil)
  }
}

## Table 1
p1 <- file.path(OUT, "Table_PhenoMantel_Final.csv")
for (i in seq_along(c("I","II","III","IV"))) {
  mgv <- c("I","II","III","IV")[i]
  add(sprintf("Table 1 Mantel MG %s", mgv),
      c(0.640, 0.566, 0.596, 0.561)[i],
      grab(p1, quote(mg == mgv), "mantel_genomic"))
  add(sprintf("Table 1 n MG %s", mgv),
      c(15, 17, 22, 14)[i],
      grab(p1, quote(mg == mgv), "n_locs"), tol = 0)
}

## Fig. 3
p3 <- file.path(OUT, "Table_DisjointMantel.csv")
for (i in seq_along(c("I","II","III","IV"))) {
  mgv <- c("I","II","III","IV")[i]
  add(sprintf("disjoint mean MG %s", mgv), c(0.496,0.461,0.425,0.400)[i],
      grab(p3, quote(mg == mgv), "disjoint_mean"))
}

## method grid: the published method must score 1
p4 <- file.path(OUT, "Table_Method_Grid_v2.csv")
if (fresh(p4)) {
  g <- fread(p4)
  add("published method ARI = 1 in all MGs", 5,
      g[method_id == "correlation|ward.D2|variance" &
        abs(ari_at_published_k - 1) < 1e-8, .N], tol = 0)
}


## ---- checks added in v1.1.0 ------------------------------------------------
## Each check is skipped (NA) when the script that writes the table was not run.
pn <- file.path(OUT, "permutation_null", "Table_PermNull_classification.csv")
if (fresh(pn)) {
  d <- fread(pn)
  add("permutation: no permuted silhouette >= observed (all MGs)", 0,
      sum(d$n_null_ge_obs), tol = 0)
}
cv <- file.path(OUT, "Table_County_CV_Summary.csv")
if (fresh(cv)) {
  d <- fread(cv)
  for (i in 1:5) {
    mgv <- c("0","I","II","III","IV")[i]
    add(sprintf("county LOO accuracy MG %s", mgv), c(0.668, 0.789, 0.727, 0.640, 0.603)[i],
        d[mg == mgv & features == "all_87" & scheme == "LOO", mean], tol = 0.03)
  }
}
yd <- file.path(OUT, "Table_YieldDivergence.csv")
for (i in 1:4) {
  mgv <- c("I","II","III","IV")[i]
  add(sprintf("yield divergence MG %s (kg/ha)", mgv), c(303, 292, 334, 364)[i],
      grab(yd, quote(mg == mgv & variant == "div_adjdev"), "mean_kgha"), tol = 1)
}
la <- file.path(OUT, "Table_Latitude_LOO.csv")
for (i in 1:4) {
  mgv <- c("I","II","III","IV")[i]
  add(sprintf("latitude LOO accuracy MG %s", mgv), c(0.846, 0.615, 0.448, 0.684)[i],
      grab(la, quote(mg == mgv & method == "logistic"), "loo_accuracy"), tol = 0.01)
}
rw <- file.path(OUT, "Table_RawYield_Response.csv")
for (i in 1:4) {
  mgv <- c("I","II","III","IV")[i]
  add(sprintf("raw-yield response, Mantel r with genomic distance MG %s", mgv),
      c(0.168, 0.129, 0.210, 0.185)[i],
      grab(rw, quote(mg == mgv), "mantel_genomic_raw"))
}
dm <- file.path(OUT, "Table_DiscriminatingMarkers_Null.csv")
if (fresh(dm)) {
  d <- fread(dm)
  for (i in 1:4) {
    mgv <- c("I","II","III","IV")[i]
    add(sprintf("markers differing between types MG %s", mgv), c(987, 1435, 555, 899)[i],
        d[mg == mgv, observed], tol = 0)
  }
  add("marker counts: largest permuted count below observed (MG I-IV)", 4,
      d[mg != "0", sum(null_max < observed)], tol = 0)
}
ct <- file.path(OUT, "Table_Covariates_by_Type.csv")
if (fresh(ct)) {
  d <- fread(ct)
  add("growing degree days, share of variation between types MG I", 0.533,
      d[mg == "I" & covariate == "GDD_mean", share_between_types], tol = 0.005)
  add("vapor pressure deficit, smallest share between types over MGs", 0.181,
      d[covariate == "VPD_mean", min(share_between_types)], tol = 0.005)
}
tp <- file.path(OUT, "Table_TypePrediction_Summary.csv")
if (fresh(tp)) {
  d <- fread(tp)
  for (i in 1:4) {
    mgv <- c("I","II","III","IV")[i]
    add(sprintf("held-out trials, same-type correlation, same year, MG %s", mgv), c(0.504, 0.341, 0.317, 0.309)[i],
        d[mg == mgv & source == "same_year", r_type], tol = 0.01)
  }
  add("envirotyping Mantel r MG II (86 covariates)", 0.259, fread(file.path(OUT, "Table_PhenoMantel_Final.csv"))[mg == "II", mantel_enviro], tol = 0.002)
}
b7 <- file.path(OUT, "Table_S7_BLUP_Comparison.csv")
if (fresh(b7)) {
  d <- fread(b7)
  add("mixed-model silhouette reduction, mean over MGs (%)", 38,
      d[input == "mm_resid", mean(sil_reduction_pct)], tol = 0.6)
}

chk <- rbindlist(chk)
fwrite(chk, file.path(OUT, "Table_Reproducibility_Check.csv"))

cat("\n--- reproducibility ---\n")
print(chk, nrows = 60)
cat(sprintf("\n%d pass, %d fail, %d not run in this session (of %d checks)\n",
            sum(chk$pass %in% TRUE), sum(chk$pass %in% FALSE), sum(is.na(chk$pass)), nrow(chk)))
if (any(is.na(chk$pass))) { cat("\nNOT RUN (script skipped or failed):\n"); print(chk[is.na(pass)]) }
if (any(!chk$pass, na.rm = TRUE)) {
  cat("\nFAILING:\n"); print(chk[pass == FALSE])
  cat("\nA failure here means the deposit does not regenerate the paper.\n")
} else if (any(is.na(chk$pass))) {
  cat("\nNo failures, but some checks were not run. Switch on the skipped scripts in RUN and rerun before pushing.\n")
} else {
  cat("\nThe deposit regenerates every checked number. Safe to push.\n")
}
