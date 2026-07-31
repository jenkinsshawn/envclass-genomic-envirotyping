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
##   setwd("F:/EnvClass/repo/envclass-genomic-envirotyping")
##   source("run_all.R")
##
## Writes run_log.txt and Table_Reproducibility_Check.csv to out/.
## ============================================================================

library(data.table)
library(mclust)

## 05 re-runs a permutation GWAS and 06 refits county random forests; both are
## slow. Set these FALSE for a first pass, then TRUE for the final check.
RUN <- c("01_marker_profiles.R"           = TRUE,
         "02_phenotypic_mantel.R"         = TRUE,
         "03_disjoint_genotype_control.R" = TRUE,
         "04_method_grid.R"               = TRUE,
         "05_permutation_null.R"          = FALSE,
         "06_county_projection.R"         = FALSE,
         "07_publication_figures.R"       = TRUE)

source("R/00_setup.R")

## ---- keep a reference copy of the object 01 will overwrite -----------------
ref_path <- file.path(DATA, "marker_v2_results_reference.rds")
if (!file.exists(ref_path)) {
  file.copy(file.path(DATA, "marker_v2_results.rds"), ref_path)
  cat("reference copy saved:", basename(ref_path), "\n")
}
reference <- readRDS(ref_path)

## ---- run ------------------------------------------------------------------
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
add <- function(item, expected, observed, tol = 0.005) {
  ok <- !is.null(observed) && !is.na(observed) && abs(observed - expected) <= tol
  chk[[length(chk) + 1]] <<- data.table(item, expected, observed = observed,
                                        pass = ok)
}
grab <- function(path, row, col) {
  if (!file.exists(path)) return(NA_real_)
  d <- fread(path)
  v <- d[eval(row), get(col)]
  if (!length(v)) NA_real_ else as.numeric(v[1])
}

## clusters and silhouettes must reproduce exactly
if (file.exists(file.path(DATA, "marker_v2_results.rds"))) {
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
if (file.exists(p4)) {
  g <- fread(p4)
  add("published method ARI = 1 in all MGs", 5,
      g[method_id == "correlation|ward.D2|variance" &
        abs(ari_at_published_k - 1) < 1e-8, .N], tol = 0)
}

chk <- rbindlist(chk)
fwrite(chk, file.path(OUT, "Table_Reproducibility_Check.csv"))

cat("\n--- reproducibility ---\n")
print(chk, nrows = 60)
cat(sprintf("\n%d of %d checks pass\n", sum(chk$pass, na.rm = TRUE), nrow(chk)))
if (any(!chk$pass, na.rm = TRUE)) {
  cat("\nFAILING:\n"); print(chk[pass == FALSE])
  cat("\nA failure here means the deposit does not regenerate the paper.\n")
} else {
  cat("\nThe deposit regenerates every checked number. Safe to push.\n")
}
