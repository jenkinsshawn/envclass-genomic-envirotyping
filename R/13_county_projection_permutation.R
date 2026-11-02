if (!exists('ROOT')) source('R/00_setup.R')

## ============================================================================
## 13_county_projection_permutation.R
##
## Label-permutation null for the county random forest accuracy.
## Model: all 87 county covariates without the NASS yield-history features
## (the model behind the county map and Fig. S2), as in
## 12_county_projection_v2.R. For each maturity group the environment-type
## labels of the training counties are permuted N_PERM times, and leave-one-out
## (LOO) and leave-one-state-out (LOSO) accuracy are recomputed. The p value is
## (1 + number of permutations with accuracy >= observed) / (1 + N_PERM), with
## the observed accuracy taken as the 25-seed mean from
## out/Table_County_CV_Summary.csv (run script 12 first).
##
## Usage from the repository root:
##   Rscript R/13_county_projection_permutation.R 1000     # about 5 min on 30 cores
##   Rscript R/13_county_projection_permutation.R 20       # quick test
## Output: out/Table_County_Permutation.csv
## ============================================================================

suppressPackageStartupMessages({ library(data.table); library(randomForest); library(parallel) })

args    <- commandArgs(trailingOnly = TRUE)
N_PERM  <- if (length(args) >= 1) as.integer(args[1]) else 1000L
N_CORES <- max(1L, detectCores() - 2L)
NTREE_CV <- 1000L
YIELD_FEATURES <- c("yield_mean", "yield_trend", "yield_cv", "yield_max", "yield_recent")
t0 <- Sys.time()

## ---- training sets (same construction as script 12) ------------------------
v2 <- readRDS(file.path(DATA, "marker_v2_results.rds"))
loc_county <- fread(file.path(DATA, "county_training_locations.csv"), colClasses = list(character = "fips"))
CACHE <- file.path(DATA, "nass_cache")
power_data <- as.data.table(readRDS(file.path(CACHE, "nasa_power_counties.rds")))
wx_eng     <- as.data.table(readRDS(file.path(CACHE, "weather_engineered.rds")))
soil_data  <- as.data.table(readRDS(file.path(CACHE, "ssurgo_multidepth.rds")))
mgmt_data  <- as.data.table(readRDS(file.path(CACHE, "nass_management.rds")))
features <- copy(power_data)
wx_cols  <- setdiff(names(wx_eng), c("fips", names(power_data)))
features <- merge(features, wx_eng[, c("fips", wx_cols), with = FALSE], by = "fips", all.x = TRUE)
features <- merge(features, soil_data, by = "fips", all.x = TRUE)
features <- merge(features, mgmt_data, by.x = "fips", by.y = "fips5", all.x = TRUE)
for (mc in intersect(c("irrig_fraction", "covercrop_fraction", "tile_drain_fraction"), names(features)))
  features[is.na(get(mc)), (mc) := 0]
features[, fips := as.character(fips)]

ts <- list()
for (mg in intersect(c("0", "I", "II", "III", "IV"), names(v2))) {
  r <- v2[[mg]]$F; if (is.null(r) || r$k < 2) next
  tl <- merge(loc_county[, .(loc_key, fips, state)],
              data.table(loc_key = names(r$clusters), cluster = r$clusters), by = "loc_key")
  tl <- tl[!is.na(fips) & !is.na(cluster)][!duplicated(fips)]
  td <- merge(features, tl[, .(fips, cluster, state)], by = "fips")
  fc <- setdiff(names(td), c("fips", "cluster", "state", YIELD_FEATURES))
  fc <- fc[vapply(td[, ..fc], function(v) is.numeric(v) && sum(!is.na(v)) > nrow(td) * 0.5, logical(1))]
  for (cn in fc) { i <- is.na(td[[cn]]); if (any(i)) set(td, which(i), cn, median(td[[cn]], na.rm = TRUE)) }
  ts[[mg]] <- list(X = as.matrix(td[, ..fc]), y = factor(td$cluster), state = td$state)
}

cv_acc <- function(X, y, folds, ntree) {
  pred <- character(length(y))
  for (idx in folds) {
    rf <- randomForest(x = X[-idx, , drop = FALSE], y = droplevels(y[-idx]), ntree = ntree,
                       mtry = floor(sqrt(ncol(X))))
    pred[idx] <- as.character(predict(rf, X[idx, , drop = FALSE]))
  }
  mean(pred == as.character(y))
}

one_perm <- function(b, ts, ntree) {
  set.seed(1e5 + b)
  do.call(rbind, lapply(names(ts), function(mg) {
    X <- ts[[mg]]$X; st <- ts[[mg]]$state; yp <- sample(ts[[mg]]$y)
    data.frame(mg = mg, perm = b,
               LOO  = cv_acc(X, yp, as.list(seq_len(nrow(X))), ntree),
               LOSO = cv_acc(X, yp, lapply(unique(st), function(s) which(st == s)), ntree))
  }))
}

cl <- makeCluster(N_CORES)
clusterEvalQ(cl, library(randomForest))
clusterExport(cl, "cv_acc", envir = environment())
null <- rbindlist(parLapply(cl, seq_len(N_PERM), one_perm, ts = ts, ntree = NTREE_CV))
stopCluster(cl)

obs <- fread(file.path(OUT, "Table_County_CV_Summary.csv"))[features == "all_87", .(mg, scheme, observed = mean)]
nl  <- melt(null, id.vars = c("mg", "perm"), variable.name = "scheme", value.name = "acc")
tab <- merge(nl, obs, by = c("mg", "scheme"))[, .(
  n_perm = .N, observed = observed[1], null_mean = mean(acc), null_95th = quantile(acc, 0.95),
  p = (1 + sum(acc >= observed[1] - 1e-9)) / (1 + .N)), by = .(mg, scheme)]
fwrite(tab, file.path(OUT, "Table_County_Permutation.csv"))
print(tab, digits = 3)
cat(sprintf("Done in %.1f min (%d permutations, %d cores)\n",
            as.numeric(difftime(Sys.time(), t0, units = "mins")), N_PERM, N_CORES))
