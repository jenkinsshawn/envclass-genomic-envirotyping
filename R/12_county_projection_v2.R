if (!exists('ROOT')) source('R/00_setup.R')

## ============================================================================
## 12_county_projection_v2.R
##
## County projection of the genomic environment types, validated as the
## Methods describe it:
##   * random forest per maturity group, trained on county covariates of the
##     testing locations, WITHOUT the five NASS yield-history features
##     (87 features; the model behind the county map and Fig. S2)
##   * three feature sets: all 87 features; the top 10 by mean decrease in
##     accuracy chosen once on all training counties (selection not repeated
##     inside cross-validation, so accuracy is optimistic); and the top 10
##     chosen again inside every training fold (nested, unbiased)
##   * leave-one-out (LOO) and leave-one-state-out (LOSO) cross-validation,
##     repeated over N_SEEDS random seeds
##   * calibration: held-out class probabilities binned by confidence
##   * majority-class baseline for each maturity group
##   * county predictions from the final model and neighbor agreement within
##     100 km (county centroids, great-circle distance)
##
## Inputs (repository layout)
##   data/marker_v2_results.rds               genomic environment types
##   data/nass_cache/*.rds                    county covariates (as 06)
##   data/county_training_locations.csv       loc_key -> county fips, state
##   data/county_centroids.csv                fips, lat, lon (2,305 counties)
##   (If the two csv files are absent the script maps locations with tigris,
##    exactly as 06_county_projection.R does.)
##
## Outputs (out/)
##   Table_County_CV_Accuracy.csv        per MG x feature set x scheme x seed
##   Table_County_CV_Summary.csv         mean, sd, range over seeds + baseline
##   Table_County_Calibration.csv        confidence bins, LOO and LOSO
##   Table_County_TopFeatures.csv        top 10 features chosen once, per MG
##   Table_County_NeighborAgreement.csv  per MG map and combined MG-zone map
##   County_Predictions_v2_{all_87,top10_once}_MG{mg}.csv  fips, predicted, confidence
##   RF_Importance_v2_MG{mg}.csv         importance, 87-feature model (Fig. S2)
##
## Run time: about 2 to 4 minutes per seed on one core; seeds run in parallel.
## Usage from the repository root:
##   Rscript R/12_county_projection_v2.R            # 10 seeds
##   Rscript R/12_county_projection_v2.R 2          # quick test, 2 seeds
## ============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(randomForest)
  library(parallel)
})

args     <- commandArgs(trailingOnly = TRUE)
N_SEEDS  <- if (length(args) >= 1) as.integer(args[1]) else 25L   # 25 seeds as reported
N_CORES  <- min(N_SEEDS, max(1L, detectCores() - 2L))
TOP_K    <- 10L
NTREE_FIT <- 2000L     # importance / final model, as in 06
NTREE_CV  <- 1000L     # held-out fits, as in 06
SEED0    <- 42L
NB_KM    <- 100

YIELD_FEATURES <- c("yield_mean", "yield_trend", "yield_cv", "yield_max", "yield_recent")

cat(sprintf("10_county_projection_v2: %d seeds on %d cores\n", N_SEEDS, N_CORES))
t0 <- Sys.time()

## ---- 1. Genomic environment types and location-to-county mapping ----------
v2 <- readRDS(file.path(DATA, "marker_v2_results.rds"))

lk_file  <- file.path(DATA, "county_training_locations.csv")
cen_file <- file.path(DATA, "county_centroids.csv")

if (file.exists(lk_file) && file.exists(cen_file)) {
  loc_county <- fread(lk_file, colClasses = list(character = "fips"))
  centroids  <- fread(cen_file, colClasses = list(character = "fips"))
} else {
  ## same mapping as 06_county_projection.R
  suppressPackageStartupMessages({ library(sf); library(tigris) })
  options(tigris_use_cache = TRUE)
  nust_coords <- fread(file.path(DATA, "nust_ecw_subset.csv"),
                       select = c("Location", "Lat", "Long"))
  nust_coords[, loc_key := trimws(tolower(gsub(",.*$", "", Location)))]
  loc_coords <- nust_coords[, .(lat = median(Lat, na.rm = TRUE),
                                lon = median(Long, na.rm = TRUE)), by = loc_key]
  loc_coords <- loc_coords[!is.na(lat) & !is.na(lon)]
  us_counties <- counties(cb = TRUE, year = 2022) |> st_transform(4326)
  loc_sf <- st_as_sf(loc_coords[lat > 24 & lat < 50 & lon > -130 & lon < -60],
                     coords = c("lon", "lat"), crs = 4326)
  loc_county <- as.data.table(st_join(loc_sf, us_counties[, c("GEOID", "NAME", "STUSPS")]))
  loc_county[, geometry := NULL]
  setnames(loc_county, c("GEOID", "NAME", "STUSPS"), c("fips", "county", "state"))
  loc_county <- loc_county[!is.na(fips)]
  soy_states <- c("AL","AR","DE","GA","IA","IL","IN","KS","KY","LA","MD","MI",
                  "MN","MO","MS","NC","ND","NE","NJ","NY","OH","OK","PA","SC",
                  "SD","TN","VA","WI","WV")
  sb <- us_counties[us_counties$STUSPS %in% soy_states, ]
  xy <- st_coordinates(st_centroid(st_geometry(sb)))
  centroids <- data.table(fips = sb$GEOID, lat = xy[, 2], lon = xy[, 1])
}

## ---- 2. County feature matrix (as 06, yield history kept for reference) ----
CACHE <- file.path(DATA, "nass_cache")
power_data <- as.data.table(readRDS(file.path(CACHE, "nasa_power_counties.rds")))
wx_eng     <- as.data.table(readRDS(file.path(CACHE, "weather_engineered.rds")))
soil_data  <- as.data.table(readRDS(file.path(CACHE, "ssurgo_multidepth.rds")))
mgmt_data  <- as.data.table(readRDS(file.path(CACHE, "nass_management.rds")))
nass_yield <- as.data.table(readRDS(file.path(CACHE, "nass_soybean_yield.rds")))

features <- copy(power_data)
wx_cols  <- setdiff(names(wx_eng), c("fips", names(power_data)))
features <- merge(features, wx_eng[, c("fips", wx_cols), with = FALSE], by = "fips", all.x = TRUE)
features <- merge(features, soil_data, by = "fips", all.x = TRUE)
county_anom <- nass_yield[, {
  if (.N >= 15) {
    fit <- lm(yield ~ year)
    .(yield_mean = mean(residuals(fit)^2), yield_trend = coef(fit)[2],
      yield_cv = sd(yield) / mean(yield), yield_max = max(yield),
      yield_recent = mean(yield[year >= 2018], na.rm = TRUE))
  }
}, by = fips5]
features <- merge(features, county_anom, by.x = "fips", by.y = "fips5", all.x = TRUE)
features <- merge(features, mgmt_data, by.x = "fips", by.y = "fips5", all.x = TRUE)
for (mc in intersect(c("irrig_fraction", "covercrop_fraction", "tile_drain_fraction"), names(features)))
  features[is.na(get(mc)), (mc) := 0]
features[, fips := as.character(fips)]
cat(sprintf("Feature matrix: %d counties x %d columns\n", nrow(features), ncol(features) - 1))

## ---- 3. Training sets --------------------------------------------------------
mgs <- intersect(c("0", "I", "II", "III", "IV"), names(v2))
train_sets <- list()
for (mg in mgs) {
  r <- v2[[mg]]$F
  if (is.null(r) || r$k < 2) next
  cl  <- data.table(loc_key = names(r$clusters), cluster = r$clusters)
  tl  <- merge(loc_county[, .(loc_key, fips, state)], cl, by = "loc_key")
  tl  <- tl[!is.na(fips) & !is.na(cluster)][!duplicated(fips)]
  td  <- merge(features, tl[, .(fips, cluster, state)], by = "fips")
  fc  <- setdiff(names(td), c("fips", "cluster", "state", YIELD_FEATURES))
  fc  <- fc[vapply(td[, ..fc], function(v) is.numeric(v) && sum(!is.na(v)) > nrow(td) * 0.5, logical(1))]
  med <- vapply(fc, function(cn) median(td[[cn]], na.rm = TRUE), numeric(1))
  for (cn in fc) { i <- is.na(td[[cn]]); if (any(i)) set(td, which(i), cn, med[[cn]]) }
  train_sets[[mg]] <- list(X = as.matrix(td[, ..fc]), y = factor(td$cluster),
                           state = td$state, fips = td$fips, feat = fc, med = med)
  cat(sprintf("MG %-3s training counties %2d, states %d, types %d, features %d\n",
              mg, nrow(td), uniqueN(td$state), nlevels(factor(td$cluster)), length(fc)))
}

## ---- 4. Cross-validation ------------------------------------------------------
top_by_mda <- function(X, y, k, ntree) {
  rf <- randomForest(x = X, y = y, ntree = ntree, mtry = floor(sqrt(ncol(X))), importance = TRUE)
  imp <- importance(rf)[, "MeanDecreaseAccuracy"]
  names(sort(imp, decreasing = TRUE))[seq_len(min(k, length(imp)))]
}

cv_one_seed <- function(seed, ts, TOP_K, NTREE_FIT, NTREE_CV) {
  set.seed(seed)
  out_acc <- list(); out_prob <- list()
  for (mg in names(ts)) {
    X <- ts[[mg]]$X; y <- ts[[mg]]$y; st <- ts[[mg]]$state
    top_once <- top_by_mda(X, y, TOP_K, NTREE_FIT)
    folds <- list(LOO = as.list(seq_len(nrow(X))),
                  LOSO = lapply(sort(unique(st)), function(s) which(st == s)))
    for (scheme in names(folds)) {
      for (fs in c("all_87", "top10_once", "top10_nested")) {
        prob <- matrix(NA_real_, nrow(X), nlevels(y), dimnames = list(NULL, levels(y)))
        for (idx in folds[[scheme]]) {
          Xtr <- X[-idx, , drop = FALSE]; ytr <- droplevels(y[-idx])
          cols <- switch(fs, all_87 = colnames(X), top10_once = top_once,
                         top10_nested = top_by_mda(Xtr, ytr, TOP_K, NTREE_FIT))
          rf <- randomForest(x = Xtr[, cols, drop = FALSE], y = ytr, ntree = NTREE_CV,
                             mtry = max(1, floor(sqrt(length(cols)))))
          p <- predict(rf, X[idx, cols, drop = FALSE], type = "prob")
          prob[idx, ] <- 0
          prob[idx, colnames(p)] <- p
        }
        pred <- levels(y)[max.col(prob, ties.method = "first")]
        out_acc[[length(out_acc) + 1]] <- data.frame(
          mg = mg, features = fs, scheme = scheme, seed = seed,
          n = nrow(X), n_folds = length(folds[[scheme]]),
          accuracy = mean(pred == as.character(y)))
        out_prob[[length(out_prob) + 1]] <- data.frame(
          mg = mg, features = fs, scheme = scheme, seed = seed,
          fips = ts[[mg]]$fips, observed = as.character(y), predicted = pred,
          confidence = apply(prob, 1, max))
      }
    }
  }
  list(acc = do.call(rbind, out_acc), prob = do.call(rbind, out_prob))
}

seeds <- SEED0 + seq_len(N_SEEDS) - 1L
if (N_CORES > 1) {
  cl <- makeCluster(N_CORES)
  clusterEvalQ(cl, library(randomForest))
  clusterExport(cl, "top_by_mda", envir = environment())
  res <- parLapply(cl, seeds, cv_one_seed, ts = train_sets, TOP_K = TOP_K,
                   NTREE_FIT = NTREE_FIT, NTREE_CV = NTREE_CV)
  stopCluster(cl)
} else {
  res <- lapply(seeds, cv_one_seed, ts = train_sets, TOP_K = TOP_K,
                NTREE_FIT = NTREE_FIT, NTREE_CV = NTREE_CV)
}
acc  <- rbindlist(lapply(res, `[[`, "acc"))
prob <- rbindlist(lapply(res, `[[`, "prob"))
fwrite(prob, file.path(OUT, "Table_County_HeldOut_Predictions.csv"))
fwrite(acc, file.path(OUT, "Table_County_CV_Accuracy.csv"))

base <- rbindlist(lapply(names(train_sets), function(mg)
  data.table(mg = mg, majority_baseline = max(table(train_sets[[mg]]$y)) / length(train_sets[[mg]]$y))))
summ <- acc[, .(n = n[1], n_folds = n_folds[1], seeds = .N, mean = mean(accuracy),
                sd = sd(accuracy), min = min(accuracy), max = max(accuracy)),
            by = .(mg, features, scheme)]
summ <- merge(summ, base, by = "mg")
fwrite(summ, file.path(OUT, "Table_County_CV_Summary.csv"))

## ---- 5. Calibration ------------------------------------------------------------
prob[, bin := cut(confidence, c(0, 0.5, 0.7, 0.9, 1.0001), right = FALSE,
                  labels = c("<0.5", "0.5-0.7", "0.7-0.9", ">=0.9"))]
calib <- prob[, .(n_predictions = .N, n_per_seed = .N / uniqueN(seed),
                  mean_confidence = mean(confidence),
                  observed_accuracy = mean(predicted == observed)),
              by = .(mg, features, scheme, bin)][order(mg, features, scheme, bin)]
calib_pooled <- prob[, .(mg = "pooled", n_predictions = .N, n_per_seed = .N / uniqueN(seed),
                         mean_confidence = mean(confidence),
                         observed_accuracy = mean(predicted == observed)),
                     by = .(features, scheme, bin)]
fwrite(rbind(calib, calib_pooled, use.names = TRUE), file.path(OUT, "Table_County_Calibration.csv"))

## ---- 6. Final models, county predictions, neighbor agreement --------------------
hav <- function(la1, lo1, la2, lo2) {
  p <- pi / 180
  a <- sin((la2 - la1) * p / 2)^2 + cos(la1 * p) * cos(la2 * p) * sin((lo2 - lo1) * p / 2)^2
  2 * 6371 * asin(sqrt(a))
}
centroids <- centroids[fips %in% features$fips][order(fips)]
nC <- nrow(centroids)
Dkm <- outer(seq_len(nC), seq_len(nC), function(i, j)
  hav(centroids$lat[i], centroids$lon[i], centroids$lat[j], centroids$lon[j]))
NB <- Dkm <= NB_KM & Dkm > 0
nb_agree <- function(lab) {
  same <- outer(lab, lab, "==")
  mean(rowSums(same & NB)[rowSums(NB) > 0] / rowSums(NB)[rowSums(NB) > 0])
}

set.seed(SEED0)
pred_all <- list(); top_tab <- list(); nb_tab <- list()
for (mg in names(train_sets)) {
  ts <- train_sets[[mg]]
  for (fs in c("all_87", "top10_once")) {
    cols <- if (fs == "all_87") ts$feat else top_by_mda(ts$X, ts$y, TOP_K, NTREE_FIT)
    if (fs == "top10_once") top_tab[[mg]] <- data.table(mg = mg, rank = seq_along(cols), feature = cols)
    rf <- randomForest(x = ts$X[, cols, drop = FALSE], y = ts$y, ntree = NTREE_FIT,
                       mtry = max(1, floor(sqrt(length(cols)))), importance = TRUE)
    if (fs == "all_87") {
      im <- importance(rf)
      fwrite(data.table(feature = rownames(im), MDA = im[, "MeanDecreaseAccuracy"],
                        MDG = im[, "MeanDecreaseGini"])[order(-MDA)],
             file.path(OUT, sprintf("RF_Importance_v2_MG%s.csv", mg)))
    }
    Xp <- copy(features)[match(centroids$fips, fips)]
    for (cn in cols) { i <- is.na(Xp[[cn]]); if (any(i)) set(Xp, which(i), cn, ts$med[[cn]]) }
    pp <- predict(rf, as.matrix(Xp[, ..cols]), type = "prob")
    pr <- data.table(fips = centroids$fips, predicted = colnames(pp)[max.col(pp, "first")],
                     confidence = apply(pp, 1, max))
    pred_all[[paste(mg, fs)]] <- pr[, `:=`(mg = mg, features = fs)]
    nb_tab[[paste(mg, fs)]] <- data.table(map = paste("MG", mg), features = fs,
                                          neighbor_agreement = nb_agree(pr$predicted),
                                          share_conf_ge_0.9 = mean(pr$confidence >= 0.9))
    fwrite(pr[, .(fips, predicted, confidence)],
           file.path(OUT, sprintf("County_Predictions_v2_%s_MG%s.csv", fs, mg)))
  }
}
pred_all <- rbindlist(pred_all)
zone <- with(centroids, fcase(lat >= 45, "0", lat >= 43, "I", lat >= 40, "II", lat >= 37, "III", default = "IV"))
for (fs in c("all_87", "top10_once")) {
  pz <- pred_all[features == fs]
  lab <- vapply(seq_len(nC), function(i) {
    v <- pz[mg == zone[i] & fips == centroids$fips[i], predicted]
    if (length(v)) paste0(zone[i], "_", v) else NA_character_ }, character(1))
  nb_tab[[paste("zone", fs)]] <- data.table(map = "combined MG-zone map", features = fs,
                                            neighbor_agreement = nb_agree(lab), share_conf_ge_0.9 = NA_real_)
}
fwrite(rbindlist(top_tab), file.path(OUT, "Table_County_TopFeatures.csv"))
fwrite(rbindlist(nb_tab), file.path(OUT, "Table_County_NeighborAgreement.csv"))

## ---- 7. Console summary ---------------------------------------------------------
cat("\nAccuracy (mean over seeds) and majority-class baseline\n")
print(dcast(summ, mg + majority_baseline ~ features + scheme, value.var = "mean"), digits = 2)
cat("\nPredictions at confidence >= 0.9 (pooled over seeds)\n")
print(calib[bin == ">=0.9", .(mg, features, scheme, n_per_seed, observed_accuracy)], digits = 2)
cat("\nNeighbor agreement within 100 km\n"); print(rbindlist(nb_tab), digits = 3)
cat(sprintf("\nDone in %.1f min\n", as.numeric(difftime(Sys.time(), t0, units = "mins"))))
