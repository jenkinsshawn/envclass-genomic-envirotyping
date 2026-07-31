if (!exists('ROOT')) source('R/00_setup.R')

## ============================================================================
## NASS COUNTY EXTENSION: ALL MGs + COMBINED, v2 CLUSTERS
## ============================================================================

library(data.table)
library(randomForest)
library(sf)
library(tigris)
library(ggplot2)
options(tigris_use_cache = TRUE)
set.seed(42)

pipeline_start <- Sys.time()

cat("============================================================================\n")
cat("  NASS COUNTY EXTENSION: ALL MGs + COMBINED, v2 CLUSTERS\n")
cat("  Started:", format(pipeline_start, "%Y-%m-%d %H:%M:%S"), "\n")
cat("============================================================================\n")

CACHE <- CACHE
FIG   <- FIG
TABLE <- OUT
dir.create(FIG, showWarnings = FALSE, recursive = TRUE)
dir.create(TABLE, showWarnings = FALSE, recursive = TRUE)

## ============================================================================
## PHASE 1: LOAD v2 CLUSTERS AND LOCATION DATA
## ============================================================================

cat("\n[1] Loading clusters and location data...\n")

v2 <- readRDS(file.path(DATA, "marker_v2_results.rds"))

cat("  v2 Method F clusters:\n")
for (mg_val in names(v2)) {
  r <- v2[[mg_val]]$F
  if (is.null(r)) next
  cat(sprintf("    MG %s: k=%d, sil=%.3f, n=%d\n", mg_val, r$k, r$sil, length(r$clusters)))
}

## Location coordinates
nust_coords <- fread(file.path(DATA, "nust_ecw_subset.csv"),
                     select = c("Location", "Lat", "Long"))
nust_coords[, loc_key := trimws(tolower(gsub(",.*$", "", Location)))]
loc_coords <- nust_coords[, .(lat = median(Lat, na.rm = TRUE),
                              lon = median(Long, na.rm = TRUE)), by = loc_key]
loc_coords <- loc_coords[!is.na(lat) & !is.na(lon)]

## Map locations to counties
cat("  Loading county geometries...\n")
us_counties <- counties(cb = TRUE, year = 2022) |> st_transform(4326)
us_states   <- states(cb = TRUE, year = 2022) |> st_transform(4326)

loc_sf <- st_as_sf(loc_coords[lat > 24 & lat < 50 & lon > -130 & lon < -60],
                   coords = c("lon", "lat"), crs = 4326)
loc_county <- st_join(loc_sf, us_counties[, c("GEOID", "NAME", "STUSPS")])
loc_county_dt <- as.data.table(loc_county)
loc_county_dt[, geometry := NULL]
cat(sprintf("  Locations mapped to counties: %d\n", nrow(loc_county_dt)))

## Soybean belt
soy_states <- c("AL","AR","DE","GA","IA","IL","IN","KS","KY","LA","MD","MI",
                "MN","MO","MS","NC","ND","NE","NJ","NY","OH","OK","PA","SC",
                "SD","TN","VA","WI","WV")
soy_belt <- us_counties[us_counties$STUSPS %in% soy_states, ]
soy_outlines <- us_states[us_states$STUSPS %in% soy_states, ]
cat(sprintf("  Soybean belt counties: %d\n", nrow(soy_belt)))

## ============================================================================
## PHASE 2: LOAD CACHED COUNTY FEATURES
## ============================================================================

cat("\n[2] Loading cached county features...\n")

power_data <- readRDS(file.path(CACHE, "nasa_power_counties.rds"))
cat(sprintf("  NASA POWER: %d counties\n", nrow(power_data)))

wx_eng <- readRDS(file.path(CACHE, "weather_engineered.rds"))
cat(sprintf("  Engineered weather: %d counties, %d features\n",
            nrow(wx_eng), ncol(wx_eng) - 1))

soil_data <- readRDS(file.path(CACHE, "ssurgo_multidepth.rds"))
cat(sprintf("  SSURGO soil: %d counties, %d features\n",
            nrow(soil_data), ncol(soil_data) - 1))

nass_yield <- readRDS(file.path(CACHE, "nass_soybean_yield.rds"))
cat(sprintf("  NASS yield: %d records, %d counties\n",
            nrow(nass_yield), uniqueN(nass_yield$fips5)))

mgmt_data <- readRDS(file.path(CACHE, "nass_management.rds"))
cat(sprintf("  Management: %d counties, %d features\n",
            nrow(mgmt_data), ncol(mgmt_data) - 1))

irrig_cache <- file.path(CACHE, "nass_irrig_yield.rds")
has_irrig <- file.exists(irrig_cache)
if (has_irrig) cat("  Irrigated yield gap: available\n")

## ============================================================================
## PHASE 3: BUILD FEATURE MATRIX
## ============================================================================

cat("\n[3] Building feature matrix...\n")

features <- copy(power_data)

## Engineered weather
wx_cols <- setdiff(names(wx_eng), c("fips", names(power_data)))
features <- merge(features, wx_eng[, c("fips", wx_cols), with = FALSE],
                  by = "fips", all.x = TRUE)
cat(sprintf("  + Engineered weather: %d features\n", length(wx_cols)))

## Soil
features <- merge(features, soil_data, by = "fips", all.x = TRUE)
cat(sprintf("  + Soil: %d features\n", ncol(soil_data) - 1))

## Yield history
county_anom <- nass_yield[, {
  if (.N >= 15) {
    fit <- lm(yield ~ year)
    .(yield_mean = mean(residuals(fit)^2),
      yield_trend = coef(fit)[2],
      yield_cv = sd(yield) / mean(yield),
      yield_max = max(yield),
      yield_recent = mean(yield[year >= 2018], na.rm = TRUE))
  }
}, by = fips5]
features <- merge(features, county_anom, by.x = "fips", by.y = "fips5", all.x = TRUE)
cat(sprintf("  + Yield history: %d counties\n", nrow(county_anom)))

## Management
features <- merge(features, mgmt_data, by.x = "fips", by.y = "fips5", all.x = TRUE)

## Irrigated yield gap
if (has_irrig) {
  irrig <- readRDS(irrig_cache)
  if (length(irrig$irrigated) > 0 && nrow(irrig$irrigated) > 0 &&
      length(irrig$nonirrigated) > 0 && nrow(irrig$nonirrigated) > 0) {
    ir <- irrig$irrigated[, .(yield_irrig = mean(yield_irrig, na.rm = TRUE)), by = fips5]
    ni <- irrig$nonirrigated[, .(yield_nonirrig = mean(yield_nonirrig, na.rm = TRUE)), by = fips5]
    gap <- merge(ir, ni, by = "fips5")
    gap[, yield_gap := yield_irrig - yield_nonirrig]
    gap[, yield_gap_pct := fifelse(yield_nonirrig > 0,
                                   100 * yield_gap / yield_nonirrig, NA_real_)]
    features <- merge(features, gap[, .(fips5, yield_gap, yield_gap_pct)],
                      by.x = "fips", by.y = "fips5", all.x = TRUE)
    cat(sprintf("  + Irrigation gap: %d counties\n", nrow(gap)))
  }
}

## Impute management zeros
mgmt_zero <- intersect(c("irrig_fraction", "covercrop_fraction", "tile_drain_fraction"),
                       names(features))
for (mc in mgmt_zero) features[is.na(get(mc)), (mc) := 0]

cat(sprintf("  Final: %d counties x %d features\n", nrow(features), ncol(features) - 1))

## ============================================================================
## PHASE 4: PER-MG RF TRAINING, PREDICTION, VALIDATION
## ============================================================================

cat("\n[4] Per-MG Random Forest...\n")

nass_results <- list()

for (mg_val in names(v2)) {
  r <- v2[[mg_val]]$F
  if (is.null(r) || r$k < 2) next
  
  cat(sprintf("\n  ======== MG %s (k=%d, sil=%.3f) ========\n", mg_val, r$k, r$sil))
  
  ## Map clusters to counties
  cl_dt <- data.table(loc_key = names(r$clusters), cluster = r$clusters)
  train_locs <- merge(loc_county_dt, cl_dt, by = "loc_key")
  train_locs <- train_locs[!is.na(GEOID) & !is.na(cluster)]
  train_locs <- train_locs[!duplicated(GEOID)]
  
  if (nrow(train_locs) < 5) {
    cat("    Too few training counties, skipping\n")
    next
  }
  
  cat(sprintf("    Training counties: %d\n", nrow(train_locs)))
  for (k in sort(unique(train_locs$cluster)))
    cat(sprintf("      C%d (%d): %s\n", k, sum(train_locs$cluster == k),
                paste(train_locs[cluster == k, loc_key], collapse = ", ")))
  
  ## Build training data
  train_data <- merge(features, train_locs[, .(fips = GEOID, cluster)], by = "fips")
  meta_cols <- c("fips", "cluster")
  feat_cols <- setdiff(names(train_data), meta_cols)
  feat_cols <- feat_cols[sapply(train_data[, ..feat_cols], function(x)
    is.numeric(x) && sum(!is.na(x)) > nrow(train_data) * 0.5)]
  
  ## Median imputation
  for (col in feat_cols) {
    na_idx <- is.na(train_data[[col]])
    if (any(na_idx)) train_data[na_idx, (col) := median(train_data[[col]], na.rm = TRUE)]
  }
  
  X_train <- as.matrix(train_data[, ..feat_cols])
  y_train <- as.factor(train_data$cluster)
  cat(sprintf("    Features: %d\n", length(feat_cols)))
  
  ## Train RF
  rf <- randomForest(x = X_train, y = y_train, ntree = 2000,
                     mtry = floor(sqrt(length(feat_cols))), importance = TRUE)
  oob_acc <- 100 * (1 - rf$err.rate[2000, "OOB"])
  cat(sprintf("    OOB accuracy: %.1f%%\n", oob_acc))
  cat("    Confusion:\n")
  print(rf$confusion)
  
  ## LOO
  cat("    LOO cross-validation...\n")
  loo_pred <- character(nrow(train_data))
  for (i in seq_len(nrow(train_data))) {
    rf_loo <- randomForest(x = X_train[-i, , drop = FALSE], y = y_train[-i],
                           ntree = 1000, mtry = floor(sqrt(length(feat_cols))))
    loo_pred[i] <- as.character(predict(rf_loo, X_train[i, , drop = FALSE]))
  }
  loo_acc <- mean(loo_pred == as.character(y_train))
  cat(sprintf("    LOO accuracy: %.1f%% (%d / %d)\n",
              100 * loo_acc, sum(loo_pred == as.character(y_train)), length(y_train)))
  
  ## Variable importance
  imp <- data.table(feature = rownames(importance(rf)),
                    MDA = importance(rf)[, "MeanDecreaseAccuracy"],
                    MDG = importance(rf)[, "MeanDecreaseGini"])[order(-MDA)]
  cat("    Top 15 features:\n")
  print(imp[1:min(15, nrow(imp)), .(feature, MDA = round(MDA, 3))])
  write.csv(imp, file.path(TABLE, sprintf("RF_Importance_MG%s.csv", mg_val)),
            row.names = FALSE)
  
  ## Predict all counties
  pred_data <- copy(features)
  for (col in feat_cols) {
    na_idx <- is.na(pred_data[[col]])
    if (any(na_idx)) pred_data[na_idx, (col) := median(train_data[[col]], na.rm = TRUE)]
  }
  X_pred <- as.matrix(pred_data[, ..feat_cols])
  pred_prob <- predict(rf, newdata = X_pred, type = "prob")
  pred_data[, predicted := colnames(pred_prob)[apply(pred_prob, 1, which.max)]]
  pred_data[, confidence := apply(pred_prob, 1, max)]
  
  cat(sprintf("    Predicted: %d counties, median conf = %.2f\n",
              nrow(pred_data), median(pred_data$confidence)))
  cat("    Distribution: ")
  cat(paste(sprintf("C%s=%d", names(table(pred_data$predicted)),
                    table(pred_data$predicted)), collapse = ", "), "\n")
  write.csv(pred_data[, .(fips, predicted, confidence)],
            file.path(TABLE, sprintf("County_Predictions_MG%s.csv", mg_val)),
            row.names = FALSE)
  
  ## Yield anomaly validation
  cat("    Yield anomaly validation...\n")
  train_fips <- train_locs$GEOID
  anom_dat <- nass_yield[fips5 %in% train_fips]
  yield_valid <- NULL
  
  if (nrow(anom_dat) > 0 && uniqueN(anom_dat$fips5) >= 4) {
    anom_wide <- dcast(anom_dat, year ~ fips5, value.var = "yield", fun.aggregate = mean)
    anom_mat <- as.matrix(anom_wide[, -1])
    r_mat <- cor(anom_mat, use = "pairwise.complete.obs")
    
    fips_cl <- setNames(train_locs$cluster, train_locs$GEOID)
    shared_fips <- intersect(colnames(r_mat), names(fips_cl))
    
    if (length(shared_fips) >= 4) {
      r_sub <- r_mat[shared_fips, shared_fips]
      s_sub <- outer(fips_cl[shared_fips], fips_cl[shared_fips], "==")
      same_r <- r_sub[lower.tri(r_sub) & s_sub[lower.tri(s_sub)]]
      diff_r <- r_sub[lower.tri(r_sub) & !s_sub[lower.tri(s_sub)]]
      same_r <- same_r[!is.na(same_r)]
      diff_r <- diff_r[!is.na(diff_r)]
      
      if (length(same_r) >= 2 && length(diff_r) >= 2) {
        wt <- wilcox.test(same_r, diff_r, alternative = "greater")
        cat(sprintf("      Same: r=%.3f (n=%d), Diff: r=%.3f (n=%d), p=%.4f\n",
                    mean(same_r), length(same_r), mean(diff_r), length(diff_r), wt$p.value))
        yield_valid <- list(same = mean(same_r), diff = mean(diff_r), p = wt$p.value)
      }
    }
  }
  
  ## Per-MG Map
  colors <- c("1" = "#E41A1C", "2" = "#377EB8", "3" = "#4DAF4A",
              "4" = "#984EA3", "5" = "#FF7F00")
  map_data <- merge(soy_belt, pred_data[, .(fips, predicted, confidence)],
                    by.x = "GEOID", by.y = "fips", all.x = TRUE)
  train_sf <- st_as_sf(merge(train_locs, loc_coords, by = "loc_key")[!is.na(lat)],
                       coords = c("lon", "lat"), crs = 4326)
  
  p <- ggplot() +
    geom_sf(data = map_data, aes(fill = predicted, alpha = confidence), color = NA) +
    geom_sf(data = soy_outlines, fill = NA, color = "gray40", linewidth = 0.3) +
    geom_sf(data = train_sf, aes(shape = as.factor(cluster)),
            size = 3, color = "black") +
    scale_fill_manual(values = colors[1:r$k], name = "Env Type", na.value = "gray90") +
    scale_alpha_continuous(range = c(0.3, 1), name = "Confidence", na.value = 0.1) +
    scale_shape_manual(values = c(16, 17, 15, 18, 8)[1:r$k], name = "Training") +
    labs(title = sprintf("MG %s: Genomic Environment Types", mg_val),
         subtitle = sprintf("k=%d, sil=%.3f, OOB=%.0f%%, LOO=%.0f%%",
                            r$k, r$sil, oob_acc, 100 * loo_acc)) +
    theme_minimal()
  ggsave(file.path(FIG, sprintf("Fig_NASS_MG%s.png", mg_val)), p,
         width = 16, height = 10, dpi = 300)
  
  ## Importance plot
  p_imp <- ggplot(imp[1:min(20, nrow(imp))],
                  aes(x = reorder(feature, MDA), y = MDA)) +
    geom_col(fill = "steelblue") + coord_flip() +
    labs(title = sprintf("MG %s: RF Variable Importance", mg_val),
         x = NULL, y = "Mean Decrease Accuracy") +
    theme_minimal()
  ggsave(file.path(FIG, sprintf("Fig_Importance_MG%s.png", mg_val)), p_imp,
         width = 10, height = 8, dpi = 300)
  
  cat("    Figures saved\n")
  
  nass_results[[mg_val]] <- list(
    rf = rf, imp = imp, loo_acc = loo_acc, oob_acc = oob_acc,
    pred_data = pred_data[, .(fips, predicted, confidence)],
    train_locs = train_locs, yield_valid = yield_valid,
    k = r$k, sil = r$sil
  )
}

## ============================================================================
## PHASE 5: COMBINED CROSS-MG MAP
## ============================================================================

cat("\n[5] Combined cross-MG map...\n")

centroids <- as.data.table(st_coordinates(st_centroid(st_geometry(soy_belt))))
setnames(centroids, c("lon", "lat"))
centroids$fips <- soy_belt$GEOID

centroids[, primary_mg := fcase(
  lat >= 45, "0",
  lat >= 43 & lat < 45, "I",
  lat >= 40 & lat < 43, "II",
  lat >= 37 & lat < 40, "III",
  lat < 37, "IV"
)]

cat("  Counties per MG zone:\n")
print(centroids[, .N, by = primary_mg][order(primary_mg)])

combined <- copy(centroids)
combined[, predicted := NA_character_]
combined[, confidence := NA_real_]
combined[, mg_label := NA_character_]

for (mg_val in names(nass_results)) {
  mg_pred <- nass_results[[mg_val]]$pred_data
  mg_fips <- combined[primary_mg == mg_val, fips]
  matched <- mg_pred[fips %in% mg_fips]
  if (nrow(matched) > 0) {
    combined[fips %in% matched$fips,
             `:=`(predicted = matched$predicted,
                  confidence = matched$confidence,
                  mg_label = paste0("MG", mg_val, "_C", matched$predicted))]
  }
}

n_pred <- sum(!is.na(combined$mg_label))
cat(sprintf("  Classified: %d / %d counties\n", n_pred, nrow(combined)))
cat("  Distribution:\n")
print(combined[!is.na(mg_label), .N, by = mg_label][order(mg_label)])

mg_colors <- c(
  "MG0_C1" = "#fee5d9", "MG0_C2" = "#fcae91", "MG0_C3" = "#fb6a4a", "MG0_C4" = "#cb181d",
  "MGI_C1" = "#deebf7", "MGI_C2" = "#3182bd",
  "MGII_C1" = "#e5f5e0", "MGII_C2" = "#31a354",
  "MGIII_C1" = "#fff7bc", "MGIII_C2" = "#d95f0e", "MGIII_C3" = "#7f2704",
  "MGIV_C1" = "#f2f0f7", "MGIV_C2" = "#756bb1"
)

map_combined <- merge(soy_belt, combined[, .(fips, mg_label, confidence, primary_mg)],
                      by.x = "GEOID", by.y = "fips", all.x = TRUE)

all_train <- rbindlist(lapply(names(nass_results), function(mg_val) {
  tl <- nass_results[[mg_val]]$train_locs
  tl[, mg := mg_val]
  tl[, .(loc_key, mg, cluster, GEOID)]
}))
all_train <- merge(all_train, loc_coords, by = "loc_key")
all_train <- all_train[!is.na(lat)]
all_train_sf <- st_as_sf(all_train, coords = c("lon", "lat"), crs = 4326)

p_combined <- ggplot() +
  geom_sf(data = map_combined, aes(fill = mg_label, alpha = confidence), color = NA) +
  geom_sf(data = soy_outlines, fill = NA, color = "gray40", linewidth = 0.3) +
  geom_sf(data = all_train_sf, size = 2, color = "black", shape = 16) +
  scale_fill_manual(values = mg_colors, name = "Env Type", na.value = "gray90") +
  scale_alpha_continuous(range = c(0.3, 1), name = "Confidence", na.value = 0.1) +
  labs(title = "Genomic Environment Types Across the US Soybean Belt",
       subtitle = "Within-MG classification projected by latitude-based MG zone") +
  theme_minimal()
ggsave(file.path(FIG, "Fig_NASS_Combined.png"), p_combined,
       width = 18, height = 10, dpi = 300)
cat("  Combined map saved\n")

## MG zone reference
p_zones <- ggplot() +
  geom_sf(data = map_combined, aes(fill = primary_mg), color = NA, alpha = 0.5) +
  geom_sf(data = soy_outlines, fill = NA, color = "gray40", linewidth = 0.3) +
  scale_fill_manual(values = c("0" = "#E41A1C", "I" = "#377EB8", "II" = "#4DAF4A",
                               "III" = "#FF7F00", "IV" = "#984EA3"),
                    name = "Primary MG") +
  labs(title = "Latitude-Based Maturity Group Zones") +
  theme_minimal()
ggsave(file.path(FIG, "Fig_MG_Zones.png"), p_zones, width = 16, height = 10, dpi = 300)
cat("  MG zone map saved\n")

## ============================================================================
## PHASE 6: SUMMARY
## ============================================================================

elapsed <- difftime(Sys.time(), pipeline_start, units = "mins")

cat("\n============================================================================\n")
cat("  NASS v2 SUMMARY\n")
cat("============================================================================\n")
cat(sprintf("  Runtime: %.1f min\n", as.numeric(elapsed)))
cat(sprintf("  Feature matrix: %d counties x %d features\n",
            nrow(features), ncol(features) - 1))

cat("\n  Per-MG Results:\n")
cat(sprintf("  %-6s %4s %5s %6s %6s %6s %7s %10s\n",
            "MG", "k", "Sil", "Train", "OOB", "LOO", "Conf", "YldP"))

for (mg_val in names(nass_results)) {
  r <- nass_results[[mg_val]]
  yp <- if (!is.null(r$yield_valid)) sprintf("%.4f", r$yield_valid$p) else "N/A"
  cat(sprintf("  %-6s %4d %5.3f %6d %5.1f%% %5.1f%% %6.2f %10s\n",
              mg_val, r$k, r$sil, nrow(r$train_locs),
              r$oob_acc, 100 * r$loo_acc,
              median(r$pred_data$confidence), yp))
}

cat(sprintf("\n  Combined: %d / %d counties classified\n", n_pred, nrow(combined)))

cat("\n  v0.4 baseline comparison:\n")
cat("    OOB: 25.0%%, LOO: 30.8%%\n")
cat("    Yield validation: same=0.125, diff=0.122 (p=0.33)\n")

saveRDS(nass_results, file.path(DATA, "nass_v2_results.rds"))
cat("\n  Saved: nass_v2_results.rds\n")

cat("\n============================================================================\n")
cat("  COMPLETE\n")
cat("============================================================================\n")
