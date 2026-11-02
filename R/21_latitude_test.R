if (!exists('ROOT')) source('R/00_setup.R')

## ============================================================================
## 21_latitude_test.R
## Does latitude alone predict genomic environment type at the training
## locations? (Methods and Results, "Latitude and environment type")
##
## Reconstructed: the original script was not recovered. Coordinates are
## built exactly as in R/01_marker_profiles.R (median Lat and Long of each
## location key in data/nust_ecw_subset.csv); types are the archived Method F
## assignments in data/marker_v2_results.rds.
##
## Leave-one-out accuracy for each maturity group under:
##   logistic    : binomial logistic regression on latitude (k = 2), or
##                 multinomial logistic regression (nnet::multinom) for k > 2
##   lda         : linear discriminant analysis on latitude (MASS::lda)
##   qda         : quadratic discriminant analysis (MASS::qda); each type has
##                 its own latitude variance (sensitivity)
##   nearest_mean: assign to the type whose mean latitude (training set) is
##                 nearest
##   threshold   : best cut points on latitude chosen on the training set
##                 (k - 1 ordered cut points, types ordered by mean latitude;
##                 exhaustive search over midpoints between observed values)
##   majority    : every location assigned to the largest type of the
##                 training set (the leave-one-out version of the baseline)
## plus the full-sample majority baseline (largest type / n).
##
## Latitude overlap (two types): length of the latitude interval shared by
## the types divided by the total latitude range. For k > 2 the mean over
## pairs of types is reported.
## Same-latitude pairs: location pairs within 0.5 degrees latitude assigned
## to different types (|dlat| <= 0.5), with great-circle distance.
##
## Output: out/Table_Latitude_LOO.csv, out/Table_Latitude_Overlap.csv,
##         out/Table_Latitude_Pairs.csv
## ============================================================================

library(data.table)
library(MASS)
library(nnet)
set.seed(42)

nust_coords <- fread(file.path(DATA, 'nust_ecw_subset.csv'),
                     select = c('Location', 'Lat', 'Long'))
nust_coords[, loc_key := trimws(tolower(gsub(',.*$', '', Location)))]
loc_coords <- nust_coords[, .(lat = median(Lat, na.rm = TRUE),
                              lon = median(Long, na.rm = TRUE)), by = loc_key]
mv2 <- readRDS(file.path(DATA, 'marker_v2_results.rds'))

haversine_km <- function(lat1, lon1, lat2, lon2) {
  r <- pi / 180
  a <- sin((lat2 - lat1) * r / 2)^2 +
       cos(lat1 * r) * cos(lat2 * r) * sin((lon2 - lon1) * r / 2)^2
  2 * 6371 * asin(pmin(1, sqrt(a)))
}

fit_predict <- function(method, tr, te) {
  y <- factor(tr$type); lv <- levels(y)
  if (length(lv) == 1) return(lv)
  switch(method,
    logistic = {
      if (length(lv) == 2) {
        m <- suppressWarnings(glm(y ~ lat, data = tr, family = binomial))
        p <- predict(m, newdata = te, type = 'response')
        ifelse(p > 0.5, lv[2], lv[1])
      } else {
        m <- multinom(y ~ lat, data = tr, trace = FALSE, maxit = 1000)
        as.character(predict(m, newdata = te, type = 'class'))
      }
    },
    lda = as.character(predict(lda(y ~ lat, data = tr), newdata = te)$class),
    qda = as.character(predict(qda(y ~ lat, data = tr), newdata = te)$class),
    nearest_mean = {
      mu <- tapply(tr$lat, y, mean)
      names(mu)[which.min(abs(mu - te$lat))]
    },
    threshold = {
      mu <- sort(tapply(tr$lat, y, mean)); ord <- names(mu); k <- length(ord)
      u <- sort(unique(tr$lat)); cuts <- (head(u, -1) + tail(u, -1)) / 2
      combos <- if (k == 2) matrix(cuts, ncol = 1) else t(combn(cuts, k - 1))
      acc <- apply(combos, 1, function(cc)
        mean(ord[findInterval(tr$lat, cc) + 1] == as.character(y)))
      cc <- combos[which.max(acc), ]
      ord[findInterval(te$lat, cc) + 1]
    },
    majority = names(which.max(table(y))))
}

methods <- c('logistic', 'lda', 'qda', 'nearest_mean', 'threshold', 'majority')
loo_out <- list(); ov_out <- list(); pair_out <- list()
for (mg_val in c('0', 'I', 'II', 'III', 'IV')) {
  cl <- mv2[[mg_val]]$F$clusters
  d <- merge(data.table(loc_key = names(cl), type = as.character(cl)),
             loc_coords, by = 'loc_key')
  n <- nrow(d); k <- uniqueN(d$type)
  for (m in methods) {
    pred <- vapply(seq_len(n), function(i)
      as.character(fit_predict(m, d[-i], d[i]))[1], character(1))
    loo_out[[paste(mg_val, m)]] <- data.table(
      mg = mg_val, k = k, n_loc = n, method = m,
      loo_accuracy = mean(pred == d$type))
  }
  loo_out[[paste(mg_val, 'base')]] <- data.table(
    mg = mg_val, k = k, n_loc = n, method = 'majority_baseline_fullsample',
    loo_accuracy = max(table(d$type)) / n)

  ## overlap
  rng <- diff(range(d$lat))
  prs <- combn(sort(unique(d$type)), 2, simplify = FALSE)
  ov <- sapply(prs, function(p) {
    a <- range(d[type == p[1], lat]); b <- range(d[type == p[2], lat])
    max(0, min(a[2], b[2]) - max(a[1], b[1])) / rng
  })
  ov_out[[mg_val]] <- data.table(mg = mg_val, k = k,
    pair = sapply(prs, paste, collapse = '-'), overlap_pct = 100 * ov,
    lat_min = min(d$lat), lat_max = max(d$lat))

  ## same-latitude pairs in different types
  ij <- t(combn(n, 2))
  pp <- data.table(mg = mg_val, loc_a = d$loc_key[ij[, 1]], loc_b = d$loc_key[ij[, 2]],
                   type_a = d$type[ij[, 1]], type_b = d$type[ij[, 2]],
                   lat_a = d$lat[ij[, 1]], lat_b = d$lat[ij[, 2]],
                   dlat = abs(d$lat[ij[, 1]] - d$lat[ij[, 2]]),
                   km = haversine_km(d$lat[ij[, 1]], d$lon[ij[, 1]],
                                     d$lat[ij[, 2]], d$lon[ij[, 2]]))
  pair_out[[mg_val]] <- pp[type_a != type_b & dlat <= 0.5][order(km)]
}
loo <- rbindlist(loo_out); ovl <- rbindlist(ov_out); prs <- rbindlist(pair_out)
fwrite(loo, file.path(OUT, 'Table_Latitude_LOO.csv'))
fwrite(ovl, file.path(OUT, 'Table_Latitude_Overlap.csv'))
fwrite(prs, file.path(OUT, 'Table_Latitude_Pairs.csv'))
print(dcast(loo, mg + k + n_loc ~ method, value.var = 'loo_accuracy'))
print(ovl); print(prs[, .N, by = mg])
