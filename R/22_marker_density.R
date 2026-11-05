if (!exists('ROOT')) source('R/00_setup.R')

## ============================================================================
## 22_marker_density.R  (Online Resource 1, Table S3)
## Marker density sensitivity: random subsets of 2,000, 1,000 and 500 markers,
## 5 replicates each, compared with the full 4,625-marker partition by the
## adjusted Rand index.
##
## The only recovered density code (robustness module v0.6, section 7c)
## clusters 20 per-chromosome means of GWAS -log10(p) from the earlier v0.4
## pipeline, whose partitions (e.g. k = 6 in MG I) are not the published
## ones. This script applies the published pipeline (R/01_marker_profiles.R,
## cluster_full: variance weighting, correlation distance, Ward's D2) to the
## archived Method F profiles in data/marker_v2_results.rds.
##
## Each profile value is the correlation of one marker with the adjusted
## deviations at one location and does not depend on the other markers, so
## recomputing profiles from a marker subset is the same as selecting those
## columns of the full profile matrix.
##
## Two comparisons per subset: (a) k fixed at the published k; (b) k chosen
## again by maximum mean silhouette (k = 2..8), as in the full analysis.
##
## Output: out/Table_S1_MarkerDensity_reps.csv, out/Table_S1_MarkerDensity.csv
## ============================================================================

library(data.table)
library(cluster)
library(mclust)
set.seed(42)

mv2 <- readRDS(file.path(DATA, 'marker_v2_results.rds'))
GRID <- c(2000, 1000, 500); REPS <- 5; MAX_K <- 8

cluster_subset <- function(pf, k_fixed) {
  mv <- apply(pf, 2, var, na.rm = TRUE)
  pf <- pf[, mv > 0, drop = FALSE]
  pf_w <- sweep(pf, 2, sqrt(apply(pf, 2, var, na.rm = TRUE)), '*')
  cm <- cor(t(pf_w), use = 'pairwise.complete.obs'); cm[is.na(cm)] <- 0
  d <- as.dist(1 - cm)
  hc <- hclust(d, method = 'ward.D2')
  mk <- min(MAX_K, nrow(pf_w) - 1)
  sils <- sapply(2:mk, function(k) mean(silhouette(cutree(hc, k), d)[, 3]))
  list(fixed = cutree(hc, k_fixed), best = cutree(hc, which.max(sils) + 1),
       k_best = which.max(sils) + 1)
}

rows <- list()
for (mg_val in c('0', 'I', 'II', 'III', 'IV')) {
  res <- mv2[[mg_val]]$F
  pf <- res$profiles; ref <- res$clusters[rownames(pf)]
  ## check: the full marker set reproduces the published partition
  full <- cluster_subset(pf, res$k)
  stopifnot(adjustedRandIndex(full$fixed, ref) == 1)
  for (n_mk in GRID) for (r in seq_len(REPS)) {
    sub <- cluster_subset(pf[, sample(ncol(pf), n_mk)], res$k)
    rows[[length(rows) + 1]] <- data.table(
      mg = mg_val, k_published = res$k, n_markers = n_mk, replicate = r,
      ari_fixed_k = adjustedRandIndex(sub$fixed, ref),
      k_reselected = sub$k_best,
      ari_reselected_k = adjustedRandIndex(sub$best, ref))
  }
}
reps <- rbindlist(rows)
tab <- reps[, .(n_loc = NA_integer_, mean_ari_fixed_k = mean(ari_fixed_k),
                min_ari_fixed_k = min(ari_fixed_k),
                mean_ari_reselected_k = mean(ari_reselected_k),
                min_ari_reselected_k = min(ari_reselected_k),
                k_reselected = paste(sort(unique(k_reselected)), collapse = '/')),
            by = .(mg, k_published, n_markers)]
tab[, n_loc := sapply(mg, function(m) length(mv2[[m]]$F$clusters))]
fwrite(reps, file.path(OUT, 'Table_S1_MarkerDensity_reps.csv'))
fwrite(tab, file.path(OUT, 'Table_S1_MarkerDensity.csv'))
print(tab)

## ---- sensitivity: 100 replicates per density (5 replicates give unstable
## means because a single misplaced location moves the index by 0.1 to 0.3)
set.seed(7); rows <- list()
for (mg_val in c('0', 'I', 'II', 'III', 'IV')) {
  res <- mv2[[mg_val]]$F; pf <- res$profiles; ref <- res$clusters[rownames(pf)]
  for (n_mk in GRID) for (r in 1:100) {
    sub <- cluster_subset(pf[, sample(ncol(pf), n_mk)], res$k)
    rows[[length(rows) + 1]] <- data.table(mg = mg_val, n_markers = n_mk, replicate = r,
      ari_fixed_k = adjustedRandIndex(sub$fixed, ref),
      ari_reselected_k = adjustedRandIndex(sub$best, ref))
  }
}
r100 <- rbindlist(rows)
t100 <- r100[, .(mean_ari_fixed_k = mean(ari_fixed_k), median_ari_fixed_k = median(ari_fixed_k),
                 p10_ari_fixed_k = quantile(ari_fixed_k, 0.1, names = FALSE),
                 mean_ari_reselected_k = mean(ari_reselected_k)), by = .(mg, n_markers)]
fwrite(r100, file.path(OUT, 'Table_S1_MarkerDensity_100reps.csv'))
fwrite(t100, file.path(OUT, 'Table_S1_MarkerDensity_100reps_summary.csv'))
print(t100)
