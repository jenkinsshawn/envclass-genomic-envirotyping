if (!exists('ROOT')) source('R/00_setup.R')

## ============================================================================
## 24_era_stability.R  (Online Resource 1, Table S8)
## Retention of environment-type assignments across three eras.
##
## Reconstructed: no script or output for Table S8 was recovered; the archived
## temporal analysis (R/01_marker_profiles.R) uses 10-year sliding windows.
## Online Resource text: "the record was divided into three eras of
## approximately equal length, profiles were recomputed within each era from
## the lines tested during that era, and cluster assignments were compared
## between eras as the fraction of locations retaining their assignment."
##
## Choices made here (stated because the text does not fix them):
##   eras      : primary 1993-2003, 2004-2013, 2014-2024 (11/10/11 years);
##               sensitivity 1993-2002/2003-2013/2014-2024 and
##               1993-2003/2004-2014/2015-2024
##   lines     : genotyped lines tested in the era; adjusted deviations
##               recomputed within the era (as for the sliding windows)
##   locations : published locations of the maturity group with >= MIN_LINES
##               genotyped lines in the era (primary 30, as for the full
##               analysis; sensitivity 20, as for the sliding windows)
##   clustering: published pipeline (variance weighting, correlation
##               distance, Ward's D2) cut at the published k
##   retention : fraction of locations with the same type after the best
##               one-to-one matching of type labels (labels are arbitrary).
##               (a) between each pair of eras, locations profiled in both;
##               (b) each era against the published full-period partition.
##   chance    : retention expected with the same group sizes and random
##               labels (999 permutations), since with two types matching
##               alone guarantees at least 50%.
##
## Output: out/Table_S8_EraRetention.csv, out/Table_S8_EraRetention_detail.csv,
##         out/Table_S8_EraLineOverlap.csv
## ============================================================================

library(data.table)
library(cluster)
set.seed(42)

nust <- fread(file.path(DATA, 'nust_ecw_subset.csv'))
nust[, loc_key := trimws(tolower(gsub(',.*$', '', Location)))]
gw <- readRDS(file.path(DATA, 'geno_worker_data.rds'))
geno_df <- gw$geno_df; rm(gw)
strain_names <- names(geno_df)[-(1:3)]
geno_mat <- as.matrix(geno_df[, -(1:3)]); rownames(geno_mat) <- as.character(geno_df$marker)
genotyped <- intersect(unique(nust$Strain), strain_names)
mv2 <- readRDS(file.path(DATA, 'marker_v2_results.rds'))

compute_adj_dev <- function(dat) {
  em <- dat[, .(env_mean = mean(YieldBuA, na.rm = TRUE)), by = Environment]
  sm <- dat[, .(strain_mean = mean(YieldBuA, na.rm = TRUE)), by = Strain]
  gm <- mean(dat$YieldBuA, na.rm = TRUE)
  dat <- merge(dat, em, by = 'Environment'); dat <- merge(dat, sm, by = 'Strain')
  dat[, adj_dev := YieldBuA - strain_mean - env_mean + gm]
  dat
}
compute_profiles <- function(dt, locs, min_n = 20) {
  prof <- matrix(NA_real_, length(locs), nrow(geno_mat), dimnames = list(locs, rownames(geno_mat)))
  for (i in seq_along(locs)) {
    ld <- dt[loc_key == locs[i], .(val = mean(value, na.rm = TRUE)), by = Strain]
    sh <- intersect(ld$Strain, strain_names)
    if (length(sh) < min_n) next
    rv <- ld[match(sh, Strain), val]
    gs <- geno_mat[, match(sh, strain_names), drop = FALSE]
    gsd <- apply(gs, 1, sd, na.rm = TRUE)
    gc <- sweep(gs, 1, rowMeans(gs, na.rm = TRUE))
    cors <- as.numeric(gc %*% (rv - mean(rv)) / (length(rv) - 1)) / (gsd * sd(rv))
    cors[gsd == 0 | !is.finite(cors)] <- 0
    prof[i, ] <- cors
  }
  prof[!is.na(prof[, 1]), , drop = FALSE]
}
cluster_k <- function(pf, k) {
  pf <- pf[, apply(pf, 2, var, na.rm = TRUE) > 0, drop = FALSE]
  pf_w <- sweep(pf, 2, sqrt(apply(pf, 2, var, na.rm = TRUE)), '*')
  cm <- cor(t(pf_w), use = 'pairwise.complete.obs'); cm[is.na(cm)] <- 0
  cutree(hclust(as.dist(1 - cm), method = 'ward.D2'), k)
}
perms <- function(v) if (length(v) <= 1) list(v) else
  do.call(c, lapply(seq_along(v), function(i) lapply(perms(v[-i]), function(p) c(v[i], p))))
## fraction agreeing after the best one-to-one relabelling of b
retention <- function(a, b) {
  la <- sort(unique(c(a, b))); best <- 0
  for (p in perms(la)) best <- max(best, mean(a == p[match(b, la)]))
  best
}
chance_retention <- function(a, b, n = 999) mean(replicate(n, retention(a, sample(b))))

ERA_SETS <- list(
  primary = list(c(1993, 2003), c(2004, 2013), c(2014, 2024)),
  alt_10_11_11 = list(c(1993, 2002), c(2003, 2013), c(2014, 2024)),
  alt_11_11_10 = list(c(1993, 2003), c(2004, 2014), c(2015, 2024)))
MIN_LINES_SET <- c(30, 20)

detail <- list(); overlap <- list()
for (mg_val in c('0', 'I', 'II', 'III', 'IV')) {
  ref <- mv2[[mg_val]]$F; k <- ref$k; locs <- names(ref$clusters)
  mm <- nust[MG == mg_val & Strain %in% genotyped & !is.na(YieldBuA),
             .(YieldBuA = mean(YieldBuA)), by = .(Strain, Environment, loc_key, Year)]
  for (es in names(ERA_SETS)) for (min_l in MIN_LINES_SET) {
    eras <- ERA_SETS[[es]]; cl_era <- list(); lines_era <- list()
    for (e in seq_along(eras)) {
      ed <- mm[Year >= eras[[e]][1] & Year <= eras[[e]][2]]
      lines_era[[e]] <- unique(ed$Strain)
      if (nrow(ed) < 100) next
      ed <- compute_adj_dev(ed)
      ln <- ed[loc_key %in% locs, .(n = uniqueN(Strain)), by = loc_key]
      el <- ln[n >= min_l, loc_key]
      if (length(el) < max(5, k + 2)) next
      pf <- compute_profiles(ed[, .(Strain, loc_key, value = adj_dev)], el, min_n = min_l)
      if (nrow(pf) < max(5, k + 2)) next
      cl_era[[as.character(e)]] <- cluster_k(pf, k)
    }
    lab <- function(e) sprintf('%d-%d', eras[[e]][1], eras[[e]][2])
    ## (b) era vs published partition
    for (e in names(cl_era)) {
      a <- cl_era[[e]]; b <- ref$clusters[names(a)]
      detail[[length(detail) + 1]] <- data.table(mg = mg_val, k = k, era_set = es,
        min_lines = min_l, comparison = 'era_vs_full', era_a = lab(as.integer(e)),
        era_b = '1993-2024', n_loc = length(a), retention = retention(a, b),
        chance = chance_retention(a, b))
    }
    ## (a) between eras
    if (length(cl_era) >= 2) for (pr in combn(names(cl_era), 2, simplify = FALSE)) {
      sh <- intersect(names(cl_era[[pr[1]]]), names(cl_era[[pr[2]]]))
      if (length(sh) < 4) next
      a <- cl_era[[pr[1]]][sh]; b <- cl_era[[pr[2]]][sh]
      detail[[length(detail) + 1]] <- data.table(mg = mg_val, k = k, era_set = es,
        min_lines = min_l, comparison = 'between_eras', era_a = lab(as.integer(pr[1])),
        era_b = lab(as.integer(pr[2])), n_loc = length(sh), retention = retention(a, b),
        chance = chance_retention(a, b))
    }
    if (es == 'primary' && min_l == 30) for (pr in combn(seq_along(eras), 2, simplify = FALSE))
      overlap[[length(overlap) + 1]] <- data.table(mg = mg_val, era_a = lab(pr[1]), era_b = lab(pr[2]),
        lines_a = length(lines_era[[pr[1]]]), lines_b = length(lines_era[[pr[2]]]),
        lines_shared = length(intersect(lines_era[[pr[1]]], lines_era[[pr[2]]])))
  }
}
det <- rbindlist(detail)
summ <- det[, .(n_comparisons = .N, n_loc_range = paste(range(n_loc), collapse = '-'),
                mean_retention_pct = 100 * mean(retention),
                range_retention_pct = paste(round(100 * range(retention)), collapse = '-'),
                mean_chance_pct = 100 * mean(chance)),
            by = .(mg, k, era_set, min_lines, comparison)]
fwrite(det, file.path(OUT, 'Table_S8_EraRetention_detail.csv'))
fwrite(summ, file.path(OUT, 'Table_S8_EraRetention.csv'))
fwrite(rbindlist(overlap), file.path(OUT, 'Table_S8_EraLineOverlap.csv'))
print(summ[era_set == 'primary'])
