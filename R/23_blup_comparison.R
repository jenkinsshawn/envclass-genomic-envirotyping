if (!exists('ROOT')) source('R/00_setup.R')

## ============================================================================
## 23_blup_comparison.R  (Online Resource 1, Table S1)
## Clusters from unshrunk adjusted deviations (published) compared with
## clusters from mixed-model residuals and BLUPs.
##
## Reconstructed: no script or output for Table S1 was recovered. The mixed
## model follows the recovered v0.4 pipeline (EnvClass_Pipeline_v0.4.R,
## Stage 1): YieldBuA ~ factor(Year) + (1|Strain) + (1|Environment), REML,
## and the GxE value is the residual y - fitted (fitted includes the shrunken
## line and environment BLUPs).
##
## Variants (column 'input'):
##   adj_dev     : published input (unshrunk adjusted deviations), reference
##   mm_resid    : residuals of the v0.4 mixed model, genotyped lines
##                 (same line-environment means as R/01_marker_profiles.R)
##   mm_gxe_blup : BLUP of a line-by-location term from
##                 YieldBuA ~ factor(Year) + (1|Strain) + (1|Environment) +
##                 (1|Strain:loc_key), i.e. shrunken line-by-location effects
##   mm_resid_alllines : as mm_resid, but the model is fitted on all lines of
##                 the maturity group (genotyped or not), as in v0.4
## For each input: marker profiles at the published locations, variance
## weighting, correlation distance, Ward's D2 (as cluster_full in
## R/01_marker_profiles.R, without the subsampling consensus); k chosen by
## maximum mean silhouette (2..8). Reported: silhouette at the chosen k and at
## the published k, its reduction relative to adj_dev, and the adjusted Rand
## index with the published partition.
##
## Output: out/Table_S7_BLUP_Comparison.csv
## ============================================================================

library(data.table)
library(cluster)
library(mclust)
library(lme4)
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

## As compute_profiles() in R/01_marker_profiles.R
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

cluster_profiles <- function(pf, k_pub, max_k = 8) {
  pf <- pf[, apply(pf, 2, var, na.rm = TRUE) > 0, drop = FALSE]
  pf_w <- sweep(pf, 2, sqrt(apply(pf, 2, var, na.rm = TRUE)), '*')
  cm <- cor(t(pf_w), use = 'pairwise.complete.obs'); cm[is.na(cm)] <- 0
  d <- as.dist(1 - cm); hc <- hclust(d, method = 'ward.D2')
  mk <- min(max_k, nrow(pf_w) - 1)
  sils <- sapply(2:mk, function(k) mean(silhouette(cutree(hc, k), d)[, 3]))
  kb <- which.max(sils) + 1
  list(k = kb, sil = max(sils), sil_kpub = sils[k_pub - 1],
       cl = cutree(hc, kb), cl_kpub = cutree(hc, k_pub))
}

ctrl <- lmerControl(optimizer = 'bobyqa', optCtrl = list(maxfun = 50000),
                    calc.derivs = FALSE)
MGS <- if (exists('S7_MGS')) S7_MGS else c('0', 'I', 'II', 'III', 'IV')
rows <- list()
for (mg_val in MGS) {
  t0 <- Sys.time()
  ref <- mv2[[mg_val]]$F; locs <- names(ref$clusters); k_pub <- ref$k
  mm <- nust[MG == mg_val & Strain %in% genotyped & !is.na(YieldBuA),
             .(YieldBuA = mean(YieldBuA)), by = .(Strain, Environment, loc_key, Year)]
  mm <- compute_adj_dev(mm)
  m1 <- lmer(YieldBuA ~ factor(Year) + (1 | Strain) + (1 | Environment),
             data = mm, REML = TRUE, control = ctrl)
  mm[, mm_resid := YieldBuA - fitted(m1)]
  m2 <- lmer(YieldBuA ~ factor(Year) + (1 | Strain) + (1 | Environment) + (1 | Strain:loc_key),
             data = mm, REML = TRUE, control = ctrl)
  re <- ranef(m2)$`Strain:loc_key`
  mm[, mm_gxe_blup := re[paste(Strain, loc_key, sep = ':'), 1]]
  vc <- as.data.frame(VarCorr(m2))
  ## sensitivity: v0.4 fitted the model on all lines of the maturity group,
  ## genotyped or not; residuals of genotyped lines then enter the profiles
  ma <- nust[MG == mg_val & !is.na(YieldBuA),
             .(YieldBuA = mean(YieldBuA)), by = .(Strain, Environment, loc_key, Year)]
  m3 <- lmer(YieldBuA ~ factor(Year) + (1 | Strain) + (1 | Environment),
             data = ma, REML = TRUE, control = ctrl)
  ma[, mm_resid_alllines := YieldBuA - fitted(m3)]
  mm <- merge(mm, ma[, .(Strain, Environment, mm_resid_alllines)],
              by = c('Strain', 'Environment'), all.x = TRUE)
  sil_ref <- NA_real_
  for (inp in c('adj_dev', 'mm_resid', 'mm_gxe_blup', 'mm_resid_alllines')) {
    pf <- compute_profiles(mm[, .(Strain, loc_key, value = get(inp))], locs)
    cr <- cluster_profiles(pf, k_pub)
    rf <- ref$clusters[rownames(pf)]
    if (inp == 'adj_dev') { sil_ref <- cr$sil; sil_ref_kpub <- cr$sil_kpub }
    rows[[length(rows) + 1]] <- data.table(
      mg = mg_val, input = inp, n_loc = nrow(pf), k_published = k_pub,
      k_chosen = cr$k, sil_chosen_k = cr$sil, sil_published_k = cr$sil_kpub,
      sil_reduction_pct = 100 * (1 - cr$sil / sil_ref),
      sil_reduction_pct_published_k = 100 * (1 - cr$sil_kpub / sil_ref_kpub),
      ari_chosen_k = adjustedRandIndex(cr$cl, rf),
      ari_published_k = adjustedRandIndex(cr$cl_kpub, rf),
      var_line_by_location = vc$vcov[vc$grp == 'Strain:loc_key'],
      var_residual = vc$vcov[vc$grp == 'Residual'])
  }
  cat(sprintf('MG %s done in %.1f min\n', mg_val,
              as.numeric(difftime(Sys.time(), t0, units = 'mins'))))
}
s7 <- rbindlist(rows)
fwrite(s7, file.path(OUT, if (exists('S7_OUTFILE')) S7_OUTFILE else 'Table_S7_BLUP_Comparison.csv'))
print(s7)
