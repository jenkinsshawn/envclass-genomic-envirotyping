if (!exists('ROOT')) source('R/00_setup.R')

## ============================================================================
## GENOMIC ENVIRONMENT CLASSIFICATION v2: WITHIN-MG SIGNED MARKER PROFILES
## ============================================================================
## Method F + variance weighting (primary)
## Method G + variance weighting (sensitivity: responsive genotypes)
## Multi-trait post-hoc validation
## Temporal window analysis
## ============================================================================

library(data.table)
library(cluster)
library(vegan)
library(ggplot2)

set.seed(42)
pipeline_start <- Sys.time()

cat("============================================================================\n")
cat("  WITHIN-MG SIGNED MARKER ENVIRONMENT CLASSIFICATION v2\n")
cat("  Started:", format(pipeline_start, "%Y-%m-%d %H:%M:%S"), "\n")
cat("============================================================================\n")

PATHS <- list(
  nust_file   = file.path(DATA, "nust_ecw_subset.csv"),
  geno_worker = file.path(DATA, "geno_worker_data.rds"),
  output_dir  = DATA,
  figure_dir  = FIG,
  cache_dir   = file.path(DATA, "marker_cache_v2")
)
dir.create(PATHS$figure_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(PATHS$cache_dir, showWarnings = FALSE, recursive = TRUE)

PARAMS <- list(
  mgs              = c("0", "I", "II", "III", "IV"),
  traits           = c("WP_season", "TE_season", "HtWP_season", "HIWP_season",
                       "Height_in", "Maturity", "SeedSize", "Lodging"),
  min_loc_strains  = 30,
  min_strain_envs  = 5,
  fw_responsive    = 1.0,
  cluster_method   = "ward.D2",
  max_k            = 8,
  kw_alpha         = 0.05,
  consensus_reps   = 500,
  consensus_prop   = 0.80,
  temporal_window  = 10
)

## ============================================================================
## PHASE 1: LOAD DATA
## ============================================================================

cat("\n[1] Loading data...\n")

nust <- fread(PATHS$nust_file)
nust[, loc_key := trimws(tolower(gsub(",.*$", "", Location)))]
cat(sprintf("  NUST: %d obs\n", nrow(nust)))

geno_worker <- readRDS(PATHS$geno_worker)
geno_df <- geno_worker$geno_df
snp_info <- data.table(marker = as.character(geno_df$marker),
                       chr = geno_df$chrom, pos = geno_df$pos)
snp_info[, chr_name := paste0("Gm", sprintf("%02d", chr))]

strain_names <- names(geno_df)[-(1:3)]
geno_mat <- as.matrix(geno_df[, -(1:3)])
rownames(geno_mat) <- as.character(geno_df$marker)
marker_ids <- as.character(geno_df$marker)
n_markers <- length(marker_ids)

genotyped <- intersect(unique(nust$Strain), strain_names)
cat(sprintf("  Genotyped: %d strains, %d markers\n", length(genotyped), n_markers))

## ============================================================================
## HELPERS
## ============================================================================

compute_profiles <- function(dt, locs, markers = marker_ids, min_n = 20) {
  profiles <- matrix(NA_real_, nrow = length(locs), ncol = length(markers),
                     dimnames = list(locs, markers))
  n_per <- integer(length(locs)); names(n_per) <- locs
  
  for (i in seq_along(locs)) {
    loc <- locs[i]
    loc_dat <- dt[loc_key == loc, .(val = mean(value, na.rm = TRUE)), by = Strain]
    shared <- intersect(loc_dat$Strain, strain_names)
    n_per[i] <- length(shared)
    if (length(shared) < min_n) next
    
    rv <- loc_dat[match(shared, Strain), val]
    sc <- match(shared, strain_names)
    gs <- geno_mat[markers, sc, drop = FALSE]
    
    gsd <- apply(gs, 1, sd, na.rm = TRUE)
    rc <- rv - mean(rv, na.rm = TRUE)
    gc <- sweep(gs, 1, rowMeans(gs, na.rm = TRUE))
    cp <- gc %*% rc / (length(rv) - 1)
    rsd <- sd(rv, na.rm = TRUE)
    cors <- as.numeric(cp) / (gsd * rsd)
    cors[gsd == 0 | is.na(cors) | is.nan(cors)] <- 0
    profiles[i, ] <- cors
  }
  valid <- !is.na(profiles[, 1])
  list(profiles = profiles[valid, , drop = FALSE], n_strains = n_per[valid])
}

compute_adj_dev <- function(dat) {
  ## Compute adjusted deviations from raw data
  ## dat must have: Strain, Environment, loc_key, Year, YieldBuA
  em <- dat[, .(env_mean = mean(YieldBuA, na.rm = TRUE)), by = Environment]
  sm <- dat[, .(strain_mean = mean(YieldBuA, na.rm = TRUE)), by = Strain]
  gm <- mean(dat$YieldBuA, na.rm = TRUE)
  dat <- merge(dat, em, by = "Environment")
  dat <- merge(dat, sm, by = "Strain")
  dat[, adj_dev := YieldBuA - strain_mean - env_mean + gm]
  dat
}

cluster_full <- function(profiles, label, max_k = PARAMS$max_k) {
  if (nrow(profiles) < 5) { cat(sprintf("  [%s] Too few locations\n", label)); return(NULL) }
  
  ## Filter zero-variance markers
  mv <- apply(profiles, 2, var, na.rm = TRUE)
  pf <- profiles[, mv > 0, drop = FALSE]
  
  ## Variance weighting
  marker_weights <- sqrt(apply(pf, 2, var, na.rm = TRUE))
  pf_w <- sweep(pf, 2, marker_weights, "*")
  
  ## Correlation distance
  cm <- cor(t(pf_w), use = "pairwise.complete.obs")
  cm[is.na(cm)] <- 0
  d <- as.dist(1 - cm)
  
  ## Hierarchical clustering
  hc <- hclust(d, method = PARAMS$cluster_method)
  mk <- min(max_k, nrow(pf_w) - 1)
  if (mk < 2) return(NULL)
  
  sils <- numeric(mk - 1)
  for (k in 2:mk) { cl <- cutree(hc, k); sils[k-1] <- mean(silhouette(cl, d)[, 3]) }
  best_k <- which.max(sils) + 1
  best_sil <- max(sils)
  clusters <- cutree(hc, k = best_k)
  
  cat(sprintf("  [%s] k=%d, sil=%.3f, n=%d, markers=%d\n",
              label, best_k, best_sil, nrow(pf_w), ncol(pf_w)))
  for (k in 2:min(6, mk))
    cat(sprintf("    k=%d: %.3f%s\n", k, sils[k-1], ifelse(k==best_k," ***","")))
  for (k in sort(unique(clusters)))
    cat(sprintf("    C%d: %s\n", k, paste(names(clusters)[clusters==k], collapse=", ")))
  
  ## Consensus
  n_l <- nrow(pf_w)
  coocc <- matrix(0, n_l, n_l); samp <- matrix(0, n_l, n_l)
  for (rep in seq_len(PARAMS$consensus_reps)) {
    idx <- sample(n_l, floor(n_l * PARAMS$consensus_prop))
    sd2 <- as.dist(as.matrix(d)[idx, idx])
    sh <- hclust(sd2, method = PARAMS$cluster_method)
    sc <- cutree(sh, k = best_k)
    for (a in seq_along(idx)) for (b in seq_along(idx)) {
      i <- idx[a]; j <- idx[b]
      samp[i,j] <- samp[i,j] + 1
      if (sc[a] == sc[b]) coocc[i,j] <- coocc[i,j] + 1
    }
  }
  cons <- coocc / pmax(samp, 1); diag(cons) <- 1
  wc <- numeric(best_k)
  for (k in seq_len(best_k)) {
    m <- which(clusters == k)
    wc[k] <- if (length(m) > 1) mean(cons[m,m][lower.tri(cons[m,m])]) else 1
  }
  cat(sprintf("    Consensus: %.3f\n", mean(wc)))
  
  ## KW discriminating markers
  kw <- rbindlist(lapply(colnames(pf), function(m) {
    v <- pf[,m]
    if (length(unique(v[!is.na(v)])) < 2) return(data.table(marker=m, kw_p=NA_real_))
    p <- tryCatch(kruskal.test(v ~ as.factor(clusters))$p.value, error=function(e) NA_real_)
    data.table(marker=m, kw_p=p)
  }))
  kw[, p_adj := p.adjust(kw_p, method="BH")]
  kw <- merge(kw, snp_info, by="marker", all.x=TRUE)
  kw <- kw[order(kw_p)]
  n_sig <- sum(kw$p_adj < PARAMS$kw_alpha, na.rm=TRUE)
  cat(sprintf("    Discriminating markers: %d\n", n_sig))
  
  if (n_sig > 0) {
    chr_disc <- kw[p_adj < PARAMS$kw_alpha, .N, by=chr_name][order(-N)]
    cat("    Top chromosomes: ")
    cat(paste(head(chr_disc$chr_name, 5), collapse=", "), "\n")
  }
  
  list(k=best_k, sil=best_sil, clusters=clusters, hc=hc, dist=d,
       consensus=cons, within_consensus=wc, kw=kw, n_sig=n_sig, profiles=pf)
}

## ============================================================================
## PHASE 2: WITHIN-MG ANALYSIS
## ============================================================================

cat("\n[2] Within-MG analysis...\n")

mg_results <- list()

for (mg_val in PARAMS$mgs) {
  cat(sprintf("\n======== MG %s ========\n", mg_val))
  
  mg_dat <- nust[MG == mg_val & Strain %in% genotyped & !is.na(YieldBuA)]
  mg_means <- mg_dat[, .(YieldBuA = mean(YieldBuA, na.rm = TRUE)),
                     by = .(Strain, Environment, loc_key, Year)]
  mg_means <- compute_adj_dev(mg_means)
  
  loc_n <- mg_means[, .(n = uniqueN(Strain)), by = loc_key]
  eligible <- loc_n[n >= PARAMS$min_loc_strains, loc_key]
  cat(sprintf("  Obs: %d, Strains: %d, Locations: %d (eligible: %d)\n",
              nrow(mg_means), uniqueN(mg_means$Strain), uniqueN(mg_means$loc_key),
              length(eligible)))
  
  if (length(eligible) < 5) { cat("  Skipping\n"); next }
  
  ## FW stability
  strain_env_n <- mg_means[, .(n_envs = uniqueN(Environment)), by = Strain]
  fw_eligible <- strain_env_n[n_envs >= PARAMS$min_strain_envs, Strain]
  env_m <- mg_means[, .(em = mean(YieldBuA)), by = Environment]
  mg_means_fw <- merge(mg_means, env_m, by = "Environment", suffixes = c("", ".fw"))
  
  fw <- mg_means_fw[Strain %in% fw_eligible, {
    if (.N >= 5) .(slope = coef(lm(YieldBuA ~ em))[2]) else .(slope = NA_real_)
  }, by = Strain]
  fw <- fw[!is.na(slope)]
  responsive <- fw[slope > PARAMS$fw_responsive, Strain]
  stable <- fw[slope < PARAMS$fw_responsive, Strain]
  cat(sprintf("  FW: %d responsive, %d stable\n", length(responsive), length(stable)))
  
  ## ---- Method F: All genotypes ----
  cat("\n  --- Method F (all, var-weighted) ---\n")
  pf_F <- compute_profiles(mg_means[, .(Strain, loc_key, value = adj_dev)], eligible)
  res_F <- cluster_full(pf_F$profiles, paste("MG", mg_val, "F"))
  
  ## ---- Method G: Responsive only ----
  cat("\n  --- Method G (responsive, var-weighted) ---\n")
  pf_G <- compute_profiles(
    mg_means[Strain %in% responsive, .(Strain, loc_key, value = adj_dev)], eligible)
  res_G <- cluster_full(pf_G$profiles, paste("MG", mg_val, "G"))
  
  ## ---- Multi-trait post-hoc ----
  if (!is.null(res_F) && res_F$k >= 2) {
    cat("\n  --- Multi-trait post-hoc ---\n")
    trait_disc_counts <- list()
    
    for (trait in PARAMS$traits) {
      if (!(trait %in% names(mg_dat))) next
      tr_dat <- mg_dat[!is.na(get(trait))]
      if (nrow(tr_dat) < 100) next
      
      tr_means <- tr_dat[, .(trait_val = mean(get(trait), na.rm = TRUE)),
                         by = .(Strain, loc_key)]
      tr_sm <- tr_means[, .(sm = mean(trait_val)), by = Strain]
      tr_lm <- tr_means[, .(lm = mean(trait_val)), by = loc_key]
      tr_gm <- mean(tr_means$trait_val, na.rm = TRUE)
      tr_means <- merge(tr_means, tr_sm, by = "Strain")
      tr_means <- merge(tr_means, tr_lm, by = "loc_key")
      tr_means[, adj_dev_tr := trait_val - sm - lm + tr_gm]
      
      f_locs <- names(res_F$clusters)
      pf_tr <- compute_profiles(tr_means[, .(Strain, loc_key, value = adj_dev_tr)], f_locs)
      
      shared_locs <- intersect(rownames(pf_tr$profiles), names(res_F$clusters))
      if (length(shared_locs) < 5) next
      
      cl_sub <- res_F$clusters[shared_locs]
      tr_sub <- pf_tr$profiles[shared_locs, , drop = FALSE]
      
      n_disc <- sum(sapply(seq_len(ncol(tr_sub)), function(m) {
        v <- tr_sub[, m]
        if (length(unique(v[!is.na(v)])) < 2) return(FALSE)
        p <- tryCatch(kruskal.test(v ~ as.factor(cl_sub))$p.value, error = function(e) NA)
        !is.na(p) && p < 0.05
      }))
      
      trait_disc_counts[[trait]] <- n_disc
      cat(sprintf("    %s: %d discriminating markers\n", trait, n_disc))
    }
    
    if (length(trait_disc_counts) > 0) {
      tdc <- data.table(trait = names(trait_disc_counts),
                        n_disc = unlist(trait_disc_counts))[order(-n_disc)]
      cat("\n    Trait ranking:\n")
      print(tdc)
    }
  }
  
  ## ---- Temporal window analysis ----
  cat("\n  --- Temporal windows ---\n")
  years <- sort(unique(mg_means$Year))
  ws <- PARAMS$temporal_window
  window_starts <- years[1]:(max(years) - ws + 1)
  
  window_profiles <- list()
  for (w_start in window_starts) {
    w_end <- w_start + ws - 1
    win_dat <- mg_means[Year >= w_start & Year <= w_end,
                        .(Strain, Environment, loc_key, Year, YieldBuA)]
    if (nrow(win_dat) < 100) next
    
    ## Recompute adjusted deviations within this window
    win_dat <- compute_adj_dev(win_dat)
    
    win_locs <- win_dat[, .(n = uniqueN(Strain)), by = loc_key]
    win_eligible <- intersect(eligible, win_locs[n >= 20, loc_key])
    if (length(win_eligible) < 5) next
    
    wp <- compute_profiles(win_dat[, .(Strain, loc_key, value = adj_dev)],
                           win_eligible, min_n = 20)
    if (nrow(wp$profiles) >= 5)
      window_profiles[[as.character(w_start)]] <- wp$profiles
  }
  
  cat(sprintf("    Valid %dyr windows: %d\n", ws, length(window_profiles)))
  
  res_temporal <- NULL
  if (length(window_profiles) >= 3) {
    all_wlocs <- unique(unlist(lapply(window_profiles, rownames)))
    loc_counts <- table(unlist(lapply(window_profiles, rownames)))
    stable_locs <- names(loc_counts[loc_counts >= length(window_profiles) * 0.5])
    
    if (length(stable_locs) >= 10) {
      avg_prof <- matrix(0, nrow = length(stable_locs), ncol = n_markers,
                         dimnames = list(stable_locs, marker_ids))
      n_win_per <- integer(length(stable_locs)); names(n_win_per) <- stable_locs
      
      for (wp in window_profiles) {
        sh <- intersect(stable_locs, rownames(wp))
        avg_prof[sh, ] <- avg_prof[sh, ] + wp[sh, ]
        n_win_per[sh] <- n_win_per[sh] + 1
      }
      avg_prof <- avg_prof / pmax(n_win_per, 1)
      
      cat(sprintf("    Averaged: %d locations across %d windows\n",
                  nrow(avg_prof), length(window_profiles)))
      res_temporal <- cluster_full(avg_prof, paste("MG", mg_val, "Temporal"))
    }
  }
  
  mg_results[[mg_val]] <- list(F = res_F, G = res_G, temporal = res_temporal,
                               n_responsive = length(responsive),
                               n_stable = length(stable), n_locs = length(eligible))
}

## ============================================================================
## PHASE 3: GEOGRAPHY TESTS
## ============================================================================

cat("\n[3] Geography tests...\n")

nust_coords <- fread(PATHS$nust_file, select = c("Location","Lat","Long"))
nust_coords[, loc_key := trimws(tolower(gsub(",.*$", "", Location)))]
loc_coords <- nust_coords[, .(lat=median(Lat,na.rm=TRUE), lon=median(Long,na.rm=TRUE)), by=loc_key]

for (mg_val in names(mg_results)) {
  res <- mg_results[[mg_val]]$F
  if (is.null(res)) next
  
  cl_locs <- data.table(loc = names(res$clusters))
  cl_locs <- merge(cl_locs, loc_coords, by.x="loc", by.y="loc_key", all.x=TRUE)
  cl_locs <- cl_locs[!is.na(lat) & !is.na(lon)]
  if (nrow(cl_locs) < 5) next
  
  geo_d <- dist(cl_locs[, .(lat, lon)])
  gen_idx <- match(cl_locs$loc, rownames(res$profiles))
  gen_d <- as.dist(as.matrix(res$dist)[gen_idx, gen_idx])
  
  mt <- vegan::mantel(gen_d, geo_d, permutations=999)
  cat(sprintf("  MG %s: Mantel r=%.3f, p=%.3f\n", mg_val, mt$statistic, mt$signif))
}

## ============================================================================
## PHASE 4: FIGURES
## ============================================================================

cat("\n[4] Figures...\n")

for (mg_val in names(mg_results)) {
  res <- mg_results[[mg_val]]$F
  if (is.null(res)) next
  
  ## Dendrogram
  png(file.path(PATHS$figure_dir, sprintf("Fig_Dendro_MG%s.png", mg_val)),
      width=1400, height=700, res=150)
  plot(res$hc, labels=names(res$clusters), cex=0.6,
       main=sprintf("MG %s (k=%d, sil=%.3f, consensus=%.3f)",
                    mg_val, res$k, res$sil, mean(res$within_consensus)))
  rect.hclust(res$hc, k=res$k, border=2:(res$k+1))
  dev.off()
  
  ## MDS
  mds <- cmdscale(res$dist, k=2)
  mds_dt <- data.table(MDS1=mds[,1], MDS2=mds[,2],
                       loc=rownames(mds), cluster=as.factor(res$clusters[rownames(mds)]))
  
  p <- ggplot(mds_dt, aes(MDS1, MDS2, color=cluster)) +
    geom_point(size=3)
  if (requireNamespace("ggrepel", quietly=TRUE))
    p <- p + ggrepel::geom_text_repel(aes(label=loc), size=2.5, max.overlaps=25, show.legend=FALSE)
  p <- p + labs(title=sprintf("MG %s Environment Types (k=%d, sil=%.3f)",
                              mg_val, res$k, res$sil)) + theme_minimal()
  ggsave(file.path(PATHS$figure_dir, sprintf("Fig_MDS_MG%s.png", mg_val)),
         p, width=12, height=9, dpi=300)
  
  ## Consensus heatmap
  png(file.path(PATHS$figure_dir, sprintf("Fig_Consensus_MG%s.png", mg_val)),
      width=1000, height=900, res=150)
  ord <- order(res$clusters)
  heatmap(res$consensus[ord,ord], Rowv=NA, Colv=NA, scale="none",
          col=colorRampPalette(c("white","steelblue","darkblue"))(100),
          main=sprintf("MG %s Consensus (k=%d)", mg_val, res$k))
  dev.off()
}
cat("  Figures saved\n")

## ============================================================================
## SAVE AND SUMMARY
## ============================================================================

cat("\n[5] Saving...\n")
saveRDS(mg_results, file.path(PATHS$output_dir, "marker_v2_results.rds"))

elapsed <- difftime(Sys.time(), pipeline_start, units="mins")

cat("\n============================================================================\n")
cat("  WITHIN-MG CLASSIFICATION v2 SUMMARY\n")
cat("============================================================================\n")
cat(sprintf("  Runtime: %.1f min\n", as.numeric(elapsed)))

cat("\n  Method F (adjusted deviation, var-weighted, all genotypes):\n")
for (mg_val in names(mg_results)) {
  r <- mg_results[[mg_val]]$F
  if (!is.null(r))
    cat(sprintf("    MG %s: k=%d, sil=%.3f, consensus=%.3f, discrim=%d, n=%d\n",
                mg_val, r$k, r$sil, mean(r$within_consensus), r$n_sig,
                length(r$clusters)))
}

cat("\n  Method G (responsive genotypes, var-weighted):\n")
for (mg_val in names(mg_results)) {
  r <- mg_results[[mg_val]]$G
  if (!is.null(r))
    cat(sprintf("    MG %s: k=%d, sil=%.3f, consensus=%.3f, discrim=%d, n=%d\n",
                mg_val, r$k, r$sil, mean(r$within_consensus), r$n_sig,
                length(r$clusters)))
}

cat("\n  Temporal (10yr window averaged):\n")
for (mg_val in names(mg_results)) {
  r <- mg_results[[mg_val]]$temporal
  if (!is.null(r))
    cat(sprintf("    MG %s: k=%d, sil=%.3f, consensus=%.3f, discrim=%d, n=%d\n",
                mg_val, r$k, r$sil, mean(r$within_consensus), r$n_sig,
                length(r$clusters)))
}

cat("\n  For comparison:\n")
cat("    Chromosome-level (v0.4):     sil=0.018\n")
cat("    Cross-MG signed markers (v1): sil=0.112\n")

cat("\n============================================================================\n")
cat("  COMPLETE\n")
cat("============================================================================\n")
