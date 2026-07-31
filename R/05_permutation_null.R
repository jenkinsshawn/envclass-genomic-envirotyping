if (!exists('ROOT')) source('R/00_setup.R')

## ============================================================================
## PERMUTATION NULL: DEFINITIVE CONTROL FOR GENOMIC ENVIRONMENT CLASSIFICATION
## Shuffles location labels on GxE residuals, re-runs GWAS, tests whether
## real chromosome profiles have more structure than permuted profiles.
## Requires: full_results.rds with gxe_all saved (re-run v0.4 first)
## Runtime: ~3-4 hours (100 permutations x 124 locations / 6 cores)
## ============================================================================

library(data.table)
library(rrBLUP)
library(foreach)
library(doParallel)
library(cluster)

set.seed(42)

PATHS <- list(
  envclass_results = file.path(DATA, "full_results.rds"),
  geno_disk        = file.path(DATA, "geno_worker_data.rds"),
  nust_geno_cache  = file.path(DATA, "genotype_cache_nust.rds"),
  output_dir       = DATA,
  figure_dir       = FIG,
  table_dir        = OUT
)

PARAMS <- list(
  n_perm      = 100,
  n_cores     = 6,
  gwas_maf    = 0.05,
  min_strains = 30,
  cluster_method = "ward.D2"
)

chr_cols <- paste0("Gm", sprintf("%02d", 1:20))

cat("============================================================================\n")
cat("  PERMUTATION NULL TEST\n")
cat("  Started:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("============================================================================\n")

## ── Load data ──────────────────────────────────────────────────────────────

cat("\n[1] Loading data...\n")

full <- readRDS(PATHS$envclass_results)

if (is.null(full$gxe_all)) {
  stop("gxe_all not found in full_results.rds. Re-run v0.4 pipeline first.")
}

gxe_all <- full$gxe_all
loc_chr_mat <- full$loc_chr_mat
cl_gen <- full$cl_gen
within_mg <- full$within_mg_clusters

## SNP chromosome mapping
snp_info <- readRDS(PATHS$nust_geno_cache)$snp_info
snp_chr_map <- setNames(snp_info$Chr_Num, snp_info$SNP)

## Connected strains (from geno_disk)
gd <- readRDS(PATHS$geno_disk)
connected_strains <- colnames(gd$geno_df)[-(1:3)]

cat(sprintf("  GxE residuals: %d obs, %d environments\n",
            nrow(gxe_all), uniqueN(gxe_all$Environment)))
cat(sprintf("  Connected strains: %d\n", length(connected_strains)))

## Identify the all-years location x MG combos from the real analysis
loc_chr <- data.table(loc_mg = rownames(loc_chr_mat))
loc_chr[, c("loc", "mg") := tstrsplit(loc_mg, "_")]

cat(sprintf("  Location-MG combos to test: %d\n", nrow(loc_chr)))


## ── GWAS function ──────────────────────────────────────────────────────────

run_gwas_on_means <- function(strain_dt, geno_disk_path, min_maf) {
  gd <- readRDS(geno_disk_path)
  available <- intersect(strain_dt$Strain, colnames(gd$geno_df)[-(1:3)])
  if (length(available) < 30)
    return(list(n = length(available), status = "too_few", logp = NULL))
  
  sdt <- strain_dt[Strain %in% available]
  gsub <- gd$geno_df[, c("marker", "chrom", "pos", available), drop = FALSE]
  Ks <- gd$K[available, available]
  rownames(Ks) <- available; colnames(Ks) <- available
  
  gv <- as.matrix(gsub[, -(1:3)])
  maf <- apply(gv, 1, function(x) { p <- (mean(x, na.rm=TRUE)+1)/2; min(p,1-p) })
  gsub <- gsub[maf >= min_maf & !is.na(maf), ]
  if (nrow(gsub) < 50)
    return(list(n = length(available), status = "too_few_snps", logp = NULL))
  
  pdf <- data.frame(gid = sdt$Strain, trait = sdt$y,
                    check.names = FALSE, stringsAsFactors = FALSE)
  res <- tryCatch(
    rrBLUP::GWAS(pheno = pdf, geno = gsub, K = Ks,
                 n.PC = 0, min.MAF = min_maf, n.core = 1, plot = FALSE),
    error = function(e) NULL)
  
  if (is.null(res)) return(list(n = length(available), status = "error", logp = NULL))
  list(n = length(available), status = "success",
       logp = setNames(res$trait, res$marker))
}

logp_to_chr <- function(logp, chr_map) {
  chrs <- chr_map[names(logp)]
  valid <- !is.na(chrs) & chrs > 0 & chrs <= 20
  out <- setNames(rep(0, 20), chr_cols)
  if (sum(valid) == 0) return(out)
  cm <- tapply(logp[valid], chrs[valid], mean, na.rm = TRUE)
  out[paste0("Gm", sprintf("%02d", as.integer(names(cm))))] <- cm
  out
}


## ── Compute real-data clustering quality ───────────────────────────────────

cat("\n[2] Computing real-data clustering quality...\n")

real_metrics <- list()
for (mg_val in names(within_mg)) {
  wc <- within_mg[[mg_val]]
  if (is.null(wc) || nrow(wc) < 10) next
  
  mg_locs <- loc_chr[mg == mg_val]
  mg_mat <- loc_chr_mat[paste(mg_locs$loc, mg_locs$mg, sep = "_"), , drop = FALSE]
  mg_sc <- scale(mg_mat)
  mg_dist <- dist(mg_sc)
  k_mg <- length(unique(wc$cluster))
  
  hc <- hclust(mg_dist, method = PARAMS$cluster_method)
  cl <- cutree(hc, k = k_mg)
  sil <- mean(silhouette(cl, mg_dist)[, 3])
  
  ## Total between-cluster variance / total variance
  ss_total <- sum(mg_sc^2)
  ss_between <- sum(sapply(unique(cl), function(k) {
    nrow(mg_sc[cl == k, , drop = FALSE]) * sum(colMeans(mg_sc[cl == k, , drop = FALSE])^2)
  }))
  
  real_metrics[[mg_val]] <- data.table(
    mg = mg_val, k = k_mg, silhouette = sil,
    var_explained = ss_between / ss_total, n_locs = nrow(mg_mat)
  )
  
  cat(sprintf("  MG %s: k=%d, silhouette=%.3f, var_explained=%.3f\n",
              mg_val, k_mg, sil, ss_between / ss_total))
}


## ── Run permutation null ───────────────────────────────────────────────────

cat("\n[3] Running permutation null...\n")

perm_cache <- file.path(PATHS$output_dir, "permutation_null_results.rds")

if (file.exists(perm_cache)) {
  cat("  Loading cached permutation results...\n")
  perm_results <- readRDS(perm_cache)
} else {
  
  ## Strategy: for each permutation, shuffle loc_key labels on gxe_all
  ## (each strain's GxE residuals get randomly assigned to a different location)
  ## Then re-aggregate to strain x location means and re-run GWAS
  
  perm_sil <- matrix(NA, PARAMS$n_perm, length(real_metrics))
  colnames(perm_sil) <- names(real_metrics)
  
  perm_var <- matrix(NA, PARAMS$n_perm, length(real_metrics))
  colnames(perm_var) <- names(real_metrics)
  
  for (p in seq_len(PARAMS$n_perm)) {
    perm_start <- Sys.time()
    cat(sprintf("\n  Permutation %d / %d\n", p, PARAMS$n_perm))
    
    ## Shuffle location labels within each MG
    ## This breaks the location-specific GxE structure while preserving
    ## strain effects and location sample sizes
    gxe_perm <- copy(gxe_all)
    for (mg_val in names(within_mg)) {
      mg_idx <- gxe_perm$MG == mg_val
      if (sum(mg_idx) == 0) next
      ## Shuffle loc_key among observations within this MG
      gxe_perm[mg_idx, loc_key := sample(loc_key)]
    }
    
    ## Build task list: all-years per location x MG (same structure as real)
    tasks <- list()
    for (i in seq_len(nrow(loc_chr))) {
      loc <- loc_chr$loc[i]
      mg <- loc_chr$mg[i]
      loc_dat <- gxe_perm[loc_key == loc & MG == mg & Strain %in% connected_strains]
      if (nrow(loc_dat) == 0) next
      sa <- loc_dat[, .(y = mean(GxE_resid, na.rm = TRUE)), by = Strain]
      sa <- sa[!is.na(y)]
      if (nrow(sa) >= PARAMS$min_strains) {
        tasks[[length(tasks) + 1]] <- list(
          loc = loc, mg = mg, strain_data = sa)
      }
    }
    
    cat(sprintf("    Tasks: %d\n", length(tasks)))
    
    ## Run GWAS in parallel
    cl_par <- makeCluster(PARAMS$n_cores)
    registerDoParallel(cl_par)
    clusterEvalQ(cl_par, { library(data.table); library(rrBLUP) })
    
    gwas_perm <- foreach(
      t_idx = seq_along(tasks),
      .packages = c("data.table", "rrBLUP"),
      .export = c("run_gwas_on_means", "PARAMS", "PATHS", "tasks"),
      .errorhandling = "pass"
    ) %dopar% {
      tk <- tasks[[t_idx]]
      r <- run_gwas_on_means(tk$strain_data, PATHS$geno_disk, PARAMS$gwas_maf)
      list(loc = tk$loc, mg = tk$mg, status = r$status, logp = r$logp)
    }
    stopCluster(cl_par)
    
    ## Aggregate to chromosome level
    perm_chr <- rbindlist(lapply(gwas_perm, function(r) {
      if (inherits(r, "error") || r$status != "success") return(NULL)
      cv <- logp_to_chr(r$logp, snp_chr_map)
      data.table(loc = r$loc, mg = r$mg, t(as.matrix(cv)))
    }))
    
    ## Cluster within each MG and compute quality metrics
    for (mg_val in names(real_metrics)) {
      mg_perm <- perm_chr[mg == mg_val]
      if (nrow(mg_perm) < 5) next
      
      mg_mat_p <- as.matrix(mg_perm[, ..chr_cols])
      rownames(mg_mat_p) <- mg_perm$loc
      
      k_mg <- real_metrics[[mg_val]]$k
      if (k_mg >= nrow(mg_mat_p)) next
      
      mg_sc_p <- scale(mg_mat_p)
      ## Replace NaN/Inf from zero-variance columns after permutation
      mg_sc_p[is.nan(mg_sc_p) | is.infinite(mg_sc_p)] <- 0
      mg_dist_p <- dist(mg_sc_p)
      ## Skip if distance matrix contains NA (degenerate permutation)
      if (any(is.na(mg_dist_p) | is.nan(mg_dist_p) | is.infinite(mg_dist_p))) next
      hc_p <- hclust(mg_dist_p, method = PARAMS$cluster_method)
      cl_p <- cutree(hc_p, k = k_mg)
      
      perm_sil[p, mg_val] <- mean(silhouette(cl_p, mg_dist_p)[, 3])
      
      ss_total <- sum(mg_sc_p^2)
      ss_between <- sum(sapply(unique(cl_p), function(k) {
        nrow(mg_sc_p[cl_p == k, , drop = FALSE]) *
          sum(colMeans(mg_sc_p[cl_p == k, , drop = FALSE])^2)
      }))
      perm_var[p, mg_val] <- ss_between / ss_total
    }
    
    perm_elapsed <- difftime(Sys.time(), perm_start, units = "mins")
    cat(sprintf("    Elapsed: %.1f min\n", as.numeric(perm_elapsed)))
    
    ## Checkpoint every 10 permutations
    if (p %% 10 == 0) {
      saveRDS(list(perm_sil = perm_sil, perm_var = perm_var,
                   real_metrics = real_metrics, completed = p),
              file.path(PATHS$output_dir, "permutation_null_partial.rds"))
      cat("    Checkpointed\n")
    }
  }
  
  perm_results <- list(
    perm_sil = perm_sil,
    perm_var = perm_var,
    real_metrics = real_metrics
  )
  saveRDS(perm_results, perm_cache)
  cat("  Cached\n")
}


## ── Evaluate permutation null ──────────────────────────────────────────────

cat("\n[4] Permutation null results:\n\n")

cat("  ═══════════════════════════════════════════════════════════════\n")
cat("  PERMUTATION NULL TEST: IS LOCATION-SPECIFIC STRUCTURE REAL?\n")
cat("  ═══════════════════════════════════════════════════════════════\n\n")

for (mg_val in names(perm_results$real_metrics)) {
  real <- perm_results$real_metrics[[mg_val]]
  perm_s <- perm_results$perm_sil[, mg_val]
  perm_v <- perm_results$perm_var[, mg_val]
  perm_s <- perm_s[!is.na(perm_s)]
  perm_v <- perm_v[!is.na(perm_v)]
  
  if (length(perm_s) < 10) next
  
  p_sil <- mean(perm_s >= real$silhouette)
  p_var <- mean(perm_v >= real$var_explained)
  
  cat(sprintf("  MG %s (k=%d, n=%d locations):\n", mg_val, real$k, real$n_locs))
  cat(sprintf("    Silhouette:    real = %.3f, perm mean = %.3f (SD = %.3f), p = %.4f\n",
              real$silhouette, mean(perm_s), sd(perm_s), p_sil))
  cat(sprintf("    Var explained: real = %.3f, perm mean = %.3f (SD = %.3f), p = %.4f\n",
              real$var_explained, mean(perm_v), sd(perm_v), p_var))
  
  if (p_sil < 0.05) {
    cat("    >> SIGNIFICANT: Clustering is stronger than permutation null.\n")
    cat("    >> Location-specific genetic architecture is REAL.\n\n")
  } else {
    cat("    >> Not significant: clustering quality not distinguishable from null.\n\n")
  }
}

## Permutation null distribution figure
for (mg_val in names(perm_results$real_metrics)) {
  perm_s <- perm_results$perm_sil[, mg_val]
  perm_s <- perm_s[!is.na(perm_s)]
  if (length(perm_s) < 10) next
  
  real_sil <- perm_results$real_metrics[[mg_val]]$silhouette
  
  p_null <- ggplot(data.table(sil = perm_s), aes(x = sil)) +
    geom_histogram(bins = 30, fill = "gray70", color = "gray40") +
    geom_vline(xintercept = real_sil, color = "red", linewidth = 1.2) +
    annotate("text", x = real_sil, y = Inf, label = sprintf("Real = %.3f", real_sil),
             hjust = -0.1, vjust = 2, color = "red", size = 4) +
    labs(title = sprintf("MG %s: Permutation Null Distribution", mg_val),
         subtitle = sprintf("p = %.4f (%d permutations)",
                            mean(perm_s >= real_sil), length(perm_s)),
         x = "Mean Silhouette Width", y = "Count") +
    theme_minimal(base_size = 13)
  
  ggsave(file.path(PATHS$figure_dir, sprintf("Fig_PermNull_MG%s.png", mg_val)),
         p_null, width = 8, height = 5, dpi = 300)
}

cat("\n============================================================================\n")
cat("  PERMUTATION NULL COMPLETE\n")
cat("============================================================================\n")
