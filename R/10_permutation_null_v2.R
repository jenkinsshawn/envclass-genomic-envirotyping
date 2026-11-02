if (!exists('ROOT')) source('R/00_setup.R')

## ============================================================================
## 10_permutation_null_v2.R
## Permutation test of cluster separation for the marker response profiles.
##
## Within each profiled location, the assignment of genotypes to line yields is
## shuffled (line identity permuted within locations), the marker response
## profiles are recomputed, the profiles are clustered exactly as in
## 01_marker_profiles.R (variance weighting, 1 - Pearson correlation, Ward.D2),
## and the best mean silhouette over k = 2 to 8 is recorded. The observed
## silhouette is the one stored in data/marker_v2_results.rds.
##
## Replaces the chromosome-level GWAS permutation in 05_permutation_null.R,
## which tested the earlier (v0.4) chromosome profiles, not the marker_v2
## profiles reported in the paper.
##
## PROFILE CONSTRUCTION (env var ENVCLASS_PERM_MODE)
##   "classification"  (default) profiles built exactly as in 01: adjusted
##                     deviations computed within each maturity group from
##                     strain x environment means of genotyped lines; this is
##                     the profile set that gives the reported clusters.
##   "pooled_adjdev"   as make_permnull() in 07_publication_figures.R and the
##                     2026-07-16 console run: adjusted deviations computed
##                     once over all maturity groups and all lines, location
##                     means taken over all maturity groups.
##   "pooled_raw"      as the 2026-06-29 publication-figure buffer that drew
##                     the first Fig3_Permutation_Null.tiff: raw yield.
##   The two pooled modes do not use the profile definition in the paper; they
##   are kept only to show where the earlier null values came from.
##
## SILHOUETTE OF EACH PERMUTATION (env var ENVCLASS_NULL_K)
##   "optimal"  (default) best mean silhouette over k = 2 to 8, the same rule
##              that chose k for the observed profiles (Methods wording).
##   "observed" mean silhouette at the k chosen for the observed profiles;
##              this reproduces the recovered perm_null_v2.rds (April 2026).
##
## P VALUE: (1 + number of permuted silhouettes >= observed) / (1 + n_perm).
## With 100 permutations the smallest attainable value is 1/101 = 0.0099.
##
## Run from the repository root (see run_permutation_null.ps1):
##   Rscript R/10_permutation_null_v2.R
## Environment variables: ENVCLASS_N_PERM (default 100), ENVCLASS_PERM_MODE,
## ENVCLASS_NULL_K.
##
## Writes out/permutation_null/Table_PermNull_<mode>.csv (one row per MG),
##        out/permutation_null/PermNull_<mode>_draws.csv (every permuted value),
##        out/permutation_null/perm_null_<mode>.rds
## ============================================================================

suppressPackageStartupMessages({ library(data.table); library(cluster) })

N_PERM    <- as.integer(Sys.getenv("ENVCLASS_N_PERM", "100"))
PERM_MODE <- Sys.getenv("ENVCLASS_PERM_MODE", "classification")
stopifnot(PERM_MODE %in% c("classification", "pooled_adjdev", "pooled_raw"))
NULL_K    <- Sys.getenv("ENVCLASS_NULL_K", "optimal")
stopifnot(NULL_K %in% c("optimal", "observed"))
TAG       <- if (NULL_K == "optimal") PERM_MODE else paste0(PERM_MODE, "_kobs")
SEED      <- 42
MAX_K     <- 8
MIN_N     <- 30           # min genotyped lines per location (as in 01)
MGS       <- c("0", "I", "II", "III", "IV")
OUT_PN    <- file.path(OUT, "permutation_null")
dir.create(OUT_PN, showWarnings = FALSE, recursive = TRUE)

t_start <- Sys.time()
cat(sprintf("Permutation null: mode = %s, k = %s, n_perm = %d\n", PERM_MODE, NULL_K, N_PERM))

## ---- data ---------------------------------------------------------------
nust <- fread(file.path(DATA, "nust_ecw_subset.csv"))
nust[, loc_key := trimws(tolower(gsub(",.*$", "", Location)))]
gw <- readRDS(file.path(DATA, "geno_worker_data.rds")); gdf <- gw$geno_df
strain_names <- names(gdf)[-(1:3)]
geno_mat <- as.matrix(gdf[, -(1:3)]); rownames(geno_mat) <- as.character(gdf$marker)
rm(gw, gdf); invisible(gc())
res <- readRDS(file.path(DATA, "marker_v2_results.rds"))
genotyped <- intersect(unique(nust$Strain), strain_names)

## per-MG table of location x strain values ("value") used for the profiles
POOLED <- NULL
value_table <- function(mg_val) {
  if (PERM_MODE == "classification") {
    d <- nust[MG == mg_val & Strain %in% genotyped & !is.na(YieldBuA)]
    d <- d[, .(YieldBuA = mean(YieldBuA, na.rm = TRUE)),
           by = .(Strain, Environment, loc_key, Year)]
    em <- d[, .(env_mean = mean(YieldBuA, na.rm = TRUE)), by = Environment]
    sm <- d[, .(strain_mean = mean(YieldBuA, na.rm = TRUE)), by = Strain]
    gm <- mean(d$YieldBuA, na.rm = TRUE)
    d <- merge(merge(d, em, by = "Environment"), sm, by = "Strain")
    d[, value := YieldBuA - strain_mean - env_mean + gm]
    return(d[, .(Strain, loc_key, value)])
  }
  if (is.null(POOLED)) {
    d <- copy(nust)
    if (!"Environment" %in% names(d)) d[, Environment := paste(loc_key, Year, sep = "_")]
    if (PERM_MODE == "pooled_adjdev") {
      gm <- mean(d$YieldBuA, na.rm = TRUE)
      d[, sm := ave(YieldBuA, Strain,      FUN = function(z) mean(z, na.rm = TRUE))]
      d[, em := ave(YieldBuA, Environment, FUN = function(z) mean(z, na.rm = TRUE))]
      d[, value := YieldBuA - sm - em + gm]
    } else d[, value := YieldBuA]
    POOLED <<- d[, .(Strain, loc_key, value)]
  }
  POOLED
}

## Centered genotype block and centered response for one location. For a
## permutation 'perm' of the genotype columns (as sample(sc) in 07), the
## profile gc[, perm] %*% rc equals gc %*% rc[order(perm)], so the genotype
## block is centered once and each permutation is a matrix-vector product.
location_block <- function(dt, loc) {
  ld <- dt[loc_key == loc, .(val = mean(value, na.rm = TRUE)), by = Strain]
  sh <- intersect(ld$Strain, strain_names)
  if (length(sh) < MIN_N) return(NULL)
  rv <- ld[match(sh, Strain), val]
  gs <- geno_mat[, match(sh, strain_names), drop = FALSE]
  gc <- gs - rowMeans(gs, na.rm = TRUE)
  nn <- rowSums(!is.na(gs))
  gsd <- sqrt(rowSums(gc^2, na.rm = TRUE) / (nn - 1))
  rc <- rv - mean(rv, na.rm = TRUE)
  list(gc = gc, gsd = gsd, rc = rc, rsd = sd(rv, na.rm = TRUE), n = length(rv))
}

profile_from <- function(b, rc) {
  cors <- as.numeric(b$gc %*% rc) / ((b$n - 1) * b$gsd * b$rsd)
  cors[b$gsd == 0 | !is.finite(cors)] <- 0
  cors
}

best_sil <- function(P, max_k = MAX_K, k_fixed = NULL) {
  if (nrow(P) < 5) return(NA_real_)
  mv <- apply(P, 2, var, na.rm = TRUE); P <- P[, mv > 0, drop = FALSE]
  P <- sweep(P, 2, sqrt(apply(P, 2, var, na.rm = TRUE)), "*")
  cm <- cor(t(P), use = "pairwise.complete.obs"); cm[is.na(cm)] <- 0
  d <- as.dist(1 - cm); hc <- hclust(d, method = "ward.D2")
  if (!is.null(k_fixed)) return(mean(silhouette(cutree(hc, k_fixed), d)[, 3]))
  mk <- min(max_k, nrow(P) - 1)
  max(sapply(2:mk, function(k) mean(silhouette(cutree(hc, k), d)[, 3])))
}

## ---- run -----------------------------------------------------------------
RNGkind("Mersenne-Twister", "Inversion", "Rejection")   # R default, made explicit
set.seed(SEED)
mg_seeds <- setNames(sample.int(.Machine$integer.max, length(MGS)), MGS)

summ <- list(); draws <- list(); store <- list()
for (mg_val in MGS) {
  r <- res[[mg_val]]$F
  if (is.null(r)) next
  locs <- names(r$clusters)
  dt <- value_table(mg_val)
  blocks <- lapply(setNames(locs, locs), function(l) location_block(dt, l))
  blocks <- blocks[!vapply(blocks, is.null, logical(1))]

  ## unpermuted profiles from this construction; in classification mode the
  ## best silhouette must equal the stored one
  P0 <- do.call(rbind, lapply(blocks, function(b) profile_from(b, b$rc)))
  sil_rebuilt <- best_sil(P0)

  set.seed(mg_seeds[[mg_val]])
  null <- numeric(N_PERM)
  for (p in seq_len(N_PERM)) {
    Pp <- do.call(rbind, lapply(blocks, function(b) {
      perm <- sample.int(b$n)                       # shuffle genotype columns
      profile_from(b, b$rc[order(perm)])
    }))
    null[p] <- best_sil(Pp, k_fixed = if (NULL_K == "observed") r$k else NULL)
  }
  null <- null[is.finite(null)]
  n_ge <- sum(null >= r$sil)
  summ[[mg_val]] <- data.table(
    mg = mg_val, mode = PERM_MODE, null_k = NULL_K, n_locs = length(blocks), k = r$k,
    observed = r$sil, observed_rebuilt = sil_rebuilt,
    null_mean = mean(null), null_sd = sd(null), null_max = max(null),
    ratio = r$sil / mean(null), n_perm = length(null), n_null_ge_obs = n_ge,
    p_plus1 = (1 + n_ge) / (1 + length(null)), p_raw = n_ge / length(null))
  draws[[mg_val]] <- data.table(mg = mg_val, perm = seq_along(null), null_sil = null)
  store[[mg_val]] <- list(observed = r$sil, null_sils = null)
  cat(sprintf("  MG %-3s  locs=%2d  obs=%.3f (rebuilt %.3f)  null mean=%.4f  ratio=%.2f  p=%.4f\n",
              mg_val, length(blocks), r$sil, sil_rebuilt, mean(null),
              r$sil / mean(null), (1 + n_ge) / (1 + length(null))))
}
summ <- rbindlist(summ); draws <- rbindlist(draws)
fwrite(summ,  file.path(OUT_PN, sprintf("Table_PermNull_%s.csv", TAG)))
fwrite(draws, file.path(OUT_PN, sprintf("PermNull_%s_draws.csv", TAG)))
saveRDS(list(summary = summ, store = store, mode = PERM_MODE, null_k = NULL_K, n_perm = N_PERM,
             seed = SEED, mg_seeds = mg_seeds),
        file.path(OUT_PN, sprintf("perm_null_%s.rds", TAG)))
## ---- Online Resource Fig. S1 (classification mode, optimal k only) --------
if (PERM_MODE == "classification" && NULL_K == "optimal" &&
    requireNamespace("ggplot2", quietly = TRUE)) {
  library(ggplot2)
  draws[, MG := factor(paste("MG", mg), levels = paste("MG", MGS))]
  obs <- summ[, .(MG = factor(paste("MG", mg), levels = paste("MG", MGS)), observed)]
  figS1 <- ggplot(draws, aes(MG, null_sil)) +
    geom_violin(fill = "grey80", colour = NA, width = 0.8) +
    geom_point(data = obs, aes(MG, observed), size = 2.2, colour = "black") +
    scale_y_continuous(limits = c(0, NA), expand = expansion(mult = c(0, 0.05))) +
    labs(x = NULL, y = "Mean silhouette width") +
    theme_classic(base_size = 9, base_family = "sans") +
    theme(axis.text = element_text(colour = "black"))
  ggsave(file.path(FIG, "FigS1_Permutation_Null.tiff"), figS1, width = 129, height = 90,
         units = "mm", dpi = 300, compression = "lzw", bg = "white")
  ggsave(file.path(FIG, "FigS1_Permutation_Null.png"), figS1, width = 129, height = 90,
         units = "mm", dpi = 300, bg = "white")
}
cat(sprintf("Done in %.1f min\n", as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
