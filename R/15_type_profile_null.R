if (!exists('ROOT')) source('R/00_setup.R')

## ============================================================================
## 15_type_profile_null.R
## Correlation between the mean marker profiles of the two environment types
## (MG I, II and IV), compared with the same quantity after line yields are
## reassigned at random among the genotyped lines within each location, the
## profiles recomputed and clustered into two groups as in the main analysis.
## Because the types are found by clustering the profiles, a negative
## correlation between type means can arise from the clustering alone; the
## permutation gives the value expected from clustering profiles that carry
## no genotype-by-location association.
## The same permutations also give the number of markers whose profile values
## differ between the two clustered groups (Kruskal-Wallis, Benjamini-Hochberg
## adjusted P < 0.05), for comparison with the observed counts.
## ENVCLASS_N_PERM (default 200). Writes out/Table_TypeProfile_Null.csv and
## out/Table_DiscriminatingMarkers_Null.csv.
## ============================================================================

suppressPackageStartupMessages({ library(data.table); library(cluster) })
N_PERM <- as.integer(Sys.getenv("ENVCLASS_N_PERM", "200"))
PERM_MODE <- "classification"; MIN_N <- 30; SEED <- 42
MGS <- c("I", "II", "IV")
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


cluster2 <- function(P) {
  mv <- apply(P, 2, var, na.rm = TRUE); keep <- mv > 0
  W <- sweep(P[, keep, drop = FALSE], 2, sqrt(mv[keep]), "*")
  cm <- cor(t(W), use = "pairwise.complete.obs"); cm[is.na(cm)] <- 0
  cutree(hclust(as.dist(1 - cm), method = "ward.D2"), 2)
}
nsig <- function(P, cl) {
  R <- apply(P, 2, rank); N <- nrow(P); n1 <- sum(cl == 1); n2 <- N - n1
  R1 <- colSums(R[cl == 1, , drop = FALSE]); R2 <- N * (N + 1) / 2 - R1
  H <- 12 / (N * (N + 1)) * (R1^2 / n1 + R2^2 / n2) - 3 * (N + 1)
  sum(p.adjust(pchisq(H, 1, lower.tail = FALSE), "BH") < 0.05)
}
type_r <- function(P, cl) cor(colMeans(P[cl == 1, , drop = FALSE]), colMeans(P[cl == 2, , drop = FALSE]))

set.seed(SEED)
out <- list(); out_n <- list()
for (mg_val in MGS) {
  r <- res[[mg_val]]$F; locs <- names(r$clusters)
  dt <- value_table(mg_val)
  blocks <- lapply(setNames(locs, locs), function(l) location_block(dt, l))
  blocks <- blocks[!vapply(blocks, is.null, logical(1))]
  P0 <- do.call(rbind, lapply(blocks, function(b) profile_from(b, b$rc)))
  obs_pub <- type_r(r$profiles, r$clusters[rownames(r$profiles)])
  obs_reb <- type_r(P0, cluster2(P0))
  sims <- replicate(N_PERM, {
    Pp <- do.call(rbind, lapply(blocks, function(b) profile_from(b, sample(b$rc))))
    clp <- cluster2(Pp); c(type_r(Pp, clp), nsig(Pp, clp))
  })
  null <- sims[1, ]; nnull <- sims[2, ]
  obs_n <- nsig(P0, cluster2(P0))
  out_n[[mg_val]] <- data.table(mg = mg_val, observed = obs_n, null_median = median(nnull),
                                null_q95 = quantile(nnull, 0.95), null_max = max(nnull),
                                p = (1 + sum(nnull >= obs_n)) / (N_PERM + 1), n_perm = N_PERM)
  out[[mg_val]] <- data.table(mg = mg_val, n_loc = nrow(P0), r_observed = obs_pub, r_rebuilt = obs_reb,
                              null_mean = mean(null), null_q05 = quantile(null, 0.05), null_q95 = quantile(null, 0.95),
                              p_lower = (1 + sum(null <= obs_reb)) / (N_PERM + 1), n_perm = N_PERM)
  cat(mg_val, "done\n")
}
res_tab <- rbindlist(out); print(res_tab)
fwrite(res_tab, file.path(OUT, "Table_TypeProfile_Null.csv"))
null_n <- rbindlist(out_n); print(null_n)
fwrite(null_n, file.path(OUT, "Table_DiscriminatingMarkers_Null.csv"))
