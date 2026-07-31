if (!exists('ROOT')) source('R/00_setup.R')

## ============================================================================
## EnvClass_MethodGrid_v2.R
##
## Table_Method_Grid.csv scores 15 clustering methods against a reference that
## is not the published partition: correlation|ward.D2|variance, which the
## Methods describe as the published method, does not score ARI = 1, and
## euclidean|ward.D2|variance does. The grid k also differs from published k
## in MG I and MG III.
##
## This rebuilds the same 15 combinations against v2[[mg]]$F$clusters, the
## partition the manuscript actually reports.
##
## Two ARIs per method:
##   ari_at_published_k  cut at the published k, so it measures partition
##                       agreement alone
##   ari_at_own_k        each method choosing its own k by silhouette, which
##                       is the harder test
##
## Output: Table_Method_Grid_v2.csv
## ============================================================================

library(data.table)
library(cluster)
library(mclust)

set.seed(42)

PATHS <- list(
  marker_v2_results = file.path(DATA, "marker_v2_results.rds"),
  table_dir         = OUT
)

PARAMS <- list(mgs = c("0", "I", "II", "III", "IV"), max_k = 8)

## the same 15 combinations as the original table
GRID <- rbindlist(list(
  data.table(metric = "correlation", link = c("ward.D2","complete","average","pam"), wt = "variance"),
  data.table(metric = "euclidean",   link = c("ward.D2","complete","average","pam"), wt = "variance"),
  data.table(metric = "manhattan",   link = c("ward.D2","complete","pam"),           wt = "unweighted"),
  data.table(metric = "cosine",      link = c("ward.D2","complete","average","pam"), wt = "variance")
))

v2 <- readRDS(PATHS$marker_v2_results)

weight_profiles <- function(pf, wt) {
  mv <- apply(pf, 2, var, na.rm = TRUE)
  pf <- pf[, mv > 0, drop = FALSE]
  if (wt == "unweighted") return(pf)
  sweep(pf, 2, sqrt(apply(pf, 2, var, na.rm = TRUE)), "*")
}

make_dist <- function(m, metric) {
  if (metric == "correlation") {
    cm <- cor(t(m), use = "pairwise.complete.obs"); cm[is.na(cm)] <- 0
    return(as.dist(1 - cm))
  }
  if (metric == "cosine") {
    nrm <- sqrt(rowSums(m^2)); nrm[nrm == 0] <- 1
    mn <- m / nrm
    return(as.dist(1 - tcrossprod(mn)))
  }
  dist(m, method = metric)
}

cut_at <- function(d, link, k) {
  if (link == "pam") pam(d, k)$clustering else cutree(hclust(d, method = link), k)
}

best_k <- function(d, link) {
  n <- attr(d, "Size")
  ks <- 2:min(PARAMS$max_k, n - 1)
  s <- sapply(ks, function(k) mean(silhouette(cut_at(d, link, k), d)[, 3]))
  ks[which.max(s)]
}

out <- list()
for (mg in PARAMS$mgs) {
  res <- v2[[mg]]$F
  if (is.null(res) || is.null(res$clusters)) next
  published <- res$clusters
  k_pub <- length(unique(published))
  locs <- names(published)

  for (i in seq_len(nrow(GRID))) {
    metric <- GRID$metric[i]; link <- GRID$link[i]; wt <- GRID$wt[i]
    m <- weight_profiles(res$profiles, wt)[locs, , drop = FALSE]
    d <- make_dist(m, metric)

    cl_pub_k <- cut_at(d, link, k_pub)
    k_own    <- best_k(d, link)
    cl_own_k <- cut_at(d, link, k_own)

    out[[length(out) + 1]] <- data.table(
      mg = mg, method_id = paste(metric, link, wt, sep = "|"),
      k_published = k_pub, k_own = k_own,
      ari_at_published_k = adjustedRandIndex(published, cl_pub_k),
      ari_at_own_k       = adjustedRandIndex(published, cl_own_k),
      sil_at_published_k = mean(silhouette(cl_pub_k, d)[, 3]),
      n_locs = length(locs))
  }
}
out <- rbindlist(out)
fwrite(out, file.path(PATHS$table_dir, "Table_Method_Grid_v2.csv"))

## ---- sanity check: the published method must score 1 -----------------------
cat("\n--- published method against itself (must be ARI = 1) ---\n")
print(out[method_id == "correlation|ward.D2|variance",
          .(mg, k_published, k_own, ari_at_published_k, ari_at_own_k)])

cat("\n--- agreement with the published partition, by MG ---\n")
print(out[, .(n_methods = .N,
              med_at_k  = round(median(ari_at_published_k), 2),
              min_at_k  = round(min(ari_at_published_k), 2),
              max_at_k  = round(max(ari_at_published_k), 2),
              med_own_k = round(median(ari_at_own_k), 2)), by = mg])

cat("\n--- by distance metric, pooled over MGs ---\n")
out[, metric := sub("\\|.*$", "", method_id)]
print(out[, .(med_at_k = round(median(ari_at_published_k), 2)), by = metric])

cat("\nIf correlation|ward.D2|variance does not return 1, the profiles or the\n")
cat("cluster object are not the ones that produced the published partition.\n")
