if (!exists('ROOT')) source('R/00_setup.R')

## ============================================================================
## DISJOINT-GENOTYPE MANTEL (circularity control), all maturity groups.
## Self-contained: just run this file. Needs data.table, vegan, ggplot2.
## MRP built from genotype half A; phenotypic response distance from half B;
## the two halves share no genotypes. Locations with a constant (degenerate)
## half-sample profile are dropped. Prints a summary, writes a TIFF and a CSV.
## ============================================================================
suppressPackageStartupMessages({ library(data.table); library(vegan); library(ggplot2) })
set.seed(42)

CFG <- list(
  results_rds = file.path(DATA, "marker_v2_results.rds"),
  nust_file   = file.path(DATA, "nust_ecw_subset.csv"),
  geno_worker = file.path(DATA, "geno_worker_data.rds"),
  table_dir   = OUT,
  out_dir     = FIG,
  yield_col   = "YieldBuA",
  n_split     = 30,
  min_n_mrp   = 10,   # strains needed to build a half-A profile
  min_shared  = 6,    # half-B strains needed per location pair
  min_locs    = 6
)
dir.create(CFG$out_dir, showWarnings = FALSE, recursive = TRUE)

load_raw <- function() {
  nust <- fread(CFG$nust_file)
  nust[, loc_key := trimws(tolower(gsub(",.*$", "", Location)))]
  gw <- readRDS(CFG$geno_worker); gdf <- gw$geno_df
  strain_names <- names(gdf)[-(1:3)]
  geno_mat <- as.matrix(gdf[, -(1:3)]); rownames(geno_mat) <- as.character(gdf$marker)
  list(nust = nust, geno_mat = geno_mat, strain_names = strain_names,
       marker_ids = as.character(gdf$marker))
}

mrp_profiles <- function(dt, locs, RAW, min_n = CFG$min_n_mrp) {
  P <- matrix(NA_real_, length(locs), length(RAW$marker_ids),
              dimnames = list(locs, RAW$marker_ids))
  for (i in seq_along(locs)) {
    ld <- dt[loc_key == locs[i], .(val = mean(value, na.rm = TRUE)), by = Strain]
    sh <- intersect(ld$Strain, RAW$strain_names)
    if (length(sh) < min_n) next
    rv <- ld[match(sh, Strain), val]
    sc <- match(sh, RAW$strain_names)
    gs <- RAW$geno_mat[, sc, drop = FALSE]
    gsd <- apply(gs, 1, sd, na.rm = TRUE)
    rc  <- rv - mean(rv, na.rm = TRUE)
    gc  <- sweep(gs, 1, rowMeans(gs, na.rm = TRUE))
    cp  <- gc %*% rc / (length(rv) - 1)
    cors <- as.numeric(cp) / (gsd * sd(rv, na.rm = TRUE))
    cors[gsd == 0 | !is.finite(cors)] <- 0
    P[i, ] <- cors
  }
  P[!is.na(P[, 1]), , drop = FALSE]
}

RAW <- load_raw()
res <- readRDS(CFG$results_rds)
dt  <- copy(RAW$nust)
if (!"Environment" %in% names(dt)) dt[, Environment := paste(loc_key, Year, sep = "_")]
dt[, y := get(CFG$yield_col)]
gm <- mean(dt$y, na.rm = TRUE)
dt[, sm := ave(y, Strain,      FUN = function(z) mean(z, na.rm = TRUE))]
dt[, em := ave(y, Environment, FUN = function(z) mean(z, na.rm = TRUE))]
dt[, adj := y - sm - em + gm]
dt[, value := adj]

one_split <- function(locs) {
  half <- sample(RAW$strain_names, length(RAW$strain_names) %/% 2)
  A <- setdiff(RAW$strain_names, half); B <- half
  Pa <- mrp_profiles(dt[Strain %in% A], locs, RAW)
  if (nrow(Pa) < CFG$min_locs) return(c(r = NA_real_, n = NA_real_))
  sub <- dt[loc_key %in% rownames(Pa) & Strain %in% B,
            .(a = mean(adj, na.rm = TRUE)), by = .(Strain, loc_key)]
  w  <- dcast(sub, Strain ~ loc_key, value.var = "a")
  am <- as.matrix(w[, -1]); rownames(am) <- w$Strain
  L  <- colnames(am)
  M  <- matrix(NA_real_, length(L), length(L), dimnames = list(L, L))
  for (i in seq_along(L)) for (j in i:length(L)) {
    a <- am[, i]; b <- am[, j]; ok <- is.finite(a) & is.finite(b)
    if (sum(ok) >= CFG$min_shared && sd(a[ok]) > 0 && sd(b[ok]) > 0)
      M[i, j] <- M[j, i] <- 1 - cor(a[ok], b[ok])
  }
  diag(M) <- 0
  keep <- intersect(rownames(Pa), rownames(M))
  Ms <- M[keep, keep, drop = FALSE]
  while (nrow(Ms) >= CFG$min_locs && any(is.na(Ms)))
    Ms <- Ms[-which.max(rowSums(is.na(Ms))), -which.max(rowSums(is.na(Ms))), drop = FALSE]
  keep2 <- rownames(Ms)[rowSums(is.na(Ms)) == 0]
  if (length(keep2) < CFG$min_locs) return(c(r = NA_real_, n = NA_real_))
  Pk <- Pa[keep2, , drop = FALSE]
  mv <- apply(Pk, 2, var, na.rm = TRUE); Pk <- Pk[, mv > 0, drop = FALSE]
  Pw <- sweep(Pk, 2, sqrt(apply(Pk, 2, var, na.rm = TRUE)), "*")
  psd <- apply(Pw, 1, sd, na.rm = TRUE)
  good <- rownames(Pw)[is.finite(psd) & psd > 0]
  if (length(good) < CFG$min_locs) return(c(r = NA_real_, n = NA_real_))
  gd <- as.dist(1 - cor(t(Pw[good, , drop = FALSE])))
  pd <- as.dist(M[good, good])
  st <- tryCatch(vegan::mantel(gd, pd, permutations = 0)$statistic, error = function(e) NA_real_)
  c(r = st, n = length(good))
}

mgs <- intersect(c("I","II","III","IV"), names(res))
rows <- list(); split_long <- list()
for (mg in mgs) {
  r <- res[[mg]]$F; if (is.null(r)) next
  locs <- names(r$clusters)
  mat <- sapply(seq_len(CFG$n_split), function(i) one_split(locs))   # 2 x n_split
  rr <- mat["r", ]; nn <- mat["n", ]
  ok <- is.finite(rr); rr <- rr[ok]; nn <- nn[ok]
  rows[[mg]] <- data.frame(
    mg = mg, n_valid = length(rr),
    median_locs = if (length(nn)) as.integer(median(nn)) else NA,
    disjoint_mean = if (length(rr)) round(mean(rr), 3) else NA,
    disjoint_median = if (length(rr)) round(median(rr), 3) else NA,
    q05 = if (length(rr)) round(quantile(rr, .05), 3) else NA,
    q95 = if (length(rr)) round(quantile(rr, .95), 3) else NA)
  if (length(rr)) split_long[[mg]] <- data.table(MG = paste0("MG ", mg), split_r = rr)
}

summary_tab <- do.call(rbind, rows)
fin <- fread(file.path(CFG$table_dir, "Table_PhenoMantel_Final.csv"))
summary_tab <- merge(summary_tab, fin[, .(mg, full_sample = round(mantel_genomic, 3))],
                     by = "mg", all.x = TRUE)
cat("\n=== Disjoint-genotype Mantel (MRP half A vs phenotype half B) ===\n")
print(summary_tab, row.names = FALSE)
fwrite(summary_tab, file.path(CFG$table_dir, "Table_DisjointMantel.csv"))
cat("\nwrote", file.path(CFG$table_dir, "Table_DisjointMantel.csv"), "\n")

if (length(split_long)) {
  sl <- rbindlist(split_long)
  sl[, MG := factor(MG, levels = paste0("MG ", c("0","I","II","III","IV")))]
  fp <- fin[, .(MG = factor(paste0("MG ", mg), levels = levels(sl$MG)), full = mantel_genomic)]
  fp <- fp[MG %in% sl$MG]
  p <- ggplot(sl, aes(MG, split_r)) +
    geom_hline(yintercept = 0, linewidth = 0.4, colour = "grey60") +
    geom_boxplot(fill = "#cfd4da", colour = "grey30", width = 0.6, outlier.size = 0.5) +
    geom_point(data = fp, aes(MG, full), colour = "#2166AC", size = 2.4) +
    labs(x = NULL, y = "Mantel r (profile from genotype half A,\nphenotype from non-overlapping half B)") +
    theme_classic(base_size = 9, base_family = "sans") +
    theme(axis.text  = element_text(colour = "black"),
          axis.line  = element_line(linewidth = 0.4),
          axis.ticks = element_line(linewidth = 0.4, colour = "black"))
  ggsave(file.path(CFG$out_dir, "Fig6_Disjoint_Genotype_Mantel.tiff"), p,
         width = 129, height = 95, units = "mm", dpi = 300, compression = "lzw", bg = "white")
  cat("wrote", file.path(CFG$out_dir, "Fig6_Disjoint_Genotype_Mantel.tiff"), "\n")
} else {
  message("No valid splits in any MG; nothing written.")
}
