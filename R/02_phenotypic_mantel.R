if (!exists('ROOT')) source('R/00_setup.R')

## ============================================================================
## EnvClass_PhenoMantel_v3_6.R
## v3.5 + detach ecodist (and proxy, sf, tigris) at the very start of the
## script, in case a previous source() left them attached.
##
## This is the fix for the persistent "logical subscript too long" error:
## ecodist::mantel was masking vegan::mantel via the package: search path,
## and `package:ecodist` was surviving across source() calls.
## ============================================================================

## ---- DETACH FIRST: clean any leftover attachments from previous sources ---
detach_safely <- function(pkg) {
  ns <- paste0("package:", pkg)
  while (ns %in% search()) {
    try(detach(ns, character.only = TRUE, unload = TRUE, force = TRUE),
        silent = TRUE)
  }
}
for (pkg in c("ecodist", "proxy", "sf", "tigris")) detach_safely(pkg)

cat("Search path at start (should NOT contain ecodist, proxy, sf, tigris):\n")
print(search())


library(data.table)
library(vegan)
library(ggplot2)

set.seed(42)

PATHS <- list(
  marker_v2_results = file.path(DATA, "marker_v2_results.rds"),
  nust_pheno        = file.path(DATA, "nust_ecw_subset.csv"),
  nass_cache        = CACHE,
  table_dir         = OUT,
  figure_dir        = FIG
)

PARAMS <- list(
  mgs              = c("0", "I", "II", "III", "IV"),
  n_perm           = 9999,
  min_pair_strains = 8,
  min_loc_coverage = 50
)

cat("\n============================================================================\n")
cat("  PHENOTYPIC MANTEL v3.6\n")
cat("============================================================================\n")


clean_dist <- function(d) {
  m <- as.matrix(d)
  storage.mode(m) <- "double"
  d2 <- stats::as.dist(m)
  class(d2) <- "dist"
  d2
}


standardize_fips <- function(dt, name) {
  dt <- as.data.table(dt)
  candidates <- c("fips", "FIPS", "fips5", "FIPS5", "fips_code",
                  "GEOID", "geoid", "county_fips", "countyfips",
                  "STCOFIPS", "stcofips")
  found <- intersect(candidates, names(dt))
  if (length(found) == 0) {
    pattern_hits <- grep("(?i)fips|geoid", names(dt), value = TRUE)
    if (length(pattern_hits) > 0) found <- pattern_hits[1]
    else stop(sprintf("Cannot locate FIPS in %s. Cols: %s", name,
                      paste(head(names(dt), 10), collapse = ", ")))
  }
  cat(sprintf("    %s: using column '%s'\n", name, found[1]))
  if (found[1] != "fips") setnames(dt, found[1], "fips")
  fc <- dt$fips
  if (is.numeric(fc) || (is.character(fc) && all(grepl("^[0-9]+$", fc[!is.na(fc)])))) {
    n_int <- suppressWarnings(as.integer(fc))
    dt[, fips := ifelse(is.na(n_int), NA_character_, sprintf("%05d", n_int))]
  } else {
    dt[, fips := as.character(fc)]
  }
  dup <- duplicated(dt$fips) & !is.na(dt$fips)
  if (any(dup)) dt <- dt[!dup]
  dt
}


cat("\n[1] Loading core data...\n")
v2   <- readRDS(PATHS$marker_v2_results)
nust <- fread(PATHS$nust_pheno)
nust[, loc_key := trimws(tolower(gsub(",.*$", "", Location)))]

power_data <- standardize_fips(readRDS(file.path(PATHS$nass_cache, "nasa_power_counties.rds")), "NASA POWER")
wx_eng     <- standardize_fips(readRDS(file.path(PATHS$nass_cache, "weather_engineered.rds")), "Engineered weather")
soil_data  <- standardize_fips(readRDS(file.path(PATHS$nass_cache, "ssurgo_multidepth.rds")), "SSURGO soil")
mgmt_data  <- standardize_fips(readRDS(file.path(PATHS$nass_cache, "nass_management.rds")), "NASS management")


cat("\n[2] Mapping locations to FIPS (sf attached temporarily)...\n")
loc_coords <- nust[, .(lat = median(Lat, na.rm = TRUE),
                       lon = median(Long, na.rm = TRUE)), by = loc_key]
loc_coords <- loc_coords[!is.na(lat) & !is.na(lon)]
loc_coords <- unique(loc_coords, by = "loc_key")

suppressPackageStartupMessages({
  library(sf)
  library(tigris)
})
options(tigris_use_cache = TRUE)

us_counties <- counties(cb = TRUE, year = 2022) |> st_transform(4326)
loc_in_us   <- loc_coords[lat > 24 & lat < 50 & lon > -130 & lon < -60]
loc_sf      <- st_as_sf(loc_in_us, coords = c("lon", "lat"), crs = 4326)
loc_join    <- st_join(loc_sf, us_counties[, "GEOID"])
loc_fips    <- as.data.table(loc_join)
loc_fips[, geometry := NULL]
setnames(loc_fips, "GEOID", "fips")
loc_fips[, fips := as.character(fips)]
loc_fips <- loc_fips[!is.na(fips)]
loc_fips <- unique(loc_fips, by = "loc_key")
cat(sprintf("  %d / %d locations mapped\n", nrow(loc_fips), nrow(loc_coords)))

## Detach again now that the geo-join is done.
for (pkg in c("tigris", "sf", "proxy", "ecodist")) detach_safely(pkg)

cat("\n  Search path after detach (should NOT contain ecodist, proxy, sf):\n")
print(search())

if ("package:ecodist" %in% search() || "package:proxy" %in% search()) {
  stop("ecodist or proxy still attached; refusing to continue. ",
       "Restart R and re-run.")
}


cat("\n[3] Building feature stack...\n")
## weather_engineered repeats VPD_mean and VPD_sd from the NASA POWER table;
## keep one copy of each (as R/12 does) so no covariate is counted twice
wx_cols  <- setdiff(names(wx_eng), names(power_data))
features <- merge(power_data, wx_eng[, c("fips", wx_cols), with = FALSE], by = "fips", all = TRUE)
features <- merge(features,   soil_data, by = "fips", all = TRUE)
features <- merge(features,   mgmt_data, by = "fips", all = TRUE)
features <- unique(features, by = "fips")

loc_features <- merge(loc_fips[, .(loc_key, fips)], features, by = "fips", all.x = TRUE)
loc_features <- unique(loc_features, by = "loc_key")

feat_cols <- setdiff(names(loc_features), c("loc_key", "fips"))
keep <- vapply(loc_features[, ..feat_cols], function(x)
  is.numeric(x) && sum(is.finite(x)) >= PARAMS$min_loc_coverage, logical(1))
feat_cols <- feat_cols[keep]

feat_mat <- as.matrix(loc_features[, ..feat_cols])
rownames(feat_mat) <- as.character(loc_features$loc_key)
if (anyDuplicated(rownames(feat_mat))) stop("feat_mat dup rownames")

for (j in seq_len(ncol(feat_mat))) {
  bad <- !is.finite(feat_mat[, j])
  if (any(bad)) feat_mat[bad, j] <- median(feat_mat[, j], na.rm = TRUE)
}
zero_var <- apply(feat_mat, 2, function(x) sd(x, na.rm = TRUE) == 0)
feat_mat <- feat_mat[, !zero_var, drop = FALSE]
feat_mat <- scale(feat_mat)
cat(sprintf("  feat_mat: %d locations x %d features\n",
            nrow(feat_mat), ncol(feat_mat)))

coord_df <- as.data.frame(loc_coords)
rownames(coord_df) <- coord_df$loc_key


compute_adj_dev <- function(dat) {
  em <- dat[, .(em = mean(YieldBuA, na.rm = TRUE)), by = Environment]
  sm <- dat[, .(sm = mean(YieldBuA, na.rm = TRUE)), by = Strain]
  gm <- mean(dat$YieldBuA, na.rm = TRUE)
  dat <- merge(dat, em, by = "Environment")
  dat <- merge(dat, sm, by = "Strain")
  dat[, adj_dev := YieldBuA - sm - em + gm]
  dat
}

build_strain_loc_mat <- function(dt) {
  w <- dcast(dt, Strain ~ loc_key, value.var = "value", fun.aggregate = mean)
  m <- as.matrix(w[, -1]); rownames(m) <- w$Strain
  m
}

pairwise_cor_dist <- function(m, min_pair) {
  K <- ncol(m); cn <- colnames(m)
  d <- matrix(NA_real_, K, K, dimnames = list(cn, cn))
  for (i in seq_len(K - 1)) for (j in (i + 1):K) {
    a <- m[, i]; b <- m[, j]; keep <- !is.na(a) & !is.na(b)
    if (sum(keep) >= min_pair && sd(a[keep]) > 0 && sd(b[keep]) > 0) {
      d[i, j] <- d[j, i] <- 1 - cor(a[keep], b[keep])
    }
  }
  diag(d) <- 0
  d
}


cat("\n[4] Per-MG analysis...\n")

results <- list()

run_one_mg <- function(mg_val) {
  cat(sprintf("\n========= MG %s =========\n", mg_val))

  res <- v2[[mg_val]]$F
  if (is.null(res)) { cat("  No marker-level result; skipping\n"); return(NULL) }

  prof <- res$profiles
  if (!is.matrix(prof)) prof <- as.matrix(prof)
  rownames(prof) <- as.character(rownames(prof))

  mv <- apply(prof, 2, var, na.rm = TRUE)
  mv[is.na(mv)] <- 0
  prof <- prof[, mv > 0, drop = FALSE]

  w    <- sqrt(apply(prof, 2, var, na.rm = TRUE))
  pf_w <- sweep(prof, 2, w, "*")
  cm_g <- cor(t(pf_w), use = "pairwise.complete.obs")
  cm_g[is.na(cm_g)] <- 0
  d_gen_full <- 1 - cm_g

  mg_dat <- nust[MG == mg_val & !is.na(YieldBuA)]
  mg_means <- mg_dat[, .(YieldBuA = mean(YieldBuA, na.rm = TRUE)),
                     by = .(Strain, Environment, loc_key, Year)]
  mg_means <- compute_adj_dev(mg_means)

  loc_in_prof <- intersect(unique(mg_means$loc_key), rownames(prof))
  m_adj <- build_strain_loc_mat(
    mg_means[loc_key %in% loc_in_prof, .(Strain, loc_key, value = adj_dev)])
  d_pheno_full <- pairwise_cor_dist(m_adj, PARAMS$min_pair_strains)

  shared <- Reduce(intersect, list(
    rownames(prof),
    rownames(d_pheno_full),
    rownames(feat_mat),
    rownames(coord_df)
  ))
  shared <- as.character(shared)
  if (length(shared) < 8) { cat("  too few; skip\n"); return(NULL) }

  m_p <- d_pheno_full[shared, shared, drop = FALSE]
  while (any(is.na(m_p)) && nrow(m_p) > 5) {
    na_counts <- rowSums(is.na(m_p))
    worst <- which.max(na_counts)
    m_p <- m_p[-worst, -worst, drop = FALSE]
  }
  shared <- as.character(rownames(m_p))
  cat(sprintf("  shared after prune: %d\n", length(shared)))
  if (length(shared) < 8) { cat("  too few; skip\n"); return(NULL) }

  d_gen   <- clean_dist(stats::as.dist(d_gen_full[shared, shared, drop = FALSE]))
  d_pheno <- clean_dist(stats::as.dist(m_p))
  d_env   <- clean_dist(stats::dist(feat_mat[shared, , drop = FALSE]))
  d_geo   <- clean_dist(stats::dist(as.matrix(coord_df[shared, c("lat","lon"), drop = FALSE])))

  ## Sanity: which mantel function does R see right now?
  resolved <- find("mantel")
  cat(sprintf("  find('mantel'): %s\n", paste(resolved, collapse = ", ")))

  mt_gen <- vegan::mantel(d_gen, d_pheno, permutations = PARAMS$n_perm)
  mt_env <- vegan::mantel(d_env, d_pheno, permutations = PARAMS$n_perm)
  mt_geo <- vegan::mantel(d_geo, d_pheno, permutations = PARAMS$n_perm)

  pm_gen_geo <- vegan::mantel.partial(d_gen, d_pheno, d_geo, permutations = PARAMS$n_perm)
  pm_gen_env <- vegan::mantel.partial(d_gen, d_pheno, d_env, permutations = PARAMS$n_perm)
  pm_env_gen <- vegan::mantel.partial(d_env, d_pheno, d_gen, permutations = PARAMS$n_perm)
  pm_env_geo <- vegan::mantel.partial(d_env, d_pheno, d_geo, permutations = PARAMS$n_perm)

  mrm_fit <- ecodist::MRM(d_pheno ~ d_gen + d_env + d_geo, nperm = PARAMS$n_perm)

  cat(sprintf("  n: %d  | gen r=%.3f (p=%.4f)  env r=%.3f (p=%.4f)  geo r=%.3f (p=%.4f)\n",
              length(shared),
              mt_gen$statistic, mt_gen$signif,
              mt_env$statistic, mt_env$signif,
              mt_geo$statistic, mt_geo$signif))
  cat(sprintf("  partial gen|env=%.3f (p=%.4f), env|gen=%.3f (p=%.4f), gen|geo=%.3f (p=%.4f)\n",
              pm_gen_env$statistic, pm_gen_env$signif,
              pm_env_gen$statistic, pm_env_gen$signif,
              pm_gen_geo$statistic, pm_gen_geo$signif))
  cat(sprintf("  ratio gen/env: %.1fx\n",
              if (mt_env$statistic > 0) mt_gen$statistic / mt_env$statistic else NA))

  data.table(
    mg = mg_val, n_locs = length(shared),
    mantel_genomic = mt_gen$statistic, mantel_genomic_p = mt_gen$signif,
    mantel_enviro  = mt_env$statistic, mantel_enviro_p  = mt_env$signif,
    mantel_geo     = mt_geo$statistic, mantel_geo_p     = mt_geo$signif,
    partial_gen_geo = pm_gen_geo$statistic, partial_gen_geo_p = pm_gen_geo$signif,
    partial_gen_env = pm_gen_env$statistic, partial_gen_env_p = pm_gen_env$signif,
    partial_env_gen = pm_env_gen$statistic, partial_env_gen_p = pm_env_gen$signif,
    partial_env_geo = pm_env_geo$statistic, partial_env_geo_p = pm_env_geo$signif,
    mrm_r2 = mrm_fit$r.squared[1],
    mrm_beta_gen = mrm_fit$coef[2, 1], mrm_p_gen = mrm_fit$coef[2, 2],
    mrm_beta_env = mrm_fit$coef[3, 1], mrm_p_env = mrm_fit$coef[3, 2],
    mrm_beta_geo = mrm_fit$coef[4, 1], mrm_p_geo = mrm_fit$coef[4, 2],
    ratio = if (mt_env$statistic > 0) mt_gen$statistic / mt_env$statistic else NA_real_
  )
}


for (mg_val in PARAMS$mgs) {
  res_dt <- tryCatch(run_one_mg(mg_val), error = function(e) {
    cat(sprintf("  >>> ERROR in MG %s: %s\n", mg_val, conditionMessage(e)))
    NULL
  })
  if (!is.null(res_dt)) results[[mg_val]] <- res_dt
}

if (length(results) == 0) stop("No MG produced results.")

out <- rbindlist(results, fill = TRUE)
fwrite(out, file.path(PATHS$table_dir, "Table_PhenoMantel_Final.csv"))


plot_dt <- rbindlist(lapply(results, function(r) {
  data.table(mg = r$mg,
             measure = c("Genomic", "Envirotyping", "Geographic"),
             r = c(r$mantel_genomic, r$mantel_enviro, r$mantel_geo),
             p = c(r$mantel_genomic_p, r$mantel_enviro_p, r$mantel_geo_p))
}))
plot_dt[, measure := factor(measure,
                            levels = c("Genomic", "Envirotyping", "Geographic"))]
plot_dt[, sig := fcase(p < 0.001, "***",
                       p < 0.01,  "**",
                       p < 0.05,  "*",
                       default    = "ns")]
plot_dt[, mg := factor(mg, levels = rev(c("0", "I", "II", "III", "IV")))]

p_forest <- ggplot(plot_dt, aes(x = r, y = mg, color = measure)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
  geom_point(size = 4, position = position_dodge(width = 0.55)) +
  geom_text(aes(label = sig), vjust = -1.0, size = 3.2, show.legend = FALSE,
            position = position_dodge(width = 0.55)) +
  scale_color_manual(values = c("Genomic"      = "#D62728",
                                "Envirotyping" = "#1F77B4",
                                "Geographic"   = "#2CA02C")) +
  labs(x = "Mantel r vs phenotypic response distance",
       y = "Maturity group", color = NULL) +
  theme_minimal(base_size = 13) + theme(legend.position = "top")

ggsave(file.path(PATHS$figure_dir, "Fig_PhenoMantel_Final.png"),
       p_forest, width = 8, height = 5, dpi = 300)

ratio_dt <- out[, .(mg, ratio)][!is.na(ratio)]
if (nrow(ratio_dt) > 0) {
  ratio_dt[, mg := factor(mg, levels = c("0", "I", "II", "III", "IV"))]
  p_ratio <- ggplot(ratio_dt, aes(x = mg, y = ratio)) +
    geom_col(fill = "#D62728", alpha = 0.85) +
    geom_text(aes(label = sprintf("%.1fx", ratio)),
              vjust = -0.4, size = 4.5) +
    geom_hline(yintercept = 1, linetype = "dashed", color = "grey50") +
    labs(x = "Maturity group",
         y = "Genomic / envirotyping Mantel r ratio") +
    theme_minimal(base_size = 13) +
    expand_limits(y = max(ratio_dt$ratio, na.rm = TRUE) * 1.15)
  ggsave(file.path(PATHS$figure_dir, "Fig_PhenoMantel_Final_Ratio.png"),
         p_ratio, width = 7, height = 5, dpi = 300)
}


cat("\n[5] Headline summary:\n\n")
print(out[, .(mg, n_locs,
              gen_r = round(mantel_genomic, 3), gen_p = mantel_genomic_p,
              env_r = round(mantel_enviro, 3),  env_p = mantel_enviro_p,
              gen_ctrl_env = round(partial_gen_env, 3),
              env_ctrl_gen = round(partial_env_gen, 3),
              ratio = round(ratio, 1))])


cat("\n============================================================================\n")
cat("  COMPLETE\n")
cat("============================================================================\n")
