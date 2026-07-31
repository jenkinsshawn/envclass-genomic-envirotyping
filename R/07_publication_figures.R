if (!exists('ROOT')) source('R/00_setup.R')

## ============================================================================
## PUBLICATION FIGURES: Genomic envirotyping (NUST marker response profiles)
## Produces 300-dpi LZW TIFFs sized to Springer / TAG specifications.
## ----------------------------------------------------------------------------
## The script is in three parts:
##   PART A  fast, runs from existing output files (CSVs + marker_v2_results.rds)
##   PART B  county map, needs the 'tigris' package (downloads county shapefiles)
##   PART C  recomputed validation figures, needs the raw data files
##           (NUST_ECW.csv + geno_worker_data.rds) and takes longer
## Each figure is self-contained; comment out parts you do not need.
## ============================================================================

## ---- CONFIG: edit these paths to match your machine -------------------------
CFG <- list(
  table_dir  = OUT,          # where the Table_*.csv live
  results_rds = file.path(DATA, "marker_v2_results.rds"),
  nust_file   = file.path(DATA, "nust_ecw_subset.csv"),
  geno_worker = file.path(DATA, "geno_worker_data.rds"),
  out_dir     = FIG,
  yield_col   = "YieldBuA",   # raw yield column in NUST_ECW.csv (Part C)
  trait_value = "YieldBuA",   # trait used to build marker response profiles
  min_loc_strains = 30,       # matches PARAMS$min_loc_strains
  cluster_method  = "ward.D2",
  max_k = 8,
  seed = 42
)
dir.create(CFG$out_dir, showWarnings = FALSE, recursive = TRUE)
set.seed(CFG$seed)

suppressPackageStartupMessages({
  library(data.table); library(ggplot2)
})
have <- function(p) requireNamespace(p, quietly = TRUE)
if (!have("patchwork")) message("Install 'patchwork' for multi-panel figures: install.packages('patchwork')")

## ---- shared theme + palette -------------------------------------------------
theme_pub <- function(base = 9) {
  theme_classic(base_size = base, base_family = "sans") +
    theme(
      axis.text  = element_text(colour = "black"),
      axis.line  = element_line(linewidth = 0.4),
      axis.ticks = element_line(linewidth = 0.4, colour = "black"),
      plot.title = element_text(face = "bold", size = base + 1, hjust = 0),
      legend.key.size = unit(3.5, "mm"),
      legend.background = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = base)
    )
}
PAL <- list(gen = "#2166AC", env = "#D1772F", geo = "#9AA0A6")
# Okabe-Ito qualitative palette for clusters
PAL_CL <- c("#0072B2","#D55E00","#009E73","#CC79A7","#E69F00","#56B4E9","#F0E442","#000000")

## ---- TIFF saver (uses ragg if available for crisper text) -------------------
save_tiff <- function(plot, name, w_mm, h_mm, dpi = 300) {
  path <- file.path(CFG$out_dir, name)
  if (have("ragg")) {
    ragg::agg_tiff(path, width = w_mm, height = h_mm, units = "mm", res = dpi,
                   compression = "lzw", background = "white")
    print(plot); dev.off()
  } else {
    ggsave(path, plot, device = "tiff", width = w_mm, height = h_mm,
           units = "mm", dpi = dpi, compression = "lzw", bg = "white")
  }
  message("wrote ", path)
}
sig_lab <- function(p) ifelse(p < 0.001, "***", ifelse(p < 0.01, "**",
                                                       ifelse(p < 0.05, "*", "n.s.")))

## ============================================================================
## PART A  -- figures from existing output files
## ============================================================================

## ---- FIGURE 1: Mantel headline (a raw, b partial) ---------------------------
fin <- fread(file.path(CFG$table_dir, "Table_PhenoMantel_Final.csv"))
fin[, MG := factor(paste0("MG ", mg), levels = paste0("MG ", mg))]

# panel a: raw Mantel (genomic / envirotyping / geographic)
a_long <- rbindlist(list(
  data.table(MG = fin$MG, type = "Genomic",     r = fin$mantel_genomic, p = fin$mantel_genomic_p),
  data.table(MG = fin$MG, type = "Envirotyping", r = fin$mantel_enviro,  p = fin$mantel_enviro_p),
  data.table(MG = fin$MG, type = "Geographic",   r = fin$mantel_geo,     p = fin$mantel_geo_p)
))
a_long[, type := factor(type, levels = c("Genomic","Envirotyping","Geographic"))]
a_long[, lab := sig_lab(p)]

pa <- ggplot(a_long, aes(MG, r, fill = type)) +
  geom_col(position = position_dodge(0.8), width = 0.72, colour = "white", linewidth = 0.2) +
  geom_text(aes(label = lab, colour = lab == "n.s."),
            position = position_dodge(0.8), vjust = -0.35, size = 2.6, show.legend = FALSE) +
  scale_fill_manual(values = c(Genomic = PAL$gen, Envirotyping = PAL$env, Geographic = PAL$geo), name = NULL) +
  scale_colour_manual(values = c(`TRUE` = "#666666", `FALSE` = "black")) +
  scale_y_continuous(limits = c(0, 0.72), expand = expansion(mult = c(0, 0.02))) +
  labs(title = "a", y = "Mantel r with phenotypic response distance", x = NULL) +
  theme_pub() + theme(legend.position = "top")

# panel b: partial Mantel asymmetry
b_long <- rbindlist(list(
  data.table(MG = fin$MG, type = "Genomic | envirotyping", r = fin$partial_gen_env, p = fin$partial_gen_env_p),
  data.table(MG = fin$MG, type = "Envirotyping | genomic", r = fin$partial_env_gen, p = fin$partial_env_gen_p)
))
b_long[, type := factor(type, levels = c("Genomic | envirotyping","Envirotyping | genomic"))]
b_long[, lab := sig_lab(p)]

pb <- ggplot(b_long, aes(MG, r, fill = type)) +
  geom_col(position = position_dodge(0.7), width = 0.62, colour = "white", linewidth = 0.2) +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  geom_text(aes(label = lab, colour = lab == "n.s."),
            position = position_dodge(0.7), vjust = -0.4, size = 2.6, show.legend = FALSE) +
  scale_fill_manual(values = c("Genomic | envirotyping" = PAL$gen,
                               "Envirotyping | genomic" = PAL$env), name = NULL) +
  scale_colour_manual(values = c(`TRUE` = "#666666", `FALSE` = "black")) +
  scale_y_continuous(limits = c(-0.1, 0.72)) +
  labs(title = "b", y = "Partial Mantel r", x = NULL) +
  theme_pub() + theme(legend.position = "top")

if (have("patchwork")) {
  fig1 <- patchwork::wrap_plots(pa, pb, nrow = 1)
  save_tiff(fig1, "Fig1_Mantel.tiff", w_mm = 174, h_mm = 88)
} else {
  save_tiff(pa, "Fig1a_Mantel_raw.tiff", 90, 88)
  save_tiff(pb, "Fig1b_Mantel_partial.tiff", 90, 88)
}

## ---- FIGURE 2: per-MG classification (MDS of marker-level profiles) ----------
## Replaces the chromosome-level PCA. Uses the canonical distances + clusters.
res <- readRDS(CFG$results_rds)
mds_dt <- rbindlist(lapply(names(res), function(mg) {
  r <- res[[mg]]$F
  if (is.null(r) || is.null(r$dist)) return(NULL)
  co <- cmdscale(r$dist, k = 2)
  data.table(MG = paste0("MG ", mg),
             MDS1 = co[,1], MDS2 = co[,2],
             cluster = factor(r$clusters[rownames(co)]),
             sil = r$sil, k = r$k)
}), fill = TRUE)
mds_dt[, MG := factor(MG, levels = paste0("MG ", c("0","I","II","III","IV")))]
lab_dt <- unique(mds_dt[, .(MG, k, sil)])
lab_dt[, lab := sprintf("k = %d, silhouette = %.2f", k, sil)]

fig2 <- ggplot(mds_dt, aes(MDS1, MDS2, colour = cluster)) +
  geom_point(size = 1.6, alpha = 0.9) +
  facet_wrap(~ MG, scales = "free", nrow = 2) +
  geom_text(data = lab_dt, aes(x = -Inf, y = Inf, label = lab),
            hjust = -0.05, vjust = 1.4, size = 2.4, inherit.aes = FALSE) +
  scale_colour_manual(values = PAL_CL, name = "Type") +
  labs(title = NULL, x = "MDS axis 1", y = "MDS axis 2") +
  theme_pub() + theme(legend.position = "bottom")
save_tiff(fig2, "Fig2_Classification_MDS.tiff", w_mm = 174, h_mm = 120)

## ---- FIGURE 5: RF variable importance (no-yield projection) ------------------
## Top 15 features per MG, from the 87-feature circularity-clean model.
imp <- rbindlist(lapply(c("0","I","II","III","IV"), function(mg) {
  f <- file.path(CFG$table_dir, sprintf("RF_Importance_NoYield_MG%s.csv", mg))
  if (!file.exists(f)) return(NULL)
  d <- fread(f); setnames(d, 1:2, c("feature","MDA"))
  d <- d[order(-MDA)][1:15]
  d[, MG := paste0("MG ", mg)]; d
}), fill = TRUE)
imp[, MG := factor(MG, levels = paste0("MG ", c("0","I","II","III","IV")))]

if (have("tidytext")) {
  imp[, feat_ord := tidytext::reorder_within(feature, MDA, MG)]
  fig5 <- ggplot(imp, aes(MDA, feat_ord)) +
    geom_col(fill = PAL$gen, width = 0.7) +
    facet_wrap(~ MG, scales = "free_y", nrow = 2) +
    tidytext::scale_y_reordered() +
    labs(title = NULL, x = "Mean decrease in accuracy", y = NULL) +
    theme_pub(8) + theme(axis.text.y = element_text(size = 6))
} else {
  message("Install 'tidytext' for correct per-panel bar ordering; using global order.")
  fig5 <- ggplot(imp, aes(MDA, reorder(feature, MDA))) +
    geom_col(fill = PAL$gen, width = 0.7) +
    facet_wrap(~ MG, scales = "free_y", nrow = 2) +
    labs(title = NULL, x = "Mean decrease in accuracy", y = NULL) +
    theme_pub(8) + theme(axis.text.y = element_text(size = 6))
}
save_tiff(fig5, "Fig5_RF_Importance.tiff", w_mm = 174, h_mm = 140)

## ---- (optional) ADVANTAGE PANEL: genomic/envirotyping ratio by MG -----------
adv <- fin[, .(MG, ratio)]
fig_adv <- ggplot(adv, aes(MG, ratio)) +
  geom_col(fill = PAL$gen, width = 0.6) +
  geom_text(aes(label = sprintf("%.1f", ratio)), vjust = -0.3, size = 2.8) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(title = NULL, y = "Genomic / envirotyping Mantel ratio", x = NULL) +
  theme_pub()
save_tiff(fig_adv, "FigSx_Advantage_ratio.tiff", w_mm = 84, h_mm = 80)

## ============================================================================
## PART B  -- county projection map (needs 'tigris' + 'sf')
## ============================================================================
make_county_map <- function() {
  if (!all(sapply(c("sf","tigris","dplyr"), have))) {
    message("Skipping map: install.packages(c('sf','tigris','dplyr'))"); return(invisible())
  }
  suppressPackageStartupMessages({ library(sf); library(tigris); library(dplyr) })
  options(tigris_use_cache = TRUE)
  
  cty <- tigris::counties(cb = TRUE, resolution = "20m", year = 2021, progress_bar = FALSE)
  cty <- cty[!cty$STATEFP %in% c("02","15","60","66","69","72","78"), ]   # drop AK/HI/territories
  cty$fips <- as.integer(paste0(cty$STATEFP, cty$COUNTYFP))
  ctr <- suppressWarnings(st_coordinates(st_centroid(st_make_valid(cty))))
  cty$lat <- ctr[,2]
  
  # latitude-based MG zone (matches Methods)
  mg_of_lat <- function(lat) ifelse(lat >= 45, "0", ifelse(lat >= 43, "I",
                                                           ifelse(lat >= 40, "II", ifelse(lat >= 37, "III", "IV"))))
  cty$MGzone <- mg_of_lat(cty$lat)
  
  # read no-yield predictions per MG, keep each county's prediction for its zone
  preds <- rbindlist(lapply(c("0","I","II","III","IV"), function(mg) {
    f <- file.path(CFG$table_dir, sprintf("County_Predictions_NoYield_MG%s.csv", mg))
    if (!file.exists(f)) return(NULL)
    d <- fread(f); d[, MGzone := mg]
    d[, env_type := paste0("MG", mg, "_C", predicted)]; d
  }), fill = TRUE)
  
  m <- merge(cty, as.data.frame(preds[, .(fips, MGzone, env_type, confidence)]),
             by = c("fips","MGzone"))
  m <- st_transform(m, 5070)   # CONUS Albers
  
  # NUST testing locations (analysis set) as black reference points
  nust_coords <- fread(CFG$nust_file, select = c("Location", "Lat", "Long"))
  nust_coords[, loc_key := trimws(tolower(gsub(",.*$", "", Location)))]
  loc_coords <- nust_coords[, .(lat = median(Lat, na.rm = TRUE),
                                lon = median(Long, na.rm = TRUE)), by = loc_key]
  used <- tryCatch(
    unique(unlist(lapply(readRDS(CFG$results_rds), function(r)
      if (!is.null(r$F)) names(r$F$clusters)))),
    error = function(e) loc_coords$loc_key)
  pts_sf <- st_transform(
    st_as_sf(loc_coords[loc_key %in% used & is.finite(lat) & is.finite(lon)],
             coords = c("lon", "lat"), crs = 4326), 5070)
  
  pmap <- ggplot(m) +
    geom_sf(aes(fill = env_type, alpha = confidence), colour = NA) +
    geom_sf(data = pts_sf, shape = 21, fill = "black", colour = "white",
            stroke = 0.15, size = 1.2, inherit.aes = FALSE) +
    scale_alpha(range = c(0.45, 1), name = "Confidence") +
    labs(title = NULL, fill = "Environment type") +
    theme_void(base_family = "sans") +
    theme(legend.key.size = unit(3.5, "mm"), legend.text = element_text(size = 6),
          legend.title = element_text(size = 7))
  save_tiff(pmap, "Fig4_County_Projection.tiff", w_mm = 174, h_mm = 120)
}
make_county_map()

## ============================================================================
## PART C  -- recomputed validation figures (needs raw data; slower)
##   Fig 3  marker-level permutation null
##   Fig 6  disjoint-genotype Mantel (circularity control)
## ============================================================================
load_raw <- function() {
  nust <- fread(CFG$nust_file)
  nust[, loc_key := trimws(tolower(gsub(",.*$", "", Location)))]
  gw <- readRDS(CFG$geno_worker); gdf <- gw$geno_df
  strain_names <- names(gdf)[-(1:3)]
  geno_mat <- as.matrix(gdf[, -(1:3)]); rownames(geno_mat) <- as.character(gdf$marker)
  list(nust = nust, geno_mat = geno_mat, strain_names = strain_names,
       marker_ids = as.character(gdf$marker))
}

# marker response profiles for a set of locations (optionally permuting geno-pheno link)
mrp_profiles <- function(dt, locs, RAW, permute = FALSE, min_n = CFG$min_loc_strains) {
  P <- matrix(NA_real_, length(locs), length(RAW$marker_ids), dimnames = list(locs, RAW$marker_ids))
  for (i in seq_along(locs)) {
    ld <- dt[loc_key == locs[i], .(val = mean(value, na.rm = TRUE)), by = Strain]
    sh <- intersect(ld$Strain, RAW$strain_names)
    if (length(sh) < min_n) next
    rv <- ld[match(sh, Strain), val]
    sc <- match(sh, RAW$strain_names)
    if (permute) sc <- sample(sc)                    # break geno-pheno link within location
    gs <- RAW$geno_mat[, sc, drop = FALSE]
    gsd <- apply(gs, 1, sd, na.rm = TRUE)
    rc  <- rv - mean(rv, na.rm = TRUE)
    gc  <- sweep(gs, 1, rowMeans(gs, na.rm = TRUE))
    cp  <- gc %*% rc / (length(rv) - 1)
    cors <- as.numeric(cp) / (gsd * sd(rv, na.rm = TRUE))
    cors[gsd == 0 | !is.finite(cors)] <- 0
    P[i, ] <- cors
  }
  P[!is.na(P[,1]), , drop = FALSE]
}

best_sil <- function(P, max_k = CFG$max_k) {
  if (nrow(P) < 5) return(NA_real_)
  mv <- apply(P, 2, var, na.rm = TRUE); P <- P[, mv > 0, drop = FALSE]
  P <- sweep(P, 2, sqrt(apply(P, 2, var, na.rm = TRUE)), "*")   # variance weighting
  cm <- cor(t(P), use = "pairwise.complete.obs"); cm[is.na(cm)] <- 0
  d <- as.dist(1 - cm); hc <- hclust(d, method = CFG$cluster_method)
  sapply(2:max_k, function(k) mean(cluster::silhouette(cutree(hc, k), d)[, 3])) -> s
  max(s)
}

## ---- FIGURE 3: permutation null vs observed silhouette ----------------------
make_permnull <- function(n_perm = 100) {
  RAW <- load_raw()
  res <- readRDS(CFG$results_rds)
  dt  <- copy(RAW$nust)
  if (!"Environment" %in% names(dt)) dt[, Environment := paste(loc_key, Year, sep = "_")]
  gm <- mean(dt[[CFG$yield_col]], na.rm = TRUE)
  dt[, sm := ave(get(CFG$yield_col), Strain,      FUN = function(z) mean(z, na.rm = TRUE))]
  dt[, em := ave(get(CFG$yield_col), Environment, FUN = function(z) mean(z, na.rm = TRUE))]
  dt[, value := get(CFG$yield_col) - sm - em + gm]   # MRP = adjusted deviations
  out <- rbindlist(lapply(names(res), function(mg) {
    r <- res[[mg]]$F; if (is.null(r)) return(NULL)
    locs <- names(r$clusters)
    nulls <- replicate(n_perm, best_sil(mrp_profiles(dt, locs, RAW, permute = TRUE)))
    data.table(MG = paste0("MG ", mg), observed = r$sil, null = nulls)
  }), fill = TRUE)
  out[, MG := factor(MG, levels = paste0("MG ", c("0","I","II","III","IV")))]
  obs <- unique(out[, .(MG, observed)])
  cat("\n=== permutation null, per MG ===\n")
  print(out[, .(observed  = round(observed[1], 3),
                null_mean = round(mean(null, na.rm = TRUE), 3),
                ratio     = round(observed[1] / mean(null, na.rm = TRUE), 2),
                p         = mean(null >= observed[1], na.rm = TRUE),
                n_perm    = .N), by = MG], row.names = FALSE)
  
  fig3 <- ggplot(out, aes(MG, null)) +
    geom_violin(fill = "#cfd4da", colour = NA, width = 0.8) +
    geom_point(data = obs, aes(MG, observed), colour = PAL$env, size = 2.2) +
    geom_text(data = obs, aes(MG, observed, label = "observed"),
              colour = PAL$env, vjust = -1, size = 2.4) +
    labs(title = NULL, y = "Mean silhouette width", x = NULL,
         caption = sprintf("Grey: %d within-location genotype permutations per MG", n_perm)) +
    theme_pub()
  save_tiff(fig3, "Fig3_Permutation_Null.tiff", w_mm = 129, h_mm = 90)
}
make_permnull(n_perm = 100)


message("\nDone. TIFFs are in: ", CFG$out_dir)
