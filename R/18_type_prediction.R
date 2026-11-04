if (!exists('ROOT')) source('R/00_setup.R')

## ============================================================================
## 18_type_prediction.R
## Do the environment types help predict how lines rank in a new trial?
##
## Each environment (location x year) of each maturity group is held out in
## turn. Without any data from the held-out environment, we recompute the
## adjusted deviations, the marker response profiles of all locations (the
## held-out location from its other years) and the Ward clusters at the
## reported number of types, and take the type of the held-out location.
## We then predict each line's yield in the held-out environment, centered
## on the environment mean, as its mean centered yield at a set of other
## locations, and record the Pearson correlation between predicted and
## observed values over the lines tested there.
##
## Sets of other locations compared (the held-out location is never used):
##   all        every other location
##   type       other locations of the same type
##   other      other locations of the other types
##   geographic the same number of locations as 'type', nearest in distance
##   profile    the same number of locations, most similar in profile
##   covariate  the same number of locations, nearest in county covariates
##              (US locations with covariates only)
##   random     the same number of locations drawn at random (N_RAND draws)
## Two sources of data at those locations:
##   same_year   the trials of the held-out year (the usual multi-location
##               trial of one season)
##   other_years the trials of all other years (a line's history)
## Correlations are compared over the same lines for all deterministic sets.
## The random sets give the null for 'type': a set of the same size chosen
## without regard to type. P for type against random is the share of N_RAND
## replicates whose mean correlation over environments is at least the
## observed one; type against all, geographic, profile and covariate sets
## are compared by paired Wilcoxon tests over environments.
## Also recorded: mean centered yield of the fifth of lines ranked highest by
## each predictor (realized selection differential).
##
## ENVCLASS_N_RAND (default 1000), ENVCLASS_CORES (default cores - 2, max 16)
## Writes out/Table_TypePrediction_Env.csv (one row per held-out environment
## and data source) and out/Table_TypePrediction_Summary.csv.
## Run time: about 10 to 30 minutes on 16 cores.
## ============================================================================

suppressPackageStartupMessages({ library(data.table); library(parallel) })
N_RAND    <- as.integer(Sys.getenv("ENVCLASS_N_RAND", "1000"))
N_CORES   <- as.integer(Sys.getenv("ENVCLASS_CORES", min(16, max(1, detectCores() - 2))))
MGS_RUN   <- strsplit(Sys.getenv("ENVCLASS_MGS", "0,I,II,III,IV"), ",")[[1]]
MAX_ENVS  <- as.integer(Sys.getenv("ENVCLASS_MAX_ENVS", "0"))   # >0 only for a quick test
MIN_LINES <- 10      # lines with predictions in a held-out environment
MIN_LOC   <- 30      # genotyped lines needed to profile a location (as R/01)
SEED      <- 2026
t_start <- Sys.time()

nust <- fread(file.path(DATA, "nust_ecw_subset.csv"))
nust[, loc_key := trimws(tolower(gsub(",.*$", "", Location)))]
gw  <- readRDS(file.path(DATA, "geno_worker_data.rds")); gdf <- gw$geno_df
strain_names <- names(gdf)[-(1:3)]
G <- as.matrix(gdf[, -(1:3)]); rownames(G) <- as.character(gdf$marker)
v2 <- readRDS(file.path(DATA, "marker_v2_results.rds"))
loc_xy <- nust[, .(lat = median(Lat, na.rm = TRUE), lon = median(Long, na.rm = TRUE)), by = loc_key]

## county covariates of the testing locations (as R/12: no yield history)
CACHE <- file.path(DATA, "nass_cache")
pw <- as.data.table(readRDS(file.path(CACHE, "nasa_power_counties.rds")))
we <- as.data.table(readRDS(file.path(CACHE, "weather_engineered.rds")))
so <- as.data.table(readRDS(file.path(CACHE, "ssurgo_multidepth.rds")))
mm <- as.data.table(readRDS(file.path(CACHE, "nass_management.rds")))
fx <- function(d) { n <- grep("(?i)^fips", names(d), value = TRUE)[1]; setnames(d, n, "fips")
                    d[, fips := sprintf("%05d", as.integer(fips))]; unique(d, by = "fips") }
pw <- fx(pw); we <- fx(we); so <- fx(so); mm <- fx(mm)
feat <- merge(pw, we[, c("fips", setdiff(names(we), names(pw))), with = FALSE], by = "fips", all.x = TRUE)
feat <- merge(feat, so, by = "fips", all.x = TRUE); feat <- merge(feat, mm, by = "fips", all.x = TRUE)
for (mc in intersect(c("irrig_fraction", "covercrop_fraction", "tile_drain_fraction"), names(feat)))
  feat[is.na(get(mc)), (mc) := 0]
lk <- fread(file.path(DATA, "county_training_locations.csv"), colClasses = list(character = "fips"))
loc_cov <- merge(unique(lk[, .(loc_key, fips)], by = "loc_key"), feat, by = "fips")

hav <- function(lat1, lon1, lat2, lon2) {
  r <- pi / 180; a <- sin((lat2 - lat1) * r / 2)^2 + cos(lat1 * r) * cos(lat2 * r) * sin((lon2 - lon1) * r / 2)^2
  12742 * asin(pmin(1, sqrt(a)))
}

profile_of <- function(dd) {                       # dd: Strain, value
  dd <- dd[, .(val = mean(value)), by = Strain]
  sc <- match(dd$Strain, strain_names); ok <- !is.na(sc)
  r <- suppressWarnings(cor(t(G[, sc[ok], drop = FALSE]), dd$val[ok]))
  r[is.na(r)] <- 0; as.numeric(r)
}

one_env <- function(E, mg_means, k, locs_all, covm) {
  Lh <- mg_means[Environment == E, loc_key][1]; Yh <- mg_means[Environment == E, Year][1]
  obs <- mg_means[Environment == E, .(Strain, e)]
  dd  <- mg_means[Environment != E]
  ## adjusted deviations without the held-out environment
  dd[, em := mean(YieldBuA), by = Environment]; dd[, sm := mean(YieldBuA), by = Strain]
  dd[, adj := YieldBuA - sm - em + mean(YieldBuA)]
  nloc <- dd[, .(n = uniqueN(Strain)), by = loc_key][n >= MIN_LOC & loc_key %in% locs_all]
  if (!(Lh %in% nloc$loc_key) || nrow(nloc) < k + 2) return(NULL)
  locs <- nloc$loc_key
  Pm <- do.call(rbind, lapply(locs, function(l) profile_of(dd[loc_key == l, .(Strain, value = adj)])))
  rownames(Pm) <- locs
  mv <- apply(Pm, 2, var); Pm <- Pm[, mv > 0, drop = FALSE]
  Pw <- sweep(Pm, 2, sqrt(apply(Pm, 2, var)), "*")
  cm <- cor(t(Pw)); cm[is.na(cm)] <- 0
  cl <- cutree(hclust(as.dist(1 - cm), method = "ward.D2"), k = k)
  ## centered yields at other locations for the lines of E
  lines <- obs$Strain
  res <- list()
  for (src in c("same_year", "other_years")) {
    sub <- if (src == "same_year") mg_means[Year == Yh & loc_key != Lh & Strain %in% lines]
           else mg_means[Year != Yh & loc_key != Lh & Strain %in% lines]
    if (nrow(sub) == 0) next
    X <- dcast(sub[, .(e = mean(e)), by = .(Strain, loc_key)], Strain ~ loc_key, value.var = "e")
    Xm <- as.matrix(X[, -1]); rownames(Xm) <- X$Strain
    O <- intersect(colnames(Xm), setdiff(locs, Lh))   # other profiled locations with data
    Tt <- O[cl[O] == cl[Lh]]; Ot <- setdiff(O, Tt); nt <- length(Tt)
    if (nt < 1 || length(Ot) < 1) next
    dgeo <- hav(loc_xy[loc_key == Lh, lat], loc_xy[loc_key == Lh, lon],
                loc_xy[match(O, loc_key), lat], loc_xy[match(O, loc_key), lon]); names(dgeo) <- O
    Gs <- names(sort(dgeo))[seq_len(nt)]
    Ps <- names(sort(cm[Lh, O], decreasing = TRUE))[seq_len(nt)]
    Cs <- NULL
    if (!is.null(covm) && Lh %in% rownames(covm)) {
      Oc <- intersect(O, rownames(covm))
      if (length(Oc) >= nt) {
        dc <- sqrt(colSums((t(covm[Oc, , drop = FALSE]) - covm[Lh, ])^2)); Cs <- names(sort(dc))[seq_len(nt)]
      }
    }
    pr <- function(S) { v <- rowMeans(Xm[, S, drop = FALSE], na.rm = TRUE); v[is.nan(v)] <- NA; v }
    y  <- obs$e[match(rownames(Xm), obs$Strain)]
    PRED <- list(all = pr(O), type = pr(Tt), other = pr(Ot), geographic = pr(Gs), profile = pr(Ps))
    if (!is.null(Cs)) PRED$covariate <- pr(Cs)
    ok <- Reduce(`&`, lapply(PRED[c("all", "type", "other", "geographic", "profile")], function(v) !is.na(v))) & !is.na(y)
    if (sum(ok) < MIN_LINES) next
    rr  <- vapply(PRED, function(v) { o2 <- ok & !is.na(v); if (sum(o2) < MIN_LINES) NA_real_ else cor(v[o2], y[o2]) }, 0)
    top <- vapply(PRED, function(v) { o2 <- ok & !is.na(v); if (sum(o2) < MIN_LINES) return(NA_real_)
                                      q <- quantile(v[o2], 0.8); mean(y[o2][v[o2] >= q]) - mean(y[o2]) }, 0)
    rnd <- vapply(seq_len(N_RAND), function(b) { S <- sample(O, nt); v <- pr(S); o2 <- ok & !is.na(v)
                                                 if (sum(o2) < MIN_LINES) NA_real_ else cor(v[o2], y[o2]) }, 0)
    ## overlap of the type set with the nearest-location sets, and the parts
    ## of the type and nearest sets that differ (locations in one set only)
    TnG <- setdiff(Tt, Gs); GnT <- setdiff(Gs, Tt)
    r_pair <- function(A, B) {
      if (length(A) == 0 || length(B) == 0) return(c(NA_real_, NA_real_))
      va <- pr(A); vb <- pr(B); o2 <- ok & !is.na(va) & !is.na(vb)
      if (sum(o2) < MIN_LINES) return(c(NA_real_, NA_real_))
      c(cor(va[o2], y[o2]), cor(vb[o2], y[o2]))
    }
    rd <- r_pair(TnG, GnT)
    res[[src]] <- list(row = data.table(environment = E, location = Lh, year = Yh, source = src,
                         lines = sum(ok), n_other = length(O), n_type = nt,
                         r_all = rr["all"], r_type = rr["type"], r_other = rr["other"],
                         r_geographic = rr["geographic"], r_profile = rr["profile"],
                         r_covariate = if ("covariate" %in% names(rr)) rr["covariate"] else NA_real_,
                         r_random_mean = mean(rnd, na.rm = TRUE),
                         sel_all = top["all"], sel_type = top["type"], sel_geographic = top["geographic"],
                         sel_profile = top["profile"],
                         overlap_geographic = length(intersect(Tt, Gs)) / nt,
                         overlap_profile = length(intersect(Tt, Ps)) / nt,
                         km_type = mean(dgeo[Tt]), km_geographic = mean(dgeo[Gs]), km_all = mean(dgeo),
                         n_type_only = length(TnG), r_type_only = rd[1], r_geographic_only = rd[2]),
                       rnd = rnd)
  }
  res
}

rows <- list(); summ <- list()
if (N_CORES > 1) {
  cl_par <- makeCluster(N_CORES)
  clusterEvalQ(cl_par, suppressPackageStartupMessages(library(data.table)))
  clusterExport(cl_par, c("G", "strain_names", "loc_xy", "hav", "profile_of", "one_env",
                          "N_RAND", "MIN_LINES", "MIN_LOC"), envir = environment())
}
for (mg in MGS_RUN) {
  r <- v2[[mg]]$F; if (is.null(r) || r$k < 2) next
  k <- r$k; locs_all <- names(r$clusters)
  mg_means <- nust[MG == mg & Strain %in% strain_names & !is.na(YieldBuA),
                   .(YieldBuA = mean(YieldBuA)), by = .(Strain, Environment, loc_key, Year)]
  mg_means[, e := YieldBuA - mean(YieldBuA), by = Environment]
  ## covariate matrix, standardized over the locations of this maturity group
  cv <- loc_cov[loc_key %in% locs_all]; cc <- setdiff(names(cv), c("fips", "loc_key"))
  cc <- cc[vapply(cc, function(z) is.numeric(cv[[z]]) && sum(!is.na(cv[[z]])) > 0.5 * nrow(cv), logical(1))]
  covm <- as.matrix(cv[, ..cc]); rownames(covm) <- cv$loc_key
  for (j in seq_len(ncol(covm))) { b <- is.na(covm[, j]); if (any(b)) covm[b, j] <- median(covm[, j], na.rm = TRUE) }
  covm <- covm[, apply(covm, 2, sd) > 0, drop = FALSE]; covm <- scale(covm)
  envs <- mg_means[loc_key %in% locs_all, .(n = .N), by = Environment][n >= MIN_LINES, Environment]
  if (MAX_ENVS > 0) { set.seed(SEED); envs <- sample(envs, min(MAX_ENVS, length(envs))) }
  cat(sprintf("MG %s: k = %d, %d held-out environments\n", mg, k, length(envs)))
  f_env <- function(E) tryCatch(one_env(E, mg_means, k, locs_all, covm), error = function(e) list(err = conditionMessage(e)))
  if (N_CORES > 1) {
    clusterExport(cl_par, c("mg_means", "k", "locs_all", "covm"), envir = environment())
    clusterSetRNGStream(cl_par, SEED)
    out <- parLapply(cl_par, envs, f_env)
  } else {
    set.seed(SEED); out <- lapply(envs, f_env)
  }
  errs <- vapply(out, function(o) if (!is.null(o$err)) o$err else NA_character_, "")
  if (any(!is.na(errs))) cat(sprintf("  %d environments failed; first error: %s\n", sum(!is.na(errs)), errs[!is.na(errs)][1]))
  out <- out[!vapply(out, function(o) is.null(o) || !is.null(o$err), logical(1))]
  cat(sprintf("  %d environments with results\n", length(out)))
  for (src in c("same_year", "other_years")) {
    o <- Filter(Negate(is.null), lapply(out, `[[`, src)); if (length(o) == 0) next
    R <- rbindlist(lapply(o, `[[`, "row")); R[, mg := mg]
    RN <- do.call(rbind, lapply(o, `[[`, "rnd"))           # environments x draws
    rows[[paste(mg, src)]] <- R
    null_mean <- colMeans(RN, na.rm = TRUE)
    wt <- function(a, b) { d <- R[[a]] - R[[b]]; d <- d[!is.na(d)]; if (length(d) < 5) NA_real_ else wilcox.test(d)$p.value }
    summ[[paste(mg, src)]] <- data.table(
      mg = mg, source = src, environments = nrow(R), lines_median = median(R$lines),
      n_type_median = median(R$n_type), n_other_median = median(R$n_other),
      r_all = mean(R$r_all), r_type = mean(R$r_type), r_other = mean(R$r_other),
      r_geographic = mean(R$r_geographic), r_profile = mean(R$r_profile),
      r_covariate = mean(R$r_covariate, na.rm = TRUE), covariate_envs = sum(!is.na(R$r_covariate)),
      r_random = mean(null_mean), r_random_q95 = quantile(null_mean, 0.95),
      p_type_vs_random = (1 + sum(null_mean >= mean(R$r_type))) / (length(null_mean) + 1),
      p_type_vs_all = wt("r_type", "r_all"), p_type_vs_other = wt("r_type", "r_other"),
      p_type_vs_geographic = wt("r_type", "r_geographic"), p_profile_vs_geographic = wt("r_profile", "r_geographic"),
      p_type_vs_covariate = wt("r_type", "r_covariate"),
      sel_all = mean(R$sel_all, na.rm = TRUE), sel_type = mean(R$sel_type, na.rm = TRUE),
      sel_geographic = mean(R$sel_geographic, na.rm = TRUE), sel_profile = mean(R$sel_profile, na.rm = TRUE),
      p_sel_type_vs_all = wt("sel_type", "sel_all"), p_sel_type_vs_geographic = wt("sel_type", "sel_geographic"),
      overlap_geographic = mean(R$overlap_geographic), overlap_profile = mean(R$overlap_profile),
      km_type = mean(R$km_type), km_geographic = mean(R$km_geographic), km_all = mean(R$km_all),
      disjoint_envs = sum(!is.na(R$r_type_only)),
      r_type_only = mean(R$r_type_only, na.rm = TRUE), r_geographic_only = mean(R$r_geographic_only, na.rm = TRUE),
      p_type_only_vs_geographic_only = wt("r_type_only", "r_geographic_only"),
      n_rand = N_RAND)
  }
  cat(sprintf("  done, %.1f min\n", as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
}
if (N_CORES > 1) stopCluster(cl_par)
env_tab <- rbindlist(rows, fill = TRUE); sum_tab <- rbindlist(summ, fill = TRUE)
fwrite(env_tab, file.path(OUT, "Table_TypePrediction_Env.csv"))
fwrite(sum_tab, file.path(OUT, "Table_TypePrediction_Summary.csv"))
print(sum_tab[, .(mg, source, environments, r_all = round(r_all, 3), r_type = round(r_type, 3),
                  r_other = round(r_other, 3), r_geo = round(r_geographic, 3), r_prof = round(r_profile, 3),
                  r_rand = round(r_random, 3), p_rand = signif(p_type_vs_random, 2),
                  p_all = signif(p_type_vs_all, 2), p_geo = signif(p_type_vs_geographic, 2))])
print(sum_tab[, .(mg, source, overlap_geo = round(overlap_geographic, 2), km_type = round(km_type),
                  km_geo = round(km_geographic), disjoint_envs, r_type_only = round(r_type_only, 3),
                  r_geo_only = round(r_geographic_only, 3), p_disjoint = signif(p_type_only_vs_geographic_only, 2))])
cat(sprintf("Done in %.1f min\n", as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
