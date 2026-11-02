## ============================================================================
## 08_shared_lines_comparison.R
## Compares the marker-profile environment types with groupings that rely on
## lines shared between locations.
##   Part 1. Cluster locations directly on the phenotypic response distance
##           (1 - correlation of adjusted yield deviations over shared lines,
##           at least 8 shared lines), same method as the profile clusters
##           (Ward.D2, silhouette-chosen k over 2..8), compare by adjusted
##           Rand index with a label-permutation reference.
##   Part 2. Count profiled location pairs sharing < 8 and 0 lines, list the
##           locations that response-distance clustering cannot place, show
##           where the profiles place them, and check placements indirectly
##           through third locations that share lines with both.
##   Part 3. Principal-component fallback for the factor-analytic grouping
##           (MG II and III): eigen-decomposition of the pairwise correlation
##           matrix of line performance, Ward clustering of the loadings.
## Definitions are copied from R/01_marker_profiles.R and R/02_phenotypic_mantel.R.
## Run from the repository root:  Rscript R/08_shared_lines_comparison.R
## ============================================================================

if (!exists('ROOT')) source('R/00_setup.R')
suppressPackageStartupMessages({ library(data.table); library(cluster) })
set.seed(42)
t0 <- Sys.time()

OUT8 <- file.path(OUT, 'shared_lines'); dir.create(OUT8, showWarnings = FALSE, recursive = TRUE)
MGS <- c('0', 'I', 'II', 'III', 'IV')
MIN_PAIR <- 8          # paper's minimum number of shared lines per location pair
MAX_K <- 8             # same as PARAMS$max_k in 01_marker_profiles.R
N_PERM <- 9999
N_PERM_SMALL <- 999

## ---- helpers copied from 02_phenotypic_mantel.R --------------------------------
compute_adj_dev <- function(dat) {
  em <- dat[, .(em = mean(YieldBuA, na.rm = TRUE)), by = Environment]
  sm <- dat[, .(sm = mean(YieldBuA, na.rm = TRUE)), by = Strain]
  gm <- mean(dat$YieldBuA, na.rm = TRUE)
  dat <- merge(dat, em, by = 'Environment')
  dat <- merge(dat, sm, by = 'Strain')
  dat[, adj_dev := YieldBuA - sm - em + gm]
  dat
}
build_strain_loc_mat <- function(dt) {
  w <- dcast(dt, Strain ~ loc_key, value.var = 'value', fun.aggregate = mean)
  m <- as.matrix(w[, -1]); rownames(m) <- w$Strain
  m[is.nan(m)] <- NA
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

## ---- new helpers -----------------------------------------------------------------
ari <- function(x, y) {
  tab <- table(x, y); n <- sum(tab)
  sc <- sum(choose(tab, 2)); sa <- sum(choose(rowSums(tab), 2)); sb <- sum(choose(colSums(tab), 2))
  e <- sa * sb / choose(n, 2); mx <- (sa + sb) / 2
  if (mx == e) return(NA_real_)
  (sc - e) / (mx - e)
}
ari_perm <- function(x, y, n_perm = N_PERM) {
  obs <- ari(x, y)
  null <- replicate(n_perm, ari(x, sample(y)))
  list(ari = obs, null_mean = mean(null), null_q95 = unname(quantile(null, 0.95)),
       p = (1 + sum(null >= obs - 1e-12)) / (n_perm + 1))
}
ward_sil <- function(d, max_k = MAX_K) {
  hc <- hclust(d, method = 'ward.D2')
  mk <- min(max_k, attr(d, 'Size') - 1)
  sils <- sapply(2:mk, function(k) mean(silhouette(cutree(hc, k), d)[, 3]))
  list(hc = hc, k = which.max(sils) + 1, sil = max(sils), sils = setNames(sils, 2:mk))
}
mean_sil_at <- function(cl, d) mean(silhouette(cl, d)[, 3])

## ---- data -----------------------------------------------------------------------
v2   <- readRDS(file.path(DATA, 'marker_v2_results.rds'))
nust <- fread(file.path(DATA, 'nust_ecw_subset.csv'))
nust[, loc_key := trimws(tolower(gsub(',.*$', '', Location)))]
nust[, MG := as.character(MG)]
loc_coords <- nust[, .(lat = median(Lat, na.rm = TRUE), lon = median(Long, na.rm = TRUE)), by = loc_key]
loc_coords <- unique(loc_coords[!is.na(lat) & !is.na(lon)], by = 'loc_key')
## 02_phenotypic_mantel.R keeps only locations that join to a US county
## (sf/tigris, needs a download). Offline equivalent: the state code after the
## comma in Location is a US state abbreviation (drops ONT, QUE, MAN, SAS).
## The check against Table 1 location counts is printed below.
## Two location keys join a US and a Canadian site (Table_LocKey_Collisions.csv);
## the median coordinate, and so the county join, follows the site with more
## observations, so the state is taken from the most frequent Location string.
loc_state <- nust[, .N, by = .(loc_key, st = trimws(sub('^.*,', '', Location)))][order(-N)]
loc_state <- unique(loc_state, by = 'loc_key')
loc_us <- intersect(loc_coords[lat > 24 & lat < 50 & lon > -130 & lon < -60, loc_key],
                    loc_state[st %in% state.abb, loc_key])
SETS <- c(validation = 'US locations (validation set of Table 1)',
          all = 'all profiled locations, including Canada')

summ <- list(); ari_rows <- list(); pairs_rows <- list(); unplaced_rows <- list()
indirect_rows <- list(); pc_rows <- list(); assign_rows <- list(); loc_rows <- list()
store <- list()

for (mg in MGS) {
  cat(sprintf('\n==== MG %s ====\n', mg))
  res <- v2[[mg]]$F
  if (is.null(res)) next
  prof_cl <- res$clusters; prof_locs <- names(prof_cl); k_prof <- res$k
  d_prof <- as.matrix(res$dist)[prof_locs, prof_locs]

  ## response distance exactly as in 02_phenotypic_mantel.R (all lines in the MG)
  mg_means <- nust[MG == mg & !is.na(YieldBuA),
                   .(YieldBuA = mean(YieldBuA, na.rm = TRUE)),
                   by = .(Strain, Environment, loc_key, Year)]
  mg_means <- compute_adj_dev(mg_means)
  loc_in_prof <- intersect(unique(mg_means$loc_key), prof_locs)
  m_adj <- build_strain_loc_mat(mg_means[loc_key %in% loc_in_prof, .(Strain, loc_key, value = adj_dev)])
  d_ph <- pairwise_cor_dist(m_adj, MIN_PAIR)
  obs <- !is.na(m_adj) * 1
  n_sh <- crossprod(obs)                       # lines shared by each location pair
  d_ph <- d_ph[prof_locs, prof_locs]; n_sh <- n_sh[prof_locs, prof_locs]

  ## complete subset by the paper's pruning rule (drop the location with most
  ## missing pairs, repeat until no pair is missing)
  prune <- function(cands) {
    m <- d_ph[cands, cands, drop = FALSE]
    while (any(is.na(m)) && nrow(m) > 5) {
      worst <- which.max(rowSums(is.na(m)))
      m <- m[-worst, -worst, drop = FALSE]
    }
    m
  }
  m_sets <- list(validation = prune(intersect(prof_locs, loc_us)),
                 all = prune(prof_locs))
  ## Part 2 uses the all-location set, so that the only reason for exclusion
  ## is the shared-line requirement
  val_locs <- rownames(m_sets$all)
  unplaced <- setdiff(prof_locs, val_locs)
  cat(sprintf('  profiled %d (US %d), validation subset %d, all-location subset %d, not placed %d\n',
              length(prof_locs), length(intersect(prof_locs, loc_us)), nrow(m_sets$validation),
              length(val_locs), length(unplaced)))

  ## ---------------- Part 2 counts (all profiled locations) ----------------------
  ut <- upper.tri(n_sh)
  n_pairs <- sum(ut); lt8 <- sum(n_sh[ut] < MIN_PAIR); z0 <- sum(n_sh[ut] == 0)
  na_pairs <- sum(is.na(d_ph[ut]))
  pairs_rows[[mg]] <- data.table(mg = mg, profiled_locations = length(prof_locs),
    location_pairs = n_pairs, pairs_lt8_shared = lt8, pairs_zero_shared = z0,
    pairs_without_response_distance = na_pairs,
    median_shared_lines = median(n_sh[ut]),
    locations_with_any_pair_lt8 = sum(rowSums(n_sh < MIN_PAIR) > 0),
    locations_in_complete_subset = length(val_locs),
    locations_not_placed = length(unplaced))

  ## per-location table
  types <- sort(unique(prof_cl))
  for (L in prof_locs) {
    others <- setdiff(prof_locs, L)
    md <- sapply(types, function(t) { mm <- setdiff(names(prof_cl)[prof_cl == t], L)
      if (length(mm)) mean(d_prof[L, mm]) else NA_real_ })
    loc_rows[[paste(mg, L)]] <- data.table(mg = mg, location = L, profile_type = prof_cl[[L]],
      in_complete_subset = L %in% val_locs,
      lines_at_location = sum(obs[, L]),
      partners_lt8_shared = sum(n_sh[L, others] < MIN_PAIR),
      partners_zero_shared = sum(n_sh[L, others] == 0),
      nearest_type_by_profile = types[which.min(md)],
      t(setNames(round(md, 4), paste0('profile_dist_type', types))))
  }

  ## ---------------- Part 1: cluster on response distance -----------------------
  store[[mg]] <- list(d_prof = d_prof, prof_cl = prof_cl, n_sh = n_sh, d_ph = d_ph,
                      unplaced = unplaced, val_locs = val_locs, sets = list())
  for (sname in names(SETS)) {
    m_p <- m_sets[[sname]]; sl <- rownames(m_p)
    if (length(sl) < 6) next
    dd <- as.dist(m_p)
    ws <- ward_sil(dd)
    resp_k_prof <- cutree(ws$hc, k = min(k_prof, length(sl) - 1))
    resp_k_sil  <- cutree(ws$hc, k = ws$k)
    pc_pub <- prof_cl[sl]
    ## profiles re-clustered on the same subset (same method, same k)
    pws <- ward_sil(as.dist(d_prof[sl, sl]))
    prof_sub_k <- cutree(pws$hc, k = k_prof)
    cmp <- list(
      c('published profile types', 'response distance, k = profile k', 'pp_rk'),
      c('published profile types', 'response distance, silhouette k', 'pp_rs'),
      c('profiles reclustered on subset, k = profile k', 'response distance, k = profile k', 'ps_rk'))
    for (cc in cmp) {
      x <- switch(substr(cc[3], 1, 2), pp = pc_pub, ps = prof_sub_k[sl])
      y <- switch(substr(cc[3], 4, 5), rk = resp_k_prof[sl], rs = resp_k_sil[sl])
      ap <- ari_perm(x, y)
      ari_rows[[paste(mg, sname, cc[3])]] <- data.table(mg = mg, location_set = sname,
        n_locations = length(sl), profile_grouping = cc[1], response_grouping = cc[2],
        k_profile = length(unique(x)), k_response = length(unique(y)),
        sizes_profile = paste(table(x), collapse = '/'), sizes_response = paste(table(y), collapse = '/'),
        ari = ap$ari, chance_mean = ap$null_mean, chance_q95 = ap$null_q95, perm_p = ap$p)
    }
    ## Do profile types separate locations on response distance? Mean response
    ## distance within types vs between types, permutation of type labels.
    Mp <- as.matrix(dd); utp <- upper.tri(Mp)
    wb <- function(lab) { same <- outer(lab, lab, '=='); c(mean(Mp[utp & same]), mean(Mp[utp & !same])) }
    o <- wb(pc_pub); obs_diff <- o[2] - o[1]
    null_diff <- replicate(N_PERM, { z <- wb(sample(pc_pub)); z[2] - z[1] })
    summ[[paste(mg, sname)]] <- data.table(mg = mg, location_set = sname,
      profiled_locations = length(prof_locs), k_profile = k_prof, profile_silhouette = res$sil,
      locations_compared = length(sl),
      response_k_silhouette = ws$k, response_silhouette_at_its_k = ws$sil,
      response_silhouette_at_profile_k = mean_sil_at(resp_k_prof, dd),
      profile_types_silhouette_on_response_distance = mean_sil_at(pc_pub, dd),
      response_dist_within_types = o[1], response_dist_between_types = o[2],
      between_minus_within = obs_diff, between_minus_within_chance_q95 = unname(quantile(null_diff, 0.95)),
      between_minus_within_perm_p = (1 + sum(null_diff >= obs_diff - 1e-12)) / (N_PERM + 1))
    assign_rows[[paste(mg, sname)]] <- data.table(mg = mg, location_set = sname, location = sl,
      profile_type = pc_pub, response_group_k_profile = resp_k_prof[sl], response_group_k_sil = resp_k_sil[sl])
    store[[mg]]$sets[[sname]] <- list(locs = sl, m_p = m_p, resp = resp_k_prof, hc = ws$hc)
  }

  ## ---------------- Part 2: placement of locations not placed -------------------
  ## (a) nearest profile type by response distance, using only the pairs that
  ##     share >= 8 lines with locations of the complete subset
  ref <- val_locs
  near_resp <- function(L, lab) {
    md <- sapply(types, function(t) { mm <- intersect(names(lab)[lab == t], ref)
      v <- d_ph[L, mm]; v <- v[!is.na(v)]; if (length(v)) mean(v) else NA_real_ })
    if (sum(!is.na(md)) < 2) return(NA_integer_)
    types[which.min(md)]
  }
  if (length(unplaced) && length(ref) >= 6) {
    nr <- sapply(unplaced, near_resp, lab = prof_cl)
    n_ref <- sapply(unplaced, function(L) sum(!is.na(d_ph[L, ref])))
    agree <- nr == prof_cl[unplaced]
    ## chance reference: shuffle the type labels of the reference locations
    lab_ref <- prof_cl[ref]
    null_rate <- replicate(N_PERM_SMALL, {
      lab <- prof_cl; lab[ref] <- sample(lab_ref)
      nrp <- sapply(unplaced, near_resp, lab = lab)
      mean(nrp == prof_cl[unplaced], na.rm = TRUE) })
    obs_rate <- mean(agree, na.rm = TRUE)
    unplaced_rows[[mg]] <- data.table(mg = mg, location = unplaced, profile_type = prof_cl[unplaced],
      comparable_reference_locations = n_ref, nearest_type_by_response = nr, agrees = agree,
      agreement_rate = obs_rate, chance_rate_mean = mean(null_rate, na.rm = TRUE),
      perm_p = (1 + sum(null_rate >= obs_rate - 1e-12, na.rm = TRUE)) / (N_PERM_SMALL + 1))
  }

  ## (b) pairs sharing < 8 lines: compare through third locations C that share
  ##     >= 8 lines with both members. Indirect dissimilarity = 1 - correlation
  ##     of the two members' response distances to the common third locations.
  pr <- which(ut, arr.ind = TRUE)
  ind <- rbindlist(lapply(seq_len(nrow(pr)), function(r) {
    a <- prof_locs[pr[r, 1]]; b <- prof_locs[pr[r, 2]]
    cset <- setdiff(prof_locs, c(a, b))
    cset <- cset[!is.na(d_ph[a, cset]) & !is.na(d_ph[b, cset])]
    ic <- if (length(cset) >= 5 && sd(d_ph[a, cset]) > 0 && sd(d_ph[b, cset]) > 0)
      1 - cor(d_ph[a, cset], d_ph[b, cset]) else NA_real_
    data.table(mg = mg, loc_a = a, loc_b = b, shared_lines = n_sh[a, b],
      comparable_directly = !is.na(d_ph[a, b]), response_distance = d_ph[a, b],
      n_third_locations = length(cset), indirect_dissimilarity = ic,
      profile_distance = d_prof[a, b], same_profile_type = prof_cl[[a]] == prof_cl[[b]])
  }))
  indirect_rows[[mg]] <- ind

  ## ---------------- Part 3: PC fallback (MG II, III) ----------------------------
  if (mg %in% c('II', 'III')) for (sname in names(store[[mg]]$sets)) {
    S <- store[[mg]]$sets[[sname]]; sl <- S$locs
    R <- 1 - S$m_p                             # pairwise correlation of line performance
    eg <- eigen(R, symmetric = TRUE)
    lam <- pmax(eg$values, 0)
    n_neg <- sum(eg$values < -1e-10)
    for (nf in 2:3) {
      Lmat <- eg$vectors[, 1:nf, drop = FALSE] %*% diag(sqrt(lam[1:nf]))
      rownames(Lmat) <- sl
      clL <- cutree(hclust(dist(Lmat), method = 'ward.D2'), k = k_prof)
      ap <- ari_perm(prof_cl[sl], clL)
      ap2 <- ari_perm(S$resp[sl], clL)
      pc_rows[[paste(mg, sname, nf)]] <- data.table(mg = mg, location_set = sname, n_locations = length(sl),
        components = nf, variance_share = sum(lam[1:nf]) / sum(lam), negative_eigenvalues = n_neg,
        k = k_prof, sizes = paste(table(clL), collapse = '/'),
        ari_vs_profile_types = ap$ari, chance_q95 = ap$null_q95, perm_p = ap$p,
        ari_vs_response_clusters = ap2$ari)
    }
  }
}

tabs <- list(
  Table_SL_Summary = rbindlist(summ, fill = TRUE),
  Table_SL_ARI = rbindlist(ari_rows),
  Table_SL_PairCounts = rbindlist(pairs_rows),
  Table_SL_Locations = rbindlist(loc_rows, fill = TRUE),
  Table_SL_Assignments = rbindlist(assign_rows),
  Table_SL_UnplacedLocations = rbindlist(unplaced_rows, fill = TRUE),
  Table_SL_PairsIndirect = rbindlist(indirect_rows),
  Table_SL_PCFallback = rbindlist(pc_rows))
for (nm in names(tabs)) fwrite(tabs[[nm]], file.path(OUT8, paste0(nm, '.csv')))

## indirect check summary: does profile distance track the indirect dissimilarity
## for pairs that cannot be compared directly? Permutation of location labels.
ind_all <- tabs$Table_SL_PairsIndirect
ind_summ <- rbindlist(lapply(MGS, function(g) {
  x <- ind_all[ind_all$mg == g & !is.na(ind_all$indirect_dissimilarity)]
  if (!nrow(x) || is.null(store[[g]])) return(NULL)
  mg <- g
  dprof <- store[[g]]$d_prof
  locs <- rownames(dprof)
  out <- list()
  for (grp in c('not comparable directly', 'comparable directly')) {
    y <- if (grp == 'not comparable directly') x[comparable_directly == FALSE] else x[comparable_directly == TRUE]
    if (nrow(y) < 5) { out[[grp]] <- data.table(mg = mg, pairs = grp, n_pairs = nrow(y)); next }
    rho <- cor(y$profile_distance, y$indirect_dissimilarity, method = 'spearman')
    ia <- match(y$loc_a, locs); ib <- match(y$loc_b, locs)
    null <- replicate(N_PERM_SMALL, { p <- sample(length(locs))
      cor(dprof[cbind(p[ia], p[ib])], y$indirect_dissimilarity, method = 'spearman') })
    rho_dir <- if (grp == 'comparable directly')
      cor(y$response_distance, y$indirect_dissimilarity, method = 'spearman') else NA_real_
    out[[grp]] <- data.table(mg = mg, pairs = grp, n_pairs = nrow(y),
      spearman_profile_vs_indirect = rho, perm_p = (1 + sum(null >= rho)) / (N_PERM_SMALL + 1),
      spearman_direct_vs_indirect = rho_dir,
      mean_indirect_same_type = mean(y$indirect_dissimilarity[y$same_profile_type]),
      mean_indirect_diff_type = mean(y$indirect_dissimilarity[!y$same_profile_type]),
      n_same_type = sum(y$same_profile_type))
  }
  rbindlist(out, fill = TRUE)
}), fill = TRUE)
fwrite(ind_summ, file.path(OUT8, 'Table_SL_IndirectSummary.csv'))
saveRDS(store, file.path(OUT8, 'shared_lines_store.rds'))

tc <- fread(file.path(OUT, 'Table_PhenoMantel_Final.csv'))[, .(mg, n_locs)]
chk <- merge(tc, tabs$Table_SL_Summary[location_set == 'validation', .(mg, locations_compared)], by = 'mg')
cat('\nCheck against Table 1 location counts:\n'); print(chk)
options(width = 250)
print(tabs$Table_SL_Summary); print(tabs$Table_SL_ARI); print(tabs$Table_SL_PairCounts)
print(tabs$Table_SL_UnplacedLocations); print(ind_summ); print(tabs$Table_SL_PCFallback)
cat(sprintf('\nRuntime %.1f min\n', as.numeric(difftime(Sys.time(), t0, units = 'mins'))))
