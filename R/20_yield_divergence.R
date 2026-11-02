if (!exists('ROOT')) source('R/00_setup.R')

## ============================================================================
## 20_yield_divergence.R
## Yield divergence between genomic environment types and the yield penalty
## for ignoring type structure (Methods, "Yield divergence"; Results, same).
##
## Reconstructed: no script for these numbers survives. Follows the Methods
## text on the archived data and the archived cluster assignments
## (data/marker_v2_results.rds, Method F).
##
## Primary definitions (as written in Methods):
##   Divergence, per pair of types (A, B) within a maturity group:
##     for each line tested at locations of both types, the mean adjusted
##     deviation over its environments in A minus that in B; the absolute
##     value is averaged over lines. Reported in kg/ha.
##   Penalty, per year and type t:
##     yield in t of the best line in t, minus yield in t of the best line
##     over all typed locations that year; averaged over years and types.
##     "Best line" = highest mean yield (line-environment means averaged over
##     the environments of that year). Only lines tested in every type that
##     year are candidates, so the overall best line has a yield in each type.
##
## Sensitivity variants (column 'variant'):
##   div_adjdev     : primary divergence on adjusted deviations
##   div_rawyield   : same on raw line-environment mean yields
##   div_adjdev_min2: lines with >= 2 environments in each type
##   div_within_year: line-by-year means; pairs formed within a year
##   pen_allcand    : penalty with all lines as candidates (best line in t
##                    from lines tested in t; overall best from lines in all
##                    types)
##   pen_adjdev     : penalty on adjusted deviations instead of raw yield
##
## Output: out/Table_YieldDivergence.csv (+ _pairs), out/Table_YieldPenalty.csv,
##         out/Table_YieldPenalty_byYear.csv
## ============================================================================

library(data.table)
set.seed(42)

BU_TO_KGHA <- 67.25   # soybean, 60 lb per bushel: 60 * 0.45359 / 0.40469

nust <- fread(file.path(DATA, 'nust_ecw_subset.csv'))
nust[, loc_key := trimws(tolower(gsub(',.*$', '', Location)))]
gw <- readRDS(file.path(DATA, 'geno_worker_data.rds'))
genotyped <- intersect(unique(nust$Strain), names(gw$geno_df)[-(1:3)])
rm(gw)
mv2 <- readRDS(file.path(DATA, 'marker_v2_results.rds'))

## Same construction as R/01_marker_profiles.R
compute_adj_dev <- function(dat) {
  em <- dat[, .(env_mean = mean(YieldBuA, na.rm = TRUE)), by = Environment]
  sm <- dat[, .(strain_mean = mean(YieldBuA, na.rm = TRUE)), by = Strain]
  gm <- mean(dat$YieldBuA, na.rm = TRUE)
  dat <- merge(dat, em, by = 'Environment')
  dat <- merge(dat, sm, by = 'Strain')
  dat[, adj_dev := YieldBuA - strain_mean - env_mean + gm]
  dat
}

mg_typed_means <- function(mg_val) {
  mg_dat <- nust[MG == mg_val & Strain %in% genotyped & !is.na(YieldBuA)]
  mm <- mg_dat[, .(YieldBuA = mean(YieldBuA, na.rm = TRUE)),
               by = .(Strain, Environment, loc_key, Year)]
  mm <- compute_adj_dev(mm)
  cl <- mv2[[mg_val]]$F$clusters
  mm[, type := cl[loc_key]]
  mm[!is.na(type)]
}

## ---- divergence -------------------------------------------------------------
divergence_pairs <- function(mm, value, min_env = 1, by_year = FALSE) {
  grp <- if (by_year) c('Strain', 'Year', 'type') else c('Strain', 'type')
  lt <- mm[, .(v = mean(get(value)), n_env = .N), by = grp]
  lt <- lt[n_env >= min_env]
  types <- sort(unique(mm$type))
  prs <- combn(types, 2, simplify = FALSE)
  key <- if (by_year) c('Strain', 'Year') else 'Strain'
  rbindlist(lapply(prs, function(p) {
    a <- lt[type == p[1]]; b <- lt[type == p[2]]
    ab <- merge(a, b, by = key, suffixes = c('.a', '.b'))
    data.table(type_a = p[1], type_b = p[2], n_lines = uniqueN(ab$Strain),
               n_units = nrow(ab), mad_bu = mean(abs(ab$v.a - ab$v.b)))
  }))
}

## ---- penalty ----------------------------------------------------------------
## pairwise = TRUE: each pair of types is treated as its own comparison
## (years in which both types were tested), which matters only for MG 0 and
## MG III; with two types it is identical to the primary definition.
penalty_years <- function(mm, value, all_candidates = FALSE, pairwise = FALSE) {
  tc <- as.character(sort(unique(mm$type)))
  sets <- if (pairwise) combn(tc, 2, simplify = FALSE) else list(tc)
  out <- list()
  for (s in sets) for (yr in sort(unique(mm$Year))) {
    yd <- mm[Year == yr & as.character(type) %in% s]
    if (uniqueN(yd$type) < length(s)) next
    lt <- yd[, .(v = mean(get(value))), by = .(Strain, type)]
    w <- dcast(lt, Strain ~ type, value.var = 'v')
    in_all <- complete.cases(w[, ..s])
    if (sum(in_all) < 2) next
    w_all <- w[in_all]
    ov <- yd[Strain %in% w_all$Strain, .(v = mean(get(value))), by = Strain]
    overall_best <- ov$Strain[which.max(ov$v)]
    for (t in s) {
      cand <- if (all_candidates) w[!is.na(get(t))] else w_all
      tb <- cand$Strain[which.max(cand[[t]])]
      out[[length(out) + 1]] <- data.table(
        set = paste(s, collapse = '-'),
        Year = yr, type = t, n_candidates = nrow(cand), n_in_all = nrow(w_all),
        best_in_type = tb, overall_best = overall_best,
        penalty_bu = cand[[t]][cand$Strain == tb] - w_all[[t]][w_all$Strain == overall_best])
    }
  }
  rbindlist(out)
}

summ_pen <- function(py, mg_val, variant) {
  if (!nrow(py)) return(NULL)
  ## a year counts as "best line differs between types" if the type-specific
  ## best lines are not all the same line
  yd <- py[, .(differs = uniqueN(best_in_type) > 1,
               mean_pen = mean(penalty_bu), max_pen = max(penalty_bu)), by = .(set, Year)]
  data.table(mg = mg_val, variant = variant, n_years = nrow(yd),
             pct_years_best_differs = 100 * mean(yd$differs),
             mean_penalty_kgha = BU_TO_KGHA * mean(py$penalty_bu),
             max_penalty_kgha = BU_TO_KGHA * max(py$penalty_bu),
             max_year = py$Year[which.max(py$penalty_bu)],
             median_candidates = median(py$n_in_all))
}

## ---- permutation null for divergence -----------------------------------------
## A mean absolute difference between two noisy line means is positive even
## when types do not differ. Null: type labels permuted among the typed
## locations (group sizes kept), divergence recomputed. One-sided p.
N_PERM <- 999
div_null <- function(mm, n_perm = N_PERM) {
  locs <- unique(mm[, .(loc_key, type)])
  obs <- mean(divergence_pairs(mm, 'adj_dev')$mad_bu)
  nul <- replicate(n_perm, {  # pairs without shared lines are skipped
    pl <- setNames(sample(locs$type), locs$loc_key)
    m2 <- copy(mm)[, type := pl[loc_key]]
    mean(divergence_pairs(m2, 'adj_dev')$mad_bu, na.rm = TRUE)
  })
  c(observed_kgha = obs * BU_TO_KGHA, null_mean_kgha = mean(nul) * BU_TO_KGHA,
    null_p95_kgha = quantile(nul, 0.95, names = FALSE) * BU_TO_KGHA,
    excess_kgha = (obs - mean(nul)) * BU_TO_KGHA,
    p = (1 + sum(nul >= obs)) / (n_perm + 1))
}

div_out <- list(); pen_out <- list(); pen_year <- list(); null_out <- list()
for (mg_val in c('0', 'I', 'II', 'III', 'IV')) {
  mm <- mg_typed_means(mg_val)
  v <- list(
    div_adjdev      = divergence_pairs(mm, 'adj_dev'),
    div_rawyield    = divergence_pairs(mm, 'YieldBuA'),
    div_adjdev_min2 = divergence_pairs(mm, 'adj_dev', min_env = 2),
    div_within_year = divergence_pairs(mm, 'adj_dev', by_year = TRUE))
  for (nm in names(v)) div_out[[paste(mg_val, nm)]] <-
    cbind(mg = mg_val, variant = nm, v[[nm]][, mad_kgha := mad_bu * BU_TO_KGHA])
  p1 <- penalty_years(mm, 'YieldBuA')
  p2 <- penalty_years(mm, 'YieldBuA', all_candidates = TRUE)
  p3 <- penalty_years(mm, 'adj_dev')
  p4 <- penalty_years(mm, 'YieldBuA', pairwise = TRUE)
  pen_out[[mg_val]] <- rbind(summ_pen(p1, mg_val, 'pen_primary'),
                             summ_pen(p2, mg_val, 'pen_allcand'),
                             summ_pen(p3, mg_val, 'pen_adjdev'),
                             summ_pen(p4, mg_val, 'pen_pairwise'))
  pen_year[[mg_val]] <- cbind(mg = mg_val, p1)
  null_out[[mg_val]] <- cbind(mg = mg_val, as.data.table(as.list(div_null(mm, N_PERM))))
}
div_dt <- rbindlist(div_out)
## maturity-group summary: mean over type pairs
div_mg <- div_dt[, .(n_pairs = .N, mean_kgha = mean(mad_kgha),
                     min_pair_kgha = min(mad_kgha), max_pair_kgha = max(mad_kgha),
                     n_lines_min = min(n_lines)), by = .(mg, variant)]
pen_dt <- rbindlist(pen_out)

fwrite(div_dt, file.path(OUT, 'Table_YieldDivergence_pairs.csv'))
fwrite(div_mg, file.path(OUT, 'Table_YieldDivergence.csv'))
fwrite(pen_dt, file.path(OUT, 'Table_YieldPenalty.csv'))
fwrite(rbindlist(pen_year), file.path(OUT, 'Table_YieldPenalty_byYear.csv'))
fwrite(rbindlist(null_out), file.path(OUT, 'Table_YieldDivergence_Null.csv'))
print(div_mg); print(pen_dt); print(rbindlist(null_out))
