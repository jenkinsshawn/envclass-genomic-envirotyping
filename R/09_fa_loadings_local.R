## ============================================================================
## 09_fa_loadings_local.R
## Factor-analytic grouping of locations for maturity groups II and III, for
## comparison with the marker-profile environment types.
##
## Model (one per maturity group, location set and number of factors):
##   y_ij = mu + loc_j + g_ij + e_ij
##   y_ij  = mean yield (bu/acre) of line i at location j over all years and
##           environments in the archive (stage-one line-by-location means)
##   loc_j = fixed location effect
##   g_ij  = random line effect at location j with factor-analytic covariance
##           among locations, Var(g_.j) = Lambda Lambda' + Psi (sommer::fam)
##   e_ij  = residual
## Locations are then clustered (Ward.D2, as for the profiles) on
##   (a) 1 - the fitted genetic correlation between locations, and
##   (b) Euclidean distance between correlation-scale loadings
##       (loadings divided by the square root of the location's genetic
##       variance; distances are unchanged by rotation of the factors),
## and compared with the published profile types by adjusted Rand index with a
## permutation reference for chance agreement.
##
## Two location sets per maturity group:
##   validation : the locations of Table 2 (reproduced in 08_shared_lines_comparison.R)
##   all        : all profiled locations. The factor-analytic model can
##                estimate a correlation for pairs with few or no shared lines
##                through the common factors, so this set shows how it places
##                the locations that response-distance clustering cannot.
##
## Requires: R >= 4.3, data.table, cluster, sommer (>= 4.4; installed below if
## missing). Run 08_shared_lines_comparison.R first (it writes the location sets).
## Expected run time: several minutes to more than an hour per fit
## (about 2,250 to 2,530 lines x 17 to 29 locations); each fit is logged.
## Run from the repository root:
##   Rscript R/09_fa_loadings_local.R            (MG II and III, 2 and 3 factors)
##   Rscript R/09_fa_loadings_local.R --mg=II --factors=2   (single fit)
## ============================================================================

if (!exists('ROOT')) source('R/00_setup.R')
args <- commandArgs(trailingOnly = TRUE)
argval <- function(name, default) {
  a <- grep(paste0('^--', name, '='), args, value = TRUE)
  if (length(a)) strsplit(sub('^[^=]*=', '', a[1]), ',')[[1]] else default
}
MGS     <- argval('mg', c('II', 'III'))
FACTORS <- as.integer(argval('factors', c('2', '3')))
SETS    <- argval('sets', c('validation', 'all'))
N_ITER  <- as.integer(argval('iters', '100'))
N_PERM  <- 9999

for (p in c('data.table', 'cluster', 'sommer'))
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p, repos = 'https://cloud.r-project.org')
suppressPackageStartupMessages({ library(data.table); library(cluster); library(sommer) })
set.seed(42)

OUT9 <- file.path(OUT, 'shared_lines', 'fa'); dir.create(OUT9, showWarnings = FALSE, recursive = TRUE)
logf <- file.path(OUT9, 'fa_run_log.txt')
logmsg <- function(...) { m <- sprintf(...); cat(m, '\n'); cat(m, '\n', file = logf, append = TRUE) }
logmsg('Started %s; sommer %s', format(Sys.time()), as.character(packageVersion('sommer')))

ari <- function(x, y) {
  tab <- table(x, y); n <- sum(tab)
  sc <- sum(choose(tab, 2)); sa <- sum(choose(rowSums(tab), 2)); sb <- sum(choose(colSums(tab), 2))
  e <- sa * sb / choose(n, 2); mx <- (sa + sb) / 2
  if (mx == e) return(NA_real_)
  (sc - e) / (mx - e)
}
ari_perm <- function(x, y) {
  obs <- ari(x, y); null <- replicate(N_PERM, ari(x, sample(y)))
  c(ari = obs, chance_q95 = unname(quantile(null, 0.95)),
    perm_p = (1 + sum(null >= obs - 1e-12)) / (N_PERM + 1))
}
ward_k <- function(d, k) cutree(hclust(d, method = 'ward.D2'), k = k)
sil_k <- function(d, max_k = 8) {
  hc <- hclust(d, method = 'ward.D2'); mk <- min(max_k, attr(d, 'Size') - 1)
  s <- sapply(2:mk, function(k) mean(silhouette(cutree(hc, k), d)[, 3])); which.max(s) + 1
}

v2   <- readRDS(file.path(DATA, 'marker_v2_results.rds'))
assn <- fread(file.path(OUT, 'shared_lines', 'Table_SL_Assignments.csv'))
nust <- fread(file.path(DATA, 'nust_ecw_subset.csv'))
nust[, loc_key := trimws(tolower(gsub(',.*$', '', Location)))]
nust[, MG := as.character(MG)]

rows <- list(); load_rows <- list()
for (mg in MGS) {
  prof_cl <- v2[[mg]]$F$clusters; k_prof <- v2[[mg]]$F$k
  for (sname in SETS) {
    vlocs <- assn$location[assn$mg == mg & assn$location_set == 'validation']
    locs <- if (sname == 'all') names(prof_cl) else vlocs
    d <- nust[MG == mg & !is.na(YieldBuA) & loc_key %in% locs,
              .(y = mean(YieldBuA)), by = .(Strain, loc_key)]
    d <- d[Strain %in% d[, .N, by = Strain][N >= 2, Strain]]
    d[, loc_key := factor(loc_key, levels = locs)][, Strain := factor(Strain)]
    setorder(d, loc_key, Strain)
    df <- as.data.frame(d)
    for (nf in FACTORS) {
      tag <- sprintf('MG%s_%s_FA%d', mg, sname, nf)
      logmsg('%s: %d lines, %d locations, %d line-location means; fitting...',
             tag, nlevels(d$Strain), length(locs), nrow(d))
      t1 <- Sys.time()
      fit <- tryCatch(mmes(y ~ loc_key, random = ~ vsm(fam(loc_key, nf), ism(Strain)),
                           rcov = ~ units, nIters = N_ITER, verbose = TRUE, data = df),
                      error = function(e) { logmsg('%s failed: %s', tag, conditionMessage(e)); NULL })
      mins <- as.numeric(difftime(Sys.time(), t1, units = 'mins'))
      logmsg('%s: finished in %.1f min', tag, mins)
      if (is.null(fit)) next
      saveRDS(fit, file.path(OUT9, paste0(tag, '_fit.rds')))
      fi <- loadings_mmes(fit, varianceScale = FALSE)   # Sigma proportional to Lambda Lambda' + Psi
      Lam <- as.matrix(fi$loadings); Psi <- as.numeric(fi$specific)
      G <- Lam %*% t(Lam) + diag(Psi); dimnames(G) <- list(rownames(Lam), rownames(Lam))
      ## loadings_mmes returns locations in factor-level order; check names
      if (!all(rownames(Lam) %in% locs)) rownames(Lam) <- rownames(G) <- colnames(G) <- locs
      C <- cov2cor(G)
      Ls <- Lam / sqrt(diag(G))
      d_cor <- as.dist(1 - C); d_load <- dist(Ls)
      pc <- prof_cl[rownames(C)]
      for (dn in c('fitted genetic correlation', 'correlation-scale loadings')) {
        dd <- if (dn == 'fitted genetic correlation') d_cor else d_load
        for (kk in unique(c(k_prof, sil_k(dd)))) {
          cl <- ward_k(dd, kk); ap <- ari_perm(pc, cl)
          rows[[paste(tag, dn, kk)]] <- data.table(mg = mg, location_set = sname, factors = nf,
            n_lines = nlevels(d$Strain), n_locations = length(locs), minutes = mins,
            AIC = if (is.null(fit$AIC)) NA_real_ else fit$AIC,
            clustering_basis = dn, k = kk, k_is_profile_k = kk == k_prof,
            sizes = paste(table(cl), collapse = '/'),
            ari = ap['ari'], chance_q95 = ap['chance_q95'], perm_p = ap['perm_p'],
            factor_share_of_genetic_variance = sum(Lam^2) / (sum(Lam^2) + sum(Psi)))
        }
      }
      load_rows[[tag]] <- data.table(mg = mg, location_set = sname, factors = nf,
        location = rownames(Lam), profile_type = prof_cl[rownames(Lam)],
        genetic_variance = diag(G), specific_variance = Psi,
        setNames(as.data.table(Ls), paste0('scaled_loading_', seq_len(nf))))
      fwrite(rbindlist(rows), file.path(OUT9, 'Table_FA_ARI.csv'))
      fwrite(rbindlist(load_rows, fill = TRUE), file.path(OUT9, 'Table_FA_Loadings.csv'))
    }
  }
}
logmsg('Done %s', format(Sys.time()))
