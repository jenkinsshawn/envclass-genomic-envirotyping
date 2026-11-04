if (!exists('ROOT')) source('R/00_setup.R')

## ============================================================================
## 11_partial_mantel_power.R
## Power of the partial Mantel test (envirotyping | genomic) to detect an
## independent envirotyping contribution, MG I to IV.
##
## Clean version of the simulation run from the console on 2026-07-16
## ("LAST.R", recovered from RStudio history; its output was printed only).
##
## 1. The Table 1 distance matrices (genomic profile, phenotypic response,
##    envirotyping, geographic) are rebuilt by running R/02_phenotypic_mantel.R
##    from a temporary copy that stops before its first fwrite(), with one line
##    added that keeps the matrices. Nothing in out/ is overwritten by 02.
##    The rebuilt matrices are cached in out/power/table1_distances.rds.
## 2. For each MG, standardized response distance is regressed on standardized
##    genomic and envirotyping distance. Simulated response distance is
##      alpha * z(gen) + beta * z(env) + sigma * z(noise)
##    where alpha and sigma are the fitted genomic loading and residual SD and
##    beta (the envirotyping loading) is swept from 0 to 0.55. The noise is the
##    standardized Euclidean distance among the locations on 10 random normal
##    features, so it has the dependence structure of a distance matrix.
## 3. Each simulated matrix is tested with vegan::mantel.partial(env, sim, gen)
##    at 199 permutations; detection is p < 0.05 (vegan's p includes the
##    observed statistic). Power is the detection rate over B_SIM simulations.
##    The row with beta = 0 is the type I error rate.
## 4. The effect size reported is the mean simulated partial r at each beta.
##    The 80% threshold is given two ways: the mean partial r at the smallest
##    beta with power >= 0.80 (the console-script rule), and the partial r at
##    which the power curve crosses 0.80 by linear interpolation.
##
## Run from the repository root (see run_partial_mantel_power.ps1):
##   Rscript R/11_partial_mantel_power.R
## Environment variables: ENVCLASS_B_SIM (default 300), ENVCLASS_SIM_PERM
## (default 199), ENVCLASS_N_CORES (default 1; parallel over MG x beta),
## ENVCLASS_LOC_FIPS (optional path to a location -> county FIPS table, used
## instead of the sf/tigris county join; leave unset on a full R install).
## Needs vegan, ecodist, sf, tigris (02 joins locations to US counties with a
## census download on first use).
##
## Writes out/power/Table_PartialMantel_Power.csv       (MG x beta)
##        out/power/Table_PartialMantel_Power80.csv     (80% thresholds)
##        out/power/Table_PartialMantel_Table1Check.csv (rebuilt vs Table 1)
## ============================================================================

suppressPackageStartupMessages({ library(data.table); library(vegan); library(parallel) })

B_SIM    <- as.integer(Sys.getenv("ENVCLASS_B_SIM", "300"))
SIM_PERM <- as.integer(Sys.getenv("ENVCLASS_SIM_PERM", "199"))
N_CORES  <- as.integer(Sys.getenv("ENVCLASS_N_CORES", "1"))
SEED     <- 7                      # as in the 2026-07-16 console run
BETAS    <- seq(0, 0.55, by = 0.05)
N_NOISE  <- 10
MGS_D    <- c("I", "II", "III", "IV")
OUT_PW   <- file.path(OUT, "power")
dir.create(OUT_PW, showWarnings = FALSE, recursive = TRUE)
t_start <- Sys.time()
cat(sprintf("Partial Mantel power: B_SIM = %d, permutations = %d, cores = %d\n",
            B_SIM, SIM_PERM, N_CORES))

## ---- 1. Table 1 distance matrices -----------------------------------------
dist_cache <- file.path(OUT_PW, "table1_distances_v112.rds")
if (file.exists(dist_cache)) {
  DISTS <- readRDS(dist_cache)
  cat("Loaded cached Table 1 distance matrices\n")
} else {
  src02 <- readLines(file.path(ROOT, "R", "02_phenotypic_mantel.R"), warn = FALSE)
  anchor <- grep("d_geo\\s*<- clean_dist", src02)[1]
  cut    <- grep("^fwrite", src02)[1] - 1
  stopifnot(!is.na(anchor), !is.na(cut), anchor < cut)
  inj <- "  DISTS[[mg_val]] <<- list(gen = d_gen, pheno = d_pheno, env = d_env, geo = d_geo, n = length(shared), locs = shared)"
  src_run <- append(src02[1:cut], inj, after = anchor)
  LOC_FIPS <- Sys.getenv("ENVCLASS_LOC_FIPS", "")
  if (nzchar(LOC_FIPS)) {
    ## Offline mode, for machines where sf/tigris/ecodist cannot be loaded:
    ## the county join is read from a table of location -> county FIPS
    ## (data_offline/loc_fips_offline.csv) instead of st_join on the 2022
    ## cartographic county boundaries, and the MRM fit (not used here) is
    ## skipped. Check Table_PartialMantel_Table1Check.csv against Table 1.
    i0 <- grep("^suppressPackageStartupMessages\\(\\{", src_run)[1]
    i1 <- grep('^for \\(pkg in c\\("tigris", "sf"', src_run)[1]
    stopifnot(!is.na(i0), !is.na(i1), i0 < i1)
    src_run <- c(src_run[1:(i0 - 1)],
      sprintf('loc_fips <- fread("%s", colClasses = list(character = "fips"))[loc_key %%in%% loc_coords[lat > 24 & lat < 50 & lon > -130 & lon < -60, loc_key], .(loc_key, fips)]', LOC_FIPS),
      src_run[(i1 + 1):length(src_run)])
    j <- grep("mrm_fit <- ecodist::MRM", src_run)
    src_run[j] <- "  mrm_fit <- list(r.squared = NA_real_, coef = matrix(NA_real_, 4, 2))"
    cat("Offline county join from", LOC_FIPS, "\n")
  }
  tmp <- tempfile(fileext = ".R")
  writeLines(src_run, tmp)
  DISTS <- list()
  e02 <- new.env(parent = globalenv())
  e02$DISTS <- list()
  ## the injected line assigns into the first DISTS found from run_one_mg's
  ## enclosure, which is e02
  source(tmp, local = e02, echo = FALSE)
  DISTS <- e02$DISTS
  tab1 <- e02$out
  saveRDS(DISTS, dist_cache)
  fwrite(tab1[, .(mg, n_locs, mantel_genomic, mantel_enviro, partial_env_gen, partial_env_gen_p)],
         file.path(OUT_PW, "Table_PartialMantel_Table1Check.csv"))
  rm(e02)
}
DISTS <- DISTS[intersect(MGS_D, names(DISTS))]
print(data.table(mg = names(DISTS), n = sapply(DISTS, `[[`, "n"),
                 pairs = sapply(DISTS, function(d) length(d$gen))))

## ---- 2-3. simulation ---------------------------------------------------------
zv <- function(d) as.vector(scale(as.vector(d)))
sim_fits <- rbindlist(lapply(names(DISTS), function(mg) {
  D <- DISTS[[mg]]
  fit <- lm(zv(D$pheno) ~ zv(D$gen) + zv(D$env))
  data.table(mg = mg, n = D$n, alpha = unname(coef(fit)[2]),
             beta_obs = unname(coef(fit)[3]), sigma = sd(resid(fit)),
             r_gen_env = cor(as.vector(D$gen), as.vector(D$env)),
             partial_obs = vegan::mantel.partial(D$env, D$pheno, D$gen, permutations = 0)$statistic)
}))
print(sim_fits)

one_cell <- function(task) {
  D <- DISTS[[task$mg]]; n <- D$n; fp <- sim_fits[mg == task$mg]
  zg <- zv(D$gen); ze <- zv(D$env)
  mk <- function(v) { d <- D$gen; d[] <- v; d }
  o <- t(replicate(B_SIM, {
    zn <- zv(dist(matrix(rnorm(n * N_NOISE), n, N_NOISE)))
    m  <- vegan::mantel.partial(D$env, mk(fp$alpha * zg + task$beta * ze + fp$sigma * zn),
                                D$gen, permutations = SIM_PERM)
    c(m$statistic, m$signif)
  }))
  data.table(mg = task$mg, n = n, beta = task$beta,
             partial_r_mean = mean(o[, 1]), partial_r_sd = sd(o[, 1]),
             power = mean(o[, 2] < 0.05), n_sim = B_SIM)
}
tasks <- unlist(lapply(names(DISTS), function(mg)
  lapply(BETAS, function(b) list(mg = mg, beta = b))), recursive = FALSE)

if (N_CORES > 1) {
  cl <- makeCluster(N_CORES)
  clusterSetRNGStream(cl, SEED)
  clusterExport(cl, c("DISTS", "sim_fits", "zv", "B_SIM", "SIM_PERM", "N_NOISE"),
                envir = environment())
  clusterEvalQ(cl, { library(data.table); library(vegan) })
  pw <- rbindlist(parLapply(cl, tasks, one_cell))
  stopCluster(cl)
} else {
  rng_old <- RNGkind()[1]
  RNGkind("L'Ecuyer-CMRG"); set.seed(SEED)
  pw <- rbindlist(lapply(tasks, one_cell))
  RNGkind(rng_old)
}
pw[, power_se := sqrt(power * (1 - power) / n_sim)]
fwrite(pw, file.path(OUT_PW, "Table_PartialMantel_Power.csv"))

## ---- 4. 80% thresholds -------------------------------------------------------
thr <- rbindlist(lapply(split(pw, pw$mg), function(d) {
  d <- d[order(beta)]
  i <- which(d$power >= 0.8)[1]
  interp <- NA_real_
  if (!is.na(i) && i > 1) {
    a <- d[i - 1]; b <- d[i]
    w <- (0.8 - a$power) / (b$power - a$power)
    interp <- a$partial_r_mean + w * (b$partial_r_mean - a$partial_r_mean)
  }
  data.table(mg = d$mg[1], n = d$n[1], type_I = d[beta == 0, power],
             beta_at_80 = if (is.na(i)) NA_real_ else d$beta[i],
             partial_r_first_80 = if (is.na(i)) NA_real_ else d$partial_r_mean[i],
             power_there = if (is.na(i)) NA_real_ else d$power[i],
             partial_r_interp_80 = interp)
}))
thr <- merge(thr, sim_fits[, .(mg, partial_obs)], by = "mg")
fwrite(thr, file.path(OUT_PW, "Table_PartialMantel_Power80.csv"))
cat("\n-- power curve --\n"); print(dcast(pw, beta ~ mg, value.var = "power"))
cat("\n-- 80% thresholds --\n"); print(thr)
cat(sprintf("Done in %.1f min\n", as.numeric(difftime(Sys.time(), t_start, units = "mins"))))
