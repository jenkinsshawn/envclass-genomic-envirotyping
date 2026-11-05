if (!exists('ROOT')) source('R/00_setup.R')

## ============================================================================
## 14_raw_yield_response.R
## Phenotypic response distance computed from raw line means instead of
## adjusted deviations (Methods, Phenotypic response distance and split-half
## control; Online Resource 1, Supplementary Methods).
##
## Runs R/02_phenotypic_mantel.R up to the distance matrices for MG I to IV,
## on the same locations as Table 2, and adds for each maturity group:
##   - phenotypic response distance from raw line-by-location mean yields
##     (1 minus the correlation over lines tested at both locations, at least
##     8 lines), and its Mantel correlation with genomic profile distance;
##   - the mean correlation between locations of raw yields and of adjusted
##     deviations over shared lines;
##   - the Mantel correlation between the raw and the adjusted response
##     distance matrices.
##
## ENVCLASS_LOC_FIPS (optional): path to data/loc_fips_offline.csv, used in
## place of the sf/tigris county join, as in R/11.
## Writes out/Table_RawYield_Response.csv. About 1 minute.
## ============================================================================

suppressPackageStartupMessages({ library(data.table); library(vegan) })
N_PERM <- 9999
OUT_RAW <- Sys.getenv("ENVCLASS_OUT", OUT)

src02 <- readLines(file.path(ROOT, "R", "02_phenotypic_mantel.R"), warn = FALSE)
cut   <- grep("^fwrite", src02)[1] - 1
src_run <- src02[1:cut]

## raw line means at each location, built next to the adjusted deviations
a1 <- grep("d_pheno_full <- pairwise_cor_dist\\(m_adj", src_run)[1]
stopifnot(!is.na(a1))
src_run <- append(src_run, c(
  "  m_raw <- build_strain_loc_mat(",
  "    mg_means[loc_key %in% loc_in_prof, .(Strain, loc_key, value = YieldBuA)])",
  "  d_raw_full <- pairwise_cor_dist(m_raw, PARAMS$min_pair_strains)"),
  after = a1)

## after the Table 2 distance matrices exist, compute the raw-yield summaries
## and return; the remaining Table 2 tests are not repeated here
a2 <- grep("d_geo\\s*<- clean_dist", src_run)[1]
stopifnot(!is.na(a2))
src_run <- append(src_run, c(
  "  d_raw <- clean_dist(stats::as.dist(d_raw_full[shared, shared, drop = FALSE]))",
  "  stopifnot(!anyNA(d_raw))",
  sprintf("  mt_adj <- vegan::mantel(d_gen, d_pheno, permutations = %d)", N_PERM),
  sprintf("  mt_raw <- vegan::mantel(d_gen, d_raw, permutations = %d)", N_PERM),
  sprintf("  mt_rw  <- vegan::mantel(d_raw, d_pheno, permutations = %d)", N_PERM),
  "  return(data.table(mg = mg_val, n_locs = length(shared),",
  "    mean_cor_raw = mean(1 - as.vector(d_raw)), mean_cor_adj = mean(1 - as.vector(d_pheno)),",
  "    mantel_genomic_adj = mt_adj$statistic, mantel_genomic_adj_p = mt_adj$signif,",
  "    mantel_genomic_raw = mt_raw$statistic, mantel_genomic_raw_p = mt_raw$signif,",
  "    mantel_raw_vs_adj = mt_rw$statistic, mantel_raw_vs_adj_p = mt_rw$signif))"),
  after = a2)

## MG I to IV only, as in the distance analyses
i_mg <- grep("mgs\\s*=\\s*c\\(\"0\"", src_run)[1]
stopifnot(!is.na(i_mg))
src_run[i_mg] <- sub('c\\("0", ', 'c(', src_run[i_mg])

LOC_FIPS <- Sys.getenv("ENVCLASS_LOC_FIPS", "")
if (nzchar(LOC_FIPS)) {
  i0 <- grep("^suppressPackageStartupMessages\\(\\{", src_run)[1]
  i1 <- grep('^for \\(pkg in c\\("tigris", "sf"', src_run)[1]
  stopifnot(!is.na(i0), !is.na(i1), i0 < i1)
  src_run <- c(src_run[1:(i0 - 1)],
    sprintf('loc_fips <- fread("%s", colClasses = list(character = "fips"))[loc_key %%in%% loc_coords[lat > 24 & lat < 50 & lon > -130 & lon < -60, loc_key], .(loc_key, fips)]', LOC_FIPS),
    src_run[(i1 + 1):length(src_run)])
}

tmp <- tempfile(fileext = ".R"); writeLines(src_run, tmp)
e02 <- new.env(parent = globalenv())
source(tmp, local = e02, echo = FALSE)
res <- rbindlist(e02$results, fill = TRUE)
print(res)
fwrite(res, file.path(OUT_RAW, "Table_RawYield_Response.csv"))
