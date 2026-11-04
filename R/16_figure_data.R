if (!exists('ROOT')) source('R/00_setup.R')

## ============================================================================
## 16_figure_data.R
## Writes the small tables that figures_src/make_figures.py reads, to
## out/figure_data/:
##   pairs.csv        location pairs of the Table 1 analysis (MG I to IV) with
##                    genomic profile, phenotypic response, envirotyping and
##                    geographic distance (rebuilt from R/02, same locations)
##   mds.csv          two-dimensional scaling of profile distance, all MGs
##   loc_coords.csv   median latitude and longitude of each location
##   county_gdd.csv   county mean growing degree days, May to September
## ENVCLASS_LOC_FIPS (optional): data/loc_fips_offline.csv in place of the
## sf/tigris county join, as in R/11 and R/14.
## ============================================================================

suppressPackageStartupMessages(library(data.table))
FD <- file.path(OUT, "figure_data"); dir.create(FD, showWarnings = FALSE, recursive = TRUE)

src02 <- readLines(file.path(ROOT, "R", "02_phenotypic_mantel.R"), warn = FALSE)
cut <- grep("^fwrite", src02)[1] - 1; s <- src02[1:cut]
a2 <- grep("d_geo\\s*<- clean_dist", s)[1]; stopifnot(!is.na(a2))
s <- append(s, c(
  "  lab <- attr(d_gen, 'Labels'); ij <- which(lower.tri(matrix(0, length(lab), length(lab))), arr.ind = TRUE)",
  "  return(data.table(mg = mg_val, loc_a = lab[ij[,2]], loc_b = lab[ij[,1]], gen = as.vector(d_gen), pheno = as.vector(d_pheno), env = as.vector(d_env), geo = as.vector(d_geo)))"),
  after = a2)
i_mg <- grep("mgs\\s*=\\s*c\\(\"0\"", s)[1]; s[i_mg] <- sub('c\\("0", ', 'c(', s[i_mg])
LOC_FIPS <- Sys.getenv("ENVCLASS_LOC_FIPS", "")
if (nzchar(LOC_FIPS)) {
  i0 <- grep("^suppressPackageStartupMessages\\(\\{", s)[1]; i1 <- grep('^for \\(pkg in c\\("tigris", "sf"', s)[1]
  s <- c(s[1:(i0 - 1)], sprintf('loc_fips <- fread("%s", colClasses = list(character = "fips"))[loc_key %%in%% loc_coords[lat > 24 & lat < 50 & lon > -130 & lon < -60, loc_key], .(loc_key, fips)]', LOC_FIPS),
         s[(i1 + 1):length(s)])
}
tmp <- tempfile(fileext = ".R"); writeLines(s, tmp)
e <- new.env(parent = globalenv())
invisible(capture.output(source(tmp, local = e)))
pairs <- rbindlist(e$results)
fwrite(pairs, file.path(FD, "pairs.csv"))

v2 <- readRDS(file.path(DATA, "marker_v2_results.rds"))
mds <- rbindlist(lapply(c("0", "I", "II", "III", "IV"), function(m) {
  r <- v2[[m]]$F; co <- cmdscale(r$dist, k = 2)
  data.table(mg = m, location = rownames(co), MDS1 = co[, 1], MDS2 = co[, 2],
             type = r$clusters[rownames(co)], k = r$k, sil = r$sil)
}))
fwrite(mds, file.path(FD, "mds.csv"))

nust <- fread(file.path(DATA, "nust_ecw_subset.csv"), select = c("Location", "Lat", "Long"))
nust[, loc_key := trimws(tolower(gsub(",.*$", "", Location)))]
fwrite(nust[, .(lat = median(Lat, na.rm = TRUE), lon = median(Long, na.rm = TRUE)), by = loc_key],
       file.path(FD, "loc_coords.csv"))

pw <- as.data.table(readRDS(file.path(CACHE, "nasa_power_counties.rds")))
setnames(pw, grep("(?i)^fips|geoid", names(pw), value = TRUE)[1], "fips")
pw[, fips := sprintf("%05d", as.integer(fips))]
fwrite(unique(pw[, .(fips, GDD_mean)], by = "fips"), file.path(FD, "county_gdd.csv"))
cat("figure data written to", FD, "\n")
