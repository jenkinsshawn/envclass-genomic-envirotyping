if (!exists('ROOT')) source('R/00_setup.R')

## ============================================================================
## 17_covariates_by_type.R
## For each maturity group, the share of variation in each county covariate of
## the county model (R/12, 87 covariates, one county per testing location)
## that lies between the environment types: R-squared of a one-way analysis of
## variance of the covariate on type, with its P value and rank within the
## maturity group. Used to choose the two county backgrounds of Fig. 1
## (growing degree days and vapor pressure deficit, May to September).
## Builds the training counties exactly as R/12 does (sections 1 to 3).
## Writes out/Table_Covariates_by_Type.csv and out/Table_Covariates_Latitude.csv.
## ============================================================================

suppressPackageStartupMessages(library(data.table))
s12 <- readLines(file.path(ROOT, "R", "12_county_projection_v2.R"), warn = FALSE)
i0 <- grep("^suppressPackageStartupMessages\\(\\{", s12)[1]
i1 <- grep("^## ---- 4\\. Cross-validation", s12)[1] - 1
stopifnot(!is.na(i0), !is.na(i1))
s <- s12[i0:i1]
s <- s[!grepl("library\\(randomForest\\)", s)]
e <- new.env(parent = globalenv())
invisible(capture.output(eval(parse(text = s), envir = e)))

res <- rbindlist(lapply(names(e$train_sets), function(mg) {
  ts <- e$train_sets[[mg]]
  rbindlist(lapply(ts$feat, function(v) {
    y <- ts$X[, v]
    if (var(y) == 0) return(NULL)
    f <- lm(y ~ ts$y)
    data.table(mg = mg, covariate = v, counties = length(y), types = nlevels(ts$y),
               share_between_types = summary(f)$r.squared, p = anova(f)[1, 5])
  }))
}))
res[, rank := frank(-share_between_types, ties.method = "min"), by = mg]
setorder(res, mg, rank)
fwrite(res, file.path(OUT, "Table_Covariates_by_Type.csv"))
print(res[covariate %in% c("GDD_mean", "VPD_mean")])
mins <- res[, .(min_share = min(share_between_types)), by = covariate][order(-min_share)]
print(head(mins, 5))

## correlation across all counties of growing degree days and vapor pressure
## deficit with county centroid latitude (Fig. 1 legend and Results)
cl <- merge(unique(e$features[, .(fips, GDD_mean, VPD_mean)], by = "fips"),
            e$centroids[, .(fips, lat)], by = "fips")
lat_cor <- data.table(covariate = c("GDD_mean", "VPD_mean"), counties = nrow(cl),
                      r_with_latitude = c(cor(cl$GDD_mean, cl$lat, use = "complete.obs"),
                                          cor(cl$VPD_mean, cl$lat, use = "complete.obs")))
fwrite(lat_cor, file.path(OUT, "Table_Covariates_Latitude.csv"))
print(lat_cor)
