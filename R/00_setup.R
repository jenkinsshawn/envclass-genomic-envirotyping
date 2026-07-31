## Paths for the whole analysis. Set ROOT once; nothing else needs editing.

ROOT  <- normalizePath(getwd(), winslash = '/')
DATA  <- file.path(ROOT, 'data')
OUT   <- file.path(ROOT, 'out')
FIG   <- file.path(ROOT, 'figures')
CACHE <- file.path(DATA, 'nass_cache')

for (d in c(OUT, FIG)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

set.seed(42)
