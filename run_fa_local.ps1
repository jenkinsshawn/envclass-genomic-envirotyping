# run_fa_local.ps1
# Factor-analytic grouping of locations (maturity groups II and III) for the
# EnvClass shared-lines comparison. Run each step in a PowerShell window.
#
# Folder: C:\Users\sje66200\Downloads\envclass-genomic-envirotyping
# (repository v1.0.0 with R\08_shared_lines_comparison.R and R\09_fa_loadings_local.R added).
# Step 1 has already been run here; its tables are in out\shared_lines\.

$Repo = "C:\Users\sje66200\Downloads\envclass-genomic-envirotyping"
$R    = "C:\Users\sje66200\AppData\Local\Programs\R\R-4.6.1\bin\Rscript.exe"

Set-Location $Repo
New-Item -ItemType Directory -Force logs | Out-Null

# Step 0 (a few minutes, once): install the packages R 4.6.1 does not have yet.
& $R -e "for (p in c('data.table','cluster','sommer')) if (!requireNamespace(p, quietly=TRUE)) install.packages(p, repos='https://cloud.r-project.org')"

# Step 1 (under 1 minute; optional, already done): location sets and shared-line tables.
# & $R R\08_shared_lines_comparison.R *> logs\08_shared_lines.log

# Step 2 (one fit, to time it): MG II, validation locations, 2 factors.
& $R R\09_fa_loadings_local.R --mg=II --sets=validation --factors=2 *> logs\09_fa_test.log
Get-Content out\shared_lines\fa\fa_run_log.txt

# Step 3 (all 8 fits: MG II and III, both location sets, 2 and 3 factors).
# Run after checking the time for step 2; expect about 8 times that.
& $R R\09_fa_loadings_local.R *> logs\09_fa_full.log

# Results: out\shared_lines\fa\Table_FA_ARI.csv, Table_FA_Loadings.csv,
# fa_run_log.txt and one *_fit.rds per fitted model.
