# run_permnull_power_local.ps1
# Permutation test of cluster separation (R\10_permutation_null_v2.R) and
# power of the partial Mantel test (R\11_partial_mantel_power.R).
# Both full runs were completed in the analysis session (about 6 minutes and
# under 1 minute);
# this file is for rerunning them on the local R installation.
#
# Folder: the repository root, with R\10_permutation_null_v2.R and
# R\11_partial_mantel_power.R added (and R\02_phenotypic_mantel.R unchanged).

$Repo = "C:\Users\sje66200\Downloads\envclass-genomic-envirotyping"
$R    = "C:\Users\sje66200\AppData\Local\Programs\R\R-4.6.1\bin\Rscript.exe"

Set-Location $Repo
New-Item -ItemType Directory -Force logs | Out-Null

# Step 0 (once): packages. 11 runs R\02_phenotypic_mantel.R, which needs
# sf, tigris (downloads the 2022 county boundaries on first use) and ecodist.
& $R -e "for (p in c('data.table','cluster','vegan','sf','tigris','ecodist')) if (!requireNamespace(p, quietly=TRUE)) install.packages(p, repos='https://cloud.r-project.org')"

# Step 1 (about 6 minutes): permutation test, 1,000 permutations, profiles as
# in 01_marker_profiles.R, optimal silhouette over k = 2 to 8 for each
# permutation. Also writes figures\FigS1_Permutation_Null.tiff.
$env:ENVCLASS_PERM_MODE = "classification"; $env:ENVCLASS_NULL_K = "optimal"; $env:ENVCLASS_N_PERM = "1000"
& $R R\10_permutation_null_v2.R *> logs\10_permnull_classification.log

# Step 1b (optional, under 1 minute): the same with 100 permutations, the count
# in the current draft. Overwrites the step 1 tables and figure.
# $env:ENVCLASS_N_PERM = "100"; & $R R\10_permutation_null_v2.R *> logs\10_permnull_100.log

# Step 1c (optional): the earlier constructions, for comparison only.
# $env:ENVCLASS_N_PERM = "100"; $env:ENVCLASS_NULL_K = "observed"
# & $R R\10_permutation_null_v2.R *> logs\10_permnull_kobs.log          # reproduces perm_null_v2.rds
# $env:ENVCLASS_NULL_K = "optimal"; $env:ENVCLASS_PERM_MODE = "pooled_adjdev"
# & $R R\10_permutation_null_v2.R *> logs\10_permnull_pooled_adjdev.log # 07 make_permnull
Remove-Item Env:ENVCLASS_PERM_MODE, Env:ENVCLASS_NULL_K, Env:ENVCLASS_N_PERM -ErrorAction SilentlyContinue

# Step 2 (test, under 1 minute): power simulation with 10 simulations per effect size.
$env:ENVCLASS_B_SIM = "10"; $env:ENVCLASS_N_CORES = "8"
& $R R\11_partial_mantel_power.R *> logs\11_power_test.log

# Step 3 (about 4 minutes on 8 cores): full run, 300 simulations per effect
# size at 199 permutations. Reuses out\power\table1_distances.rds from step 2.
$env:ENVCLASS_B_SIM = "300"; $env:ENVCLASS_N_CORES = "8"
& $R R\11_partial_mantel_power.R *> logs\11_power_full.log
Remove-Item Env:ENVCLASS_B_SIM, Env:ENVCLASS_N_CORES -ErrorAction SilentlyContinue

# Results: out\permutation_null\Table_PermNull_classification.csv,
# out\power\Table_PartialMantel_Power.csv, Table_PartialMantel_Power80.csv,
# Table_PartialMantel_Table1Check.csv (rebuilt n and Mantel r against Table 1).
