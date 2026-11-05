if (!exists('ROOT')) source('R/00_setup.R')

## ============================================================================
## 25_table_S5_check.R  (Online Resource 1, Table S2 and annotation text)
## Checks the recovered annotation tables against the archived marker
## results and the Online Resource text, and applies the stated Gene Ontology
## term-size rule.
##
## Inputs (recovered from the author's computer, not in the v1.0.0 archive;
## place them in data/annotation/):
##   Table_S_Discriminating_Genes_a2.csv : top 30 Kruskal-Wallis markers per
##       maturity group, nearest Wm82.a2.v1 gene, distance, Arabidopsis hit
##   Table_GO_enrichment_a2.csv          : hypergeometric GO test per group,
##       produced by the recovered console script of 2026-07-16 19:46
##       (BH over terms with >= 1 foreground gene; no term-size filter)
##
## Output: out/Table_S5.csv (formatted Table S2), out/Table_S5_summary.csv,
##         out/Table_S5_GO_min10.csv
## ============================================================================

library(data.table)
if (!exists('ANN')) ANN <- file.path(DATA, 'annotation')
mv2 <- readRDS(file.path(DATA, 'marker_v2_results.rds'))
g  <- fread(file.path(ANN, 'Table_S_Discriminating_Genes_a2.csv'))
go <- fread(file.path(ANN, 'Table_GO_enrichment_a2.csv'))
g[, mg := as.character(maturity_group)]; go[, mg := as.character(mg)]
MGS <- c('0', 'I', 'II', 'III', 'IV')

## top 30 markers in the table are the top 30 of the archived Kruskal-Wallis
## results (data/marker_v2_results.rds, Method F)?
top_match <- sapply(MGS, function(m) {
  kw <- as.data.table(mv2[[m]]$F$kw)[order(kw_p)]
  setequal(as.character(head(kw$marker, 30)), as.character(g[mg == m, marker]))
})

summ <- g[, .(n_markers = .N, within_50kb = sum(within_50kb),
              median_distance_bp_within50 = as.numeric(median(distance_bp[within_50kb])),
              median_distance_bp_all = as.numeric(median(distance_bp)),
              arabidopsis_hit = sum(arabidopsis_hit != ''),
              arabidopsis_named_symbol = sum(arabidopsis_symbol != ''),
              min_BH_p = min(BH_adjusted_p)), by = mg]
summ[, top30_matches_archive := top_match[mg]]
summ[, n_discriminating_BH05 := sapply(mg, function(m) mv2[[m]]$F$n_sig)]

## Gene Ontology: apply the Online Resource rule (terms with >= 10 annotated
## background genes) and re-adjust within maturity group. Terms absent from
## the table had no foreground gene (p = 1); adding them would only make the
## adjustment more conservative.
go_all <- go[, .(mg, go, k, K, p, p_adj_as_archived = p_adj)]
go10 <- go_all[K >= 10][, p_adj_min10 := p.adjust(p, 'BH'), by = mg]
go_summ <- merge(
  go_all[, .(terms_archived = .N, sig_archived = sum(p_adj_as_archived < 0.05),
             min_padj_archived = min(p_adj_as_archived)), by = mg],
  go10[, .(terms_min10 = .N, sig_min10 = sum(p_adj_min10 < 0.05),
           min_padj_min10 = min(p_adj_min10)), by = mg], by = 'mg')
summ <- merge(summ, go_summ, by = 'mg')[match(MGS, mg)]

## multi-marker groups: >= 4 top-30 markers on one chromosome, consecutive
## markers no more than 0.5 Mb apart (reproduces the three groups in the text)
runs <- g[order(mg, chromosome, position_bp)][, {
  grp <- cumsum(c(1, diff(position_bp) > 5e5))
  .SD[, .(n = .N, start_mb = min(position_bp) / 1e6, end_mb = max(position_bp) / 1e6),
      by = .(grp = grp)]
}, by = .(mg, chromosome)][n >= 4]

s5 <- g[, .(Maturity_group = mg, Marker = marker, Chromosome = chromosome,
            Position_bp_Wm82a2 = position_bp, KW_p = signif(kruskal_wallis_p, 3),
            BH_adjusted_p = signif(BH_adjusted_p, 3), Nearest_gene = nearest_gene,
            Distance_bp = distance_bp, Within_50kb = within_50kb,
            Arabidopsis_ortholog = arabidopsis_hit, Symbol = arabidopsis_symbol,
            Description = arabidopsis_defline, GO_terms = GO_terms)]
fwrite(s5, file.path(OUT, 'Table_S5.csv'))
fwrite(summ, file.path(OUT, 'Table_S5_summary.csv'))
fwrite(go10, file.path(OUT, 'Table_S5_GO_min10.csv'))
fwrite(runs, file.path(OUT, 'Table_S5_marker_runs.csv'))
print(summ); print(runs)
