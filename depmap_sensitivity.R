if (!requireNamespace("data.table", quietly = TRUE)) install.packages("data.table", repos = "https://cloud.r-project.org")
ge <- as.data.frame(data.table::fread("depmap/CRISPRGeneEffect.csv"))
ids <- ge[[1]]; ge <- ge[, -1]; rownames(ge) <- ids
colnames(ge) <- sub(" [(].*$", "", colnames(ge)); ge <- ge[, !duplicated(colnames(ge))]
md <- read.csv("depmap/Model.csv", check.names = FALSE, stringsAsFactors = FALSE)
pdac <- intersect(md$ModelID[md$OncotreePrimaryDisease %in% "Pancreatic Adenocarcinoma"], rownames(ge))
oth <- setdiff(rownames(ge), pdac)
lin <- md$OncotreeLineage[match(rownames(ge), md$ModelID)]; names(lin) <- rownames(ge)
cat("PDAC lines:", length(pdac), "| other lines:", length(oth), "\nLineages among other lines:\n")
print(sort(table(lin[oth]), decreasing = TRUE))
nonepi <- c("Myeloid", "Lymphoid", "CNS/Brain", "Peripheral Nervous System", "Soft Tissue", "Bone", "Skin", "Eye", "Other", "Fibroblast")
epi <- oth[!(lin[oth] %in% nonepi)]
cat("\nEpithelial-type comparison lines:", length(epi), "(excluded lineages:", paste(nonepi, collapse = ", "), ")\n")

S <- read.csv("results/meta/robust_tcga_features.csv", stringsAsFactors = FALSE)
G <- intersect(S$gene, colnames(ge)); cat("Genes tested:", length(G), "\n")
cmp <- function(a, b) {
  a <- a[!is.na(a)]; b <- b[!is.na(b)]
  w <- suppressWarnings(wilcox.test(a, b))
  c(delta = median(a) - median(b), cliff = 1 - 2 * unname(w$statistic) / (length(a) * length(b)),
    p = w$p.value, frac_other = mean(b < -0.5))
}
res <- do.call(rbind, lapply(G, function(g) {
  x <- cmp(ge[pdac, g], ge[oth, g]); y <- cmp(ge[pdac, g], ge[epi, g])
  data.frame(gene = g, med_pdac = median(ge[pdac, g], na.rm = TRUE),
             delta_all = x["delta"], cliff_all = x["cliff"], p_all = x["p"], frac_other_all = x["frac_other"],
             delta_epi = y["delta"], cliff_epi = y["cliff"], p_epi = y["p"], frac_other_epi = y["frac_other"])
}))
rownames(res) <- NULL
res$fdr_all <- p.adjust(res$p_all, "BH"); res$fdr_epi <- p.adjust(res$p_epi, "BH")
res$rule_all    <- res$frac_other_all < 0.9 & res$delta_all <= -0.2 & res$fdr_all < 0.05
res$rule_epi    <- res$frac_other_epi < 0.9 & res$delta_epi <= -0.2 & res$fdr_epi < 0.05
res$rule_strict <- res$rule_all & res$rule_epi & res$delta_all <= -0.3 & res$delta_epi <= -0.3
cat("\nPDAC-selective, rule as used (all other lines):", sum(res$rule_all), "\n")
cat("PDAC-selective vs epithelial-type lines only:", sum(res$rule_epi), "\n")
cat("Strict (both comparisons, delta <= -0.3):", sum(res$rule_strict), "\n")
cat("\nGenes passing the rule as used:", paste(res$gene[res$rule_all], collapse = ", "), "\n")
cat("Genes passing vs epithelial lines:", paste(res$gene[res$rule_epi], collapse = ", "), "\n")
cat("Genes passing the strict rule:", paste(res$gene[res$rule_strict], collapse = ", "), "\n")
old <- read.csv("results/meta/robust_funnel_table.csv", stringsAsFactors = FALSE)
cat("\nConsistency check, 'PDAC-selective' in saved table:", sum(old$class == "PDAC-selective"), "\n")
num <- sapply(res, is.numeric); res[num] <- lapply(res[num], signif, 3)
write.csv(res, "results/meta/depmap_dependency_sensitivity.csv", row.names = FALSE)
cat("\nKLF5 and top rows:\n"); print(head(res[order(res$p_all), c("gene","med_pdac","delta_all","fdr_all","delta_epi","fdr_epi","cliff_all","rule_all","rule_epi","rule_strict")], 12), row.names = FALSE)
if ("KLF5" %in% res$gene) print(res[res$gene == "KLF5", ], row.names = FALSE)
