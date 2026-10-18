if (!requireNamespace("data.table", quietly = TRUE)) install.packages("data.table", repos = "https://cloud.r-project.org")
ge <- as.data.frame(data.table::fread("depmap/CRISPRGeneEffect.csv"))
ids <- ge[[1]]; ge <- ge[, -1]; rownames(ge) <- ids
colnames(ge) <- sub(" [(].*$", "", colnames(ge)); ge <- ge[, !duplicated(colnames(ge))]
md <- read.csv("depmap/Model.csv", check.names = FALSE, stringsAsFactors = FALSE)
if (!all(c("ModelID", "OncotreePrimaryDisease") %in% names(md))) stop("Model.csv columns: ", paste(names(md), collapse = ", "))
pdac <- intersect(md$ModelID[md$OncotreePrimaryDisease %in% "Pancreatic Adenocarcinoma"], rownames(ge))
if (length(pdac) == 0) { print(table(md$OncotreePrimaryDisease[md$OncotreeLineage %in% "Pancreas"])); stop("No PDAC lines found, check the disease label above") }
oth <- setdiff(rownames(ge), pdac)
cat("Models screened:", nrow(ge), "| PDAC lines:", length(pdac), "\n")
S <- read.csv("results/meta/gene_prioritization.csv", stringsAsFactors = FALSE)
G <- intersect(S$gene, colnames(ge))
cat("DEGs present in the CRISPR data:", length(G), "of", nrow(S), "\n")
res <- do.call(rbind, lapply(G, function(g) {
  a <- ge[pdac, g]; b <- ge[oth, g]
  data.frame(gene = g, med_pdac = median(a, na.rm = TRUE), med_other = median(b, na.rm = TRUE),
             frac_dep_pdac = mean(a < -0.5, na.rm = TRUE), frac_dep_other = mean(b < -0.5, na.rm = TRUE),
             p = tryCatch(wilcox.test(a, b)$p.value, error = function(e) NA))
}))
res$delta <- res$med_pdac - res$med_other
res$fdr <- p.adjust(res$p, "BH")
res$class <- ifelse(res$frac_dep_other >= 0.9, "common essential",
             ifelse(res$delta <= -0.2 & res$fdr < 0.05, "PDAC-selective",
             ifelse(res$frac_dep_pdac >= 0.5, "dependent in PDAC (not selective)", "not a dependency")))
D <- merge(S[, c("gene", "lfc", "I2", "comp", "pdac_best", "surv_n", "score")], res, by = "gene")
D <- D[order(D$med_pdac), ]
write.csv(D, "results/meta/depmap_dependency_DEGs.csv", row.names = FALSE)
cat("\nClass counts:\n"); print(table(D$class))
sl <- c("LAMC2","ITGA2","SERPINB5","ITGB6","ANLN","MELK","FGD6","ITGA3","MET","CDH3","TOP2A","TPX2","DLGAP5","DCBLD2","NT5E","VSIG1","SCEL","FER1L6","TMC5","DPCR1","SYTL2","PLA2R1","ANXA8L1","CEACAM5","CEACAM6","SLC6A14","CTSE","TMPRSS4","ANXA10")
cols <- c("gene", "lfc", "comp", "pdac_best", "med_pdac", "med_other", "frac_dep_pdac", "class")
cat("\nShortlist:\n"); print(D[D$gene %in% sl, cols], row.names = FALSE, digits = 2)
cat("\nShortlist genes not screened in DepMap:", paste(setdiff(sl, D$gene), collapse = ", "), "\n")
cat("\nTop PDAC-selective dependencies among all DEGs:\n")
print(head(D[D$class == "PDAC-selective", c(cols, "fdr")][order(D$delta[D$class == "PDAC-selective"]), ], 20), row.names = FALSE, digits = 2)
