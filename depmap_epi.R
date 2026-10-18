suppressPackageStartupMessages(library(data.table))
ge <- as.data.frame(fread("depmap/CRISPRGeneEffect.csv"))
ids <- ge[[1]]; ge <- ge[, -1]; rownames(ge) <- ids
colnames(ge) <- sub(" [(].*$", "", colnames(ge)); ge <- ge[, !duplicated(colnames(ge))]
md <- read.csv("depmap/Model.csv", check.names = FALSE, stringsAsFactors = FALSE)
pdac <- intersect(md$ModelID[md$OncotreePrimaryDisease %in% "Pancreatic Adenocarcinoma"], rownames(ge))
oth <- setdiff(rownames(ge), pdac)
E <- read.csv("results/meta/epigenetic_regulators_relaxed.csv", stringsAsFactors = FALSE)
G <- intersect(E$gene, colnames(ge))
cat("Epigenetic-annotated genes:", nrow(E), "| in CRISPR data:", length(G), "| not screened:", paste(setdiff(E$gene, G), collapse = ", "), "\n")
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
             ifelse(res$frac_dep_pdac >= 0.5, "dependent in PDAC", "not a dependency")))
D <- merge(E[, c("gene", "lfc", "I2", "pdac")], res, by = "gene")
if (file.exists("results/survival/cox_TCGA_PAAD.csv")) {
  cx <- read.csv("results/survival/cox_TCGA_PAAD.csv", stringsAsFactors = FALSE)[, c("gene", "HR_adj", "p_adj")]
  D <- merge(D, cx, by = "gene", all.x = TRUE)
}
D <- D[order(D$med_pdac), ]
write.csv(D, "results/meta/epigenetic_depmap.csv", row.names = FALSE)
print(D[, setdiff(names(D), c("p", "frac_dep_other"))], row.names = FALSE, digits = 2)
