suppressPackageStartupMessages({ library(data.table); library(org.Hs.eg.db); library(AnnotationDbi) })
S <- read.csv("results/meta/meta_robust_DEGs.csv", stringsAsFactors = FALSE)
Tf <- read.csv("results/meta/robust_tcga_features.csv", stringsAsFactors = FALSE)
Rn <- read.csv("results/rnaseq/rnaseq_limma_voom_3contrasts.csv", stringsAsFactors = FALSE)[, c("gene", "lfc_T_GTEx", "fdr_T_GTEx")]
old <- read.csv("results/meta/meta_three_datasets_DEGs.csv", stringsAsFactors = FALSE)$gene
lit <- read.csv("results/meta/meta_DEGs_pubmed_counts_alias.csv", stringsAsFactors = FALSE)[, c("gene", "pdac_best")]
ge <- as.data.frame(fread("depmap/CRISPRGeneEffect.csv")); ids <- ge[[1]]; ge <- ge[, -1]; rownames(ge) <- ids
colnames(ge) <- sub(" [(].*$", "", colnames(ge)); ge <- ge[, !duplicated(colnames(ge))]
md <- read.csv("depmap/Model.csv", check.names = FALSE, stringsAsFactors = FALSE)
pdac <- intersect(md$ModelID[md$OncotreePrimaryDisease %in% "Pancreatic Adenocarcinoma"], rownames(ge)); oth <- setdiff(rownames(ge), pdac)
Gd <- intersect(Tf$gene, colnames(ge))
dm <- do.call(rbind, lapply(Gd, function(g) { a <- ge[pdac, g]; b <- ge[oth, g]
  data.frame(gene = g, med_pdac = median(a, na.rm = TRUE), med_other = median(b, na.rm = TRUE), frac_dep_pdac = mean(a < -0.5, na.rm = TRUE),
             frac_dep_other = mean(b < -0.5, na.rm = TRUE), p_dm = tryCatch(wilcox.test(a, b)$p.value, error = function(e) NA)) }))
dm$delta <- dm$med_pdac - dm$med_other; dm$fdr_dm <- p.adjust(dm$p_dm, "BH")
dm$class <- ifelse(dm$frac_dep_other >= 0.9, "common essential", ifelse(dm$delta <= -0.2 & dm$fdr_dm < 0.05, "PDAC-selective",
            ifelse(dm$frac_dep_pdac >= 0.5, "dependent in PDAC", "not a dependency")))
go <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = c("GO:0006325", "GO:0016569", "GO:0006306", "GO:0080111", "GO:0042393", "GO:0042054"), keytype = "GOALL", columns = "SYMBOL"))
D <- merge(S, Tf, by = "gene"); D <- merge(D, Rn, by = "gene", all.x = TRUE); D <- merge(D, dm[, c("gene", "med_pdac", "med_other", "frac_dep_pdac", "delta", "class")], by = "gene", all.x = TRUE)
D <- merge(D, lit, by = "gene", all.x = TRUE)
D$new <- !(D$gene %in% old); D$dir <- ifelse(D$lfc_core > 0, "up", "down"); D$epigenetic <- D$gene %in% unique(go$SYMBOL)
D$rna_ok <- !is.na(D$lfc_T_GTEx) & sign(D$lfc_T_GTEx) == sign(D$lfc_core) & D$fdr_T_GTEx < 0.05
D <- D[order(-abs(D$lfc_core)), ]; write.csv(D, "results/meta/robust_funnel_table.csv", row.names = FALSE)
cat("Genes:", nrow(D), "| new vs earlier 402:", sum(D$new), "| in DepMap:", sum(!is.na(D$class)), "\n")
cat("\nDepMap classes:\n"); print(table(D$class, useNA = "ifany")); cat("\nRNA-seq agrees (direction, FDR<0.05):", sum(D$rna_ok), "of", nrow(D), "\n")
cols <- c("gene", "new", "lfc_core", "I2_core", "comp", "r_prolif", "r_classical", "r_basal", "HR", "p", "p_pro", "class", "pdac_best", "epigenetic")
cand <- D[D$dir == "up" & D$comp %in% c("epithelial", "unassigned") & D$r_prolif < 0.5 & D$rna_ok, ]
cat("\n1) Tumor-up, tumor-cell-like, not proliferation-linked, RNA-seq confirmed:", nrow(cand), "| new genes among them:", sum(cand$new), "\nTop 30 NEW by effect size:\n"); print(head(cand[cand$new, cols], 30), row.names = FALSE, digits = 3)
ind <- D[D$dir == "up" & !is.na(D$p) & D$p < 0.05 & D$p_pro < 0.05, ]
cat("\n2) Prognostic (TCGA p<0.05) and independent of proliferation (p<0.05):", nrow(ind), "\n"); print(head(ind[order(ind$p), cols], 25), row.names = FALSE, digits = 3)
dep <- D[D$class %in% c("PDAC-selective", "dependent in PDAC") & D$r_prolif < 0.5, ]
cat("\n3) DepMap dependent in PDAC and not proliferation-like:", nrow(dep), "\n"); print(head(dep[order(dep$med_pdac), c(cols[1:9], "med_pdac", "delta", "class")], 25), row.names = FALSE, digits = 3)
cat("\n4) Epigenetic-annotated among list 1:", sum(cand$epigenetic), "\n"); print(cand[cand$epigenetic, cols], row.names = FALSE, digits = 3)
