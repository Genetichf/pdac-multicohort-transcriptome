suppressPackageStartupMessages({ library(SummarizedExperiment); library(recount3); library(edgeR); library(org.Hs.eg.db); library(AnnotationDbi) })
S <- read.csv("results/meta/meta_DEGs_pubmed_counts_alias.csv", stringsAsFactors = FALSE)
S$meta_fdr <- S$fdr
mk <- list(stroma = c("COL1A1","COL1A2","COL3A1","DCN","LUM","FAP"), immune = c("PTPRC","CD53","LAPTM5","CORO1A","CD3E","CD68"),
           epithelial = c("KRT19","KRT7","EPCAM","KRT8","KRT18","CLDN4"), acinar = c("PRSS1","CPA1","PNLIP","CTRB1"),
           classical = c("GATA6","TFF1","AGR2","LGALS4"), basal = c("KRT5","KRT6A","TP63","S100A2"))
rt <- readRDS("results/rnaseq/rse_tcga_PAAD.rds")
assay(rt, "counts") <- compute_read_counts(rt)
rt <- rt[, rt$tcga.cgc_sample_sample_type %in% "Primary Tumor" & rt$tcga.cgc_case_histological_diagnosis %in% "Pancreas-Adenocarcinoma Ductal Type"]
rd <- rowData(rt)
lcpm <- cpm(calcNormFactors(DGEList(assay(rt, "counts"))), log = TRUE, prior.count = 1)
need <- unique(c(S$gene, unlist(mk)))
pc <- rd$gene_type %in% "protein_coding" & rd$gene_name %in% need
lcpm <- lcpm[pc, ]; sym <- rd$gene_name[pc]
mu <- rowMeans(lcpm); best <- tapply(seq_along(sym), sym, function(i) i[which.max(mu[i])])
X <- lcpm[unlist(best), ]; rownames(X) <- names(best)
sc <- sapply(mk, function(g) { g <- intersect(g, rownames(X)); colMeans(t(scale(t(X[g, , drop = FALSE])))) })
G <- intersect(S$gene, rownames(X))
r <- t(sapply(G, function(g) as.numeric(cor(X[g, ], sc, method = "spearman"))))
colnames(r) <- paste0("r_", colnames(sc))
main <- c("stroma", "immune", "epithelial", "acinar")
rm4 <- r[, paste0("r_", main)]
comp <- ifelse(apply(rm4, 1, max) >= 0.4, main[apply(rm4, 1, which.max)], "unassigned")
rc <- data.frame(gene = rownames(r), round(r[, c("r_classical", "r_basal")], 2), comp = comp)
go <- suppressMessages(AnnotationDbi::select(org.Hs.eg.db, keys = c("GO:0006325", "GO:0016569", "GO:0006306"), keytype = "GOALL", columns = "SYMBOL"))
epi <- unique(go$SYMBOL)
R <- read.csv("results/rnaseq/rnaseq_limma_voom_3contrasts.csv", stringsAsFactors = FALSE)
Cx <- read.csv("results/survival/cox_GSE28735_GSE62452.csv", stringsAsFactors = FALSE)[, c("gene", "HR_A", "p_A", "HR_B", "p_B")]
Tc <- read.csv("results/survival/cox_TCGA_PAAD.csv", stringsAsFactors = FALSE)[, c("gene", "HR", "p")]; names(Tc) <- c("gene", "HR_TCGA", "p_TCGA")
D <- merge(S[, c("gene", "lfc", "I2", "k", "meta_fdr", "pdac_best")], R[, c("gene", "lfc_T_GTEx", "fdr_T_GTEx")], by = "gene", all.x = TRUE)
D <- merge(D, Cx, by = "gene", all.x = TRUE); D <- merge(D, Tc, by = "gene", all.x = TRUE); D <- merge(D, rc, by = "gene", all.x = TRUE)
D$epigenetic <- D$gene %in% epi
D$dir <- ifelse(D$lfc > 0, "up", "down")
D$rna_ok <- !is.na(D$lfc_T_GTEx) & sign(D$lfc_T_GTEx) == sign(D$lfc) & D$fdr_T_GTEx < 0.05
D$surv_n <- rowSums(cbind(D$p_A < 0.05, D$p_B < 0.05, D$p_TCGA < 0.05), na.rm = TRUE)
D$surv_dir <- !is.na(D$HR_A) & !is.na(D$HR_B) & !is.na(D$HR_TCGA) & sign(log(D$HR_A)) == sign(log(D$HR_B)) & sign(log(D$HR_B)) == sign(log(D$HR_TCGA))
D$s_effect <- (abs(D$lfc) >= 1.5) + (abs(D$lfc) >= 2)
D$s_consist <- ifelse(D$I2 < 50, 2, ifelse(D$I2 < 75, 1, 0))
D$s_rna <- as.numeric(D$rna_ok)
D$s_tumorcell <- as.numeric(D$comp %in% c("epithelial", "unassigned"))
D$s_surv <- ifelse(D$surv_dir & D$surv_n >= 2, 2, ifelse(D$surv_n >= 1, 1, 0))
D$s_novel <- ifelse(D$pdac_best == 0, 2, ifelse(D$pdac_best <= 5, 1, 0))
D$score <- rowSums(D[, c("s_effect", "s_consist", "s_rna", "s_tumorcell", "s_surv", "s_novel")], na.rm = TRUE)
D <- D[order(-D$score, -abs(D$lfc)), ]
write.csv(D, "results/meta/gene_prioritization.csv", row.names = FALSE)
cat("Compartment calls (rows = direction in tumor):\n"); print(table(D$dir, D$comp, useNA = "ifany"))
cols <- c("gene", "lfc", "I2", "k", "comp", "pdac_best", "surv_n", "epigenetic", "score")
cat("\nTop 40 by priority score:\n"); print(head(D[, cols], 40), row.names = FALSE, digits = 3)
cat("\nEpigenetic-flagged DEGs:", sum(D$epigenetic), "| top 15 by score:\n"); print(head(D[D$epigenetic, cols], 15), row.names = FALSE, digits = 3)
cat("\nTumor-up genes with <=2 PDAC papers, with subtype correlations:\n")
print(head(D[D$dir == "up" & D$pdac_best <= 2, c("gene", "lfc", "comp", "r_classical", "r_basal", "pdac_best")], 15), row.names = FALSE, digits = 3)
