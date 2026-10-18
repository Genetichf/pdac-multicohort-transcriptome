if (!requireNamespace("edgeR", quietly = TRUE)) BiocManager::install("edgeR", update = FALSE, ask = FALSE)
suppressPackageStartupMessages({ library(SummarizedExperiment); library(recount3); library(edgeR); library(limma) })
rg <- readRDS("results/rnaseq/rse_gtex_PANCREAS.rds"); rt <- readRDS("results/rnaseq/rse_tcga_PAAD.rds")
assay(rg, "counts") <- compute_read_counts(rg); assay(rt, "counts") <- compute_read_counts(rt)
st <- rt$tcga.cgc_sample_sample_type; hs <- rt$tcga.cgc_case_histological_diagnosis
is_T <- st %in% "Primary Tumor" & hs %in% "Pancreas-Adenocarcinoma Ductal Type"
is_N <- st %in% "Solid Tissue Normal"
stopifnot(identical(rowData(rg)$gene_id, rowData(rt)$gene_id))
cnt <- cbind(assay(rg, "counts"), assay(rt, "counts")[, is_T | is_N])
grp <- factor(c(rep("GTEx", ncol(rg)), ifelse(is_T[is_T | is_N], "TCGA_T", "TCGA_N")), levels = c("GTEx", "TCGA_N", "TCGA_T"))
cat("Groups:\n"); print(table(grp))
rd <- rowData(rg)
if (!all(c("gene_name", "gene_type") %in% names(rd))) stop("rowData columns: ", paste(names(rd), collapse = ", "))
pc <- rd$gene_type == "protein_coding"
cnt <- cnt[pc, ]; sym <- rd$gene_name[pc]
keep <- filterByExpr(cnt, group = grp)
cnt <- cnt[keep, ]; sym <- sym[keep]
mu <- rowMeans(cpm(cnt, log = TRUE))
best <- tapply(seq_along(sym), sym, function(i) i[which.max(mu[i])])
cnt <- cnt[unlist(best), ]; rownames(cnt) <- names(best)
cat("Protein-coding genes tested:", nrow(cnt), "\n")
y <- calcNormFactors(DGEList(cnt))
des <- model.matrix(~ 0 + grp); colnames(des) <- levels(grp)
v <- voom(y, des)
cm <- makeContrasts(T_vs_GTEx = TCGA_T - GTEx, T_vs_N = TCGA_T - TCGA_N, N_vs_GTEx = TCGA_N - GTEx, levels = des)
fit <- eBayes(contrasts.fit(lmFit(v, des), cm), robust = TRUE)
tt <- function(co) { x <- topTable(fit, coef = co, number = Inf, sort.by = "none"); data.frame(gene = rownames(x), logFC = x$logFC, fdr = x$adj.P.Val) }
A <- tt("T_vs_GTEx"); B <- tt("T_vs_N"); C <- tt("N_vs_GTEx")
R <- data.frame(gene = A$gene, lfc_T_GTEx = A$logFC, fdr_T_GTEx = A$fdr, lfc_T_N = B$logFC, fdr_T_N = B$fdr, lfc_N_GTEx = C$logFC, fdr_N_GTEx = C$fdr)
write.csv(R, "results/rnaseq/rnaseq_limma_voom_3contrasts.csv", row.names = FALSE)
isdeg <- R$fdr_T_GTEx < 0.05 & abs(R$lfc_T_GTEx) >= 1
cat("RNA-seq DEGs, tumor vs GTEx (FDR<0.05, |log2FC|>=1):", sum(isdeg), "| up:", sum(isdeg & R$lfc_T_GTEx > 0), "| down:", sum(isdeg & R$lfc_T_GTEx < 0), "\n")
Mt <- read.csv("results/meta/meta_three_datasets_all_genes.csv", stringsAsFactors = FALSE)
S <- read.csv("results/meta/meta_three_datasets_DEGs.csv", stringsAsFactors = FALSE)
j <- merge(R, Mt[, c("gene", "lfc")], by = "gene")
cat("Genes shared with the array meta-analysis:", nrow(j), "| logFC correlation Pearson/Spearman:",
    round(cor(j$lfc_T_GTEx, j$lfc), 3), "/", round(cor(j$lfc_T_GTEx, j$lfc, method = "spearman"), 3), "\n")
k <- merge(S[, c("gene", "lfc")], R, by = "gene")
same <- sign(k$lfc) == sign(k$lfc_T_GTEx)
cat("\nArray meta-DEGs present in RNA-seq:", nrow(k), "of", nrow(S), "\n")
cat("  same direction:", sum(same), "| and FDR<0.05:", sum(same & k$fdr_T_GTEx < 0.05), "| and |log2FC|>=1:", sum(same & k$fdr_T_GTEx < 0.05 & abs(k$lfc_T_GTEx) >= 1), "\n")
cat("  also supported in tumor vs TCGA normal (n=4), same direction:", sum(sign(k$lfc) == sign(k$lfc_T_N)), "\n")
cat("  flagged as possible cohort effect (|log2FC|>=1 and FDR<0.05 between the two normal sets, same sign as tumor-vs-GTEx):",
    sum(k$fdr_N_GTEx < 0.05 & abs(k$lfc_N_GTEx) >= 1 & sign(k$lfc_N_GTEx) == sign(k$lfc_T_GTEx)), "\n")
G <- read.csv("results/survival/cox_GSE28735_GSE62452.csv", stringsAsFactors = FALSE)$gene
M <- read.csv("results/survival/cox_GSE28735_GSE62452.csv", stringsAsFactors = FALSE)
cat("\nThe 15 prognostic candidates in RNA-seq (tumor vs GTEx):\n")
print(R[R$gene %in% M$gene[M$both], c("gene", "lfc_T_GTEx", "fdr_T_GTEx")], row.names = FALSE, digits = 3)
