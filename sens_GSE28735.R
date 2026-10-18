suppressPackageStartupMessages(library(limma))
d <- "results/GSE28735/"
ex   <- read.csv(paste0(d, "expression_gene_level.csv"), row.names = 1, check.names = FALSE)
meta <- read.csv(paste0(d, "sample_annotation_used.csv"), stringsAsFactors = FALSE)
full <- read.csv(paste0(d, "GSE28735_limma_Tumor_vs_Normal.csv"))
bad  <- c("GSM711905", "GSM711957")
m <- meta[!(meta$patient %in% meta$patient[meta$gsm %in% bad]), ]
e <- as.matrix(ex[, m$gsm])
fit <- eBayes(lmFit(e, model.matrix(~ factor(m$patient) + factor(m$tissue, levels = c("Normal", "Tumor")))),
              trend = TRUE, robust = TRUE)
tt <- topTable(fit, coef = ncol(fit$design), number = Inf)
tt$gene <- rownames(tt)
j <- merge(full[, c("gene", "logFC", "adj.P.Val")], tt[, c("gene", "logFC", "adj.P.Val")],
           by = "gene", suffixes = c("_full", "_sens"))
sig <- function(l, p) abs(l) >= 1 & p < 0.05
a <- sig(j$logFC_full, j$adj.P.Val_full); b <- sig(j$logFC_sens, j$adj.P.Val_sens)
cat("Pairs used:", length(unique(m$patient)), "\n")
cat("logFC correlation full vs sensitivity:", round(cor(j$logFC_full, j$logFC_sens), 4), "\n")
cat("DEGs (FDR<0.05 & |log2FC|>=1) full:", sum(a), "| sensitivity:", sum(b), "| shared:", sum(a & b), "\n")
