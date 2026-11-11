# PRELIMINARY CHECK ONLY (69 genes, mucins included, mean scores). Not used for the results reported in the manuscript.
suppressMessages({library(Seurat)})
obj <- readRDS("sc/seurat_obj.rds")
ma <- read.csv("results/meta/module_assignments_ward_corrected.csv")
cat("Columns:", colnames(ma), "\n")
gc <- intersect(c("gene","symbol","Gene","SYMBOL"), colnames(ma))[1]
mc <- grep("module|cluster", colnames(ma), value=TRUE, ignore.case=TRUE)[1]
cat("Module column:", mc, "| values:", paste(unique(ma[[mc]]), collapse=", "), "\n")
sel <- as.character(ma[[mc]]) == "1"
g3 <- intersect(ma[[gc]][sel], rownames(obj))
cat("Corrected Module 1 genes selected:", sum(sel), "| found in scRNA:", length(g3), "\n")
obj <- AddModuleScore(obj, features=list(g3), name="gast", seed=1)
d <- subset(obj, subset = celltype == "Ductal")
key <- c("TFF1","CTSE","VSIG2","CAPN8","ARNTL2")
key <- intersect(key, rownames(d))
ex <- GetAssayData(d, layer="data")[key, , drop=FALSE]
res <- do.call(rbind, lapply(split(seq_len(ncol(d)), d$sample), function(i) {
  data.frame(ductal_cells=length(i),
    mean_score=round(mean(d$gast1[i]), 3),
    pct_score_pos=round(100*mean(d$gast1[i] > 0), 1),
    t(round(100*rowMeans(ex[, i, drop=FALSE] > 0), 1)))
}))
res$group <- ifelse(grepl("Norm", rownames(res)), "Normal", "PDAC")
print(res[order(res$group), ])
write.csv(res, "results/sc/module1_ward_corrected_per_sample.csv")
a <- res$mean_score[res$group=="PDAC"]; b <- res$mean_score[res$group=="Normal"]
cat("\nPDAC samples with mean score above the highest normal:",
    sum(a > max(b)), "of", length(a), "\n")
