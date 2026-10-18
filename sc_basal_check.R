suppressMessages(library(Seurat))
obj <- readRDS("sc/seurat_obj.rds")
d <- subset(obj, subset = celltype == "Ductal")
basal <- intersect(c("KRT5","KRT6A","TP63","S100A2","KRT17","LAMB3"), rownames(d))
classical <- intersect(c("GATA6","TFF1","AGR2","LGALS4","CEACAM6","TFF3"), rownames(d))
d <- AddModuleScore(d, features=list(basal, classical), name="sub", seed=1)
d$balance <- d$sub2 - d$sub1
ex <- GetAssayData(d, layer="data")
show <- intersect(c("KRT5","TP63","S100A2","KRT17","GATA6","TFF3","CEACAM6"), rownames(d))
res <- do.call(rbind, lapply(split(seq_len(ncol(d)), d$sample), function(i) {
  data.frame(n=length(i),
    basal_score=round(mean(d$sub1[i]), 3),
    classical_score=round(mean(d$sub2[i]), 3),
    balance=round(mean(d$balance[i]), 3),
    t(round(100*rowMeans(ex[show, i, drop=FALSE] > 0), 1)))
}))
res$group <- ifelse(grepl("Norm", rownames(res)), "Normal", "PDAC")
print(res[order(res$group), ])
write.csv(res, "results/sc/basal_classical_per_sample.csv")
