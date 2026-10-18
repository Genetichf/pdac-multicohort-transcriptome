suppressMessages({library(Seurat)})
o <- readRDS("sc2/GSE212966/seurat_obj.rds"); o <- JoinLayers(o)
genes <- c("FAR2","PI3","NHS","ELOVL6","ACACB","ARNTL2","TPX2","TFF1","AGR2","LGALS4","CTSE","ANXA10","GPX2","TSPAN8","VSIG2","CAPN8","TFF2","ZG16B","ST6GALNAC1","TMEM45B","EPS8L3","VSIG1","F5","ADAM28","ZG16B","MUC20","TCN1","CEACAM7","DUOX2","OLFM4","VNN1")
genes <- unique(genes[genes %in% rownames(o)])
m <- GetAssayData(o, assay = "RNA", layer = "counts")[genes, , drop = FALSE]
grp <- ifelse(o$celltype == "Ductal", paste0("Ductal_", o$group), "NonDuctal")
out <- do.call(rbind, lapply(unique(grp), function(g) {
  i <- which(grp == g)
  data.frame(gene = genes, group = g, n_cells = length(i), pct_expressing = round(100 * rowMeans(m[, i, drop = FALSE] > 0), 1))
}))
write.csv(out, "results/sc/atlas2_target_expression.csv", row.names = FALSE)
print(subset(out, gene %in% c("FAR2","PI3","NHS","ELOVL6")))
