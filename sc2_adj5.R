suppressMessages(library(Seurat))
obj <- readRDS("sc2/GSE212966/seurat_obj.rds")
cat("Cell-type composition per sample (% of cells):\n")
print(round(100*prop.table(table(obj$sample, obj$celltype), 1), 1))
d <- subset(obj, subset = celltype == "Ductal")
g <- ifelse(d$sample == "ADJ5", "ADJ5", ifelse(d$sample == "ADJ4", "ADJ4",
     ifelse(d$group == "Adjacent", "ADJ_others", "PDAC")))
sets <- list(
  tumor_assoc = c("S100P","CEACAM6","CEACAM5","LCN2","MMP7","KRT17","CLDN18","TSPAN1","MUC4"),
  normal_duct = c("CFTR","SLC4A4","HNF1B","DCDC2","CA2","SPP1"),
  prolif      = c("MKI67","TOP2A"),
  gastric_meta = c("GKN1","GKN2","PGC","MUC6","TFF2","LYZ","TFF3","REG4"),
  intestinal  = c("CDX2","MUC2","VIL1"))
genes <- intersect(unlist(sets), rownames(d))
X <- GetAssayData(d, layer="data")[genes, , drop=FALSE]
pct <- sapply(split(seq_len(ncol(d)), g), function(i) 100*rowMeans(X[, i, drop=FALSE] > 0))
set <- rep(names(sets), lengths(sets))[match(genes, unlist(sets))]
cat("\n% of ductal cells expressing, by group (n cells):\n")
cat("  n:", paste(names(table(g)), table(g), collapse="  "), "\n")
print(data.frame(set=set, round(pct[, c("PDAC","ADJ5","ADJ4","ADJ_others")], 1)))
cat("\nDuctal cells per cluster (clusters with >=10 ductal cells):\n")
tb <- table(cluster=d$seurat_clusters, g)
tb <- tb[rowSums(tb) >= 10, , drop=FALSE]; print(tb)
sh <- tb[, "PDAC"] / rowSums(tb)
cat("\n% of ADJ5 ductal cells in clusters where >=70% of ductal cells are tumor:",
    round(100*sum(tb[sh >= 0.7, "ADJ5"]) / sum(tb[, "ADJ5"]), 1), "\n")
