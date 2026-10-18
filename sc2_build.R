suppressMessages(library(Seurat))
d0 <- "sc2/GSE212966"; ex <- file.path(d0, "files"); dir.create(ex, showWarnings=FALSE)
if (!length(list.files(ex))) untar(file.path(d0, "GSE212966_RAW.tar"), exdir=ex)
objs <- list()
for (f in list.files(ex, pattern="matrix\\.mtx\\.gz$", full.names=TRUE)) {
  b <- sub("_matrix\\.mtx\\.gz$", "", f)
  s <- sub("^.*_(PDAC[0-9]+|ADJ[0-9]+)$", "\\1", b)
  m <- ReadMtx(f, cells=paste0(b, "_barcodes.tsv.gz"),
               features=paste0(b, "_genes.tsv.gz"), feature.column=2)
  o <- CreateSeuratObject(m, project=s, min.cells=3, min.features=200)
  o$sample <- s; o$group <- ifelse(grepl("^ADJ", s), "Adjacent", "PDAC")
  o[["pct.mt"]] <- PercentageFeatureSet(o, pattern="^MT-")
  objs[[s]] <- subset(o, subset = nFeature_RNA > 300 & pct.mt < 20)
  cat(s, "cells kept:", ncol(objs[[s]]), "\n")
}
set.seed(1)
objs <- lapply(objs, function(o) if (ncol(o) > 2500) o[, sample(colnames(o), 2500)] else o)
obj <- merge(objs[[1]], objs[-1], add.cell.ids=names(objs))
rm(objs); invisible(gc())
obj <- NormalizeData(obj, verbose=FALSE)
obj <- FindVariableFeatures(obj, nfeatures=2000, verbose=FALSE)
obj <- ScaleData(obj, verbose=FALSE); obj <- RunPCA(obj, verbose=FALSE)
obj <- IntegrateLayers(obj, method=RPCAIntegration, orig.reduction="pca",
                       new.reduction="integ", verbose=FALSE)
obj <- FindNeighbors(obj, reduction="integ", dims=1:30, verbose=FALSE)
obj <- FindClusters(obj, resolution=0.6, verbose=FALSE)
obj <- JoinLayers(obj)
mk <- list(Ductal=c("KRT19","KRT7","EPCAM","CFTR","SOX9","MMP7"),
  Acinar=c("PRSS1","CPA1","CELA3A","PNLIP","CTRB1","REG1A"),
  Fibroblast=c("COL1A1","DCN","LUM","PDGFRB","COL3A1"),
  Endothelial=c("PECAM1","VWF","CLDN5","PLVAP"),
  T_NK=c("CD3E","CD2","NKG7","CD3D"), Myeloid=c("CD68","LYZ","CD14","AIF1"),
  B_plasma=c("MS4A1","CD79A","MZB1","JCHAIN"), Endocrine=c("INS","GCG","CHGA"))
mk <- lapply(mk, function(g) intersect(g, rownames(obj)))
obj <- AddModuleScore(obj, features=mk, name="ms", seed=1)
sc <- obj@meta.data[, grep("^ms[0-9]", colnames(obj@meta.data))]
colnames(sc) <- names(mk)
cl <- aggregate(sc, list(cluster=obj$seurat_clusters), mean)
lab <- setNames(names(mk)[max.col(cl[, -1])], cl$cluster)
obj$celltype <- unname(lab[as.character(obj$seurat_clusters)])
cat("\nCells by type and group:\n"); print(table(obj$celltype, obj$group))
cat("\nDuctal cells per sample:\n"); print(table(obj$sample[obj$celltype == "Ductal"]))
saveRDS(obj, file.path(d0, "seurat_obj.rds"))
