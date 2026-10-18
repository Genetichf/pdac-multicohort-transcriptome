suppressMessages({library(GEOquery); library(Seurat); library(ggplot2)})
gsm <- c(PDAC_1="GSM4710689",PDAC_2="GSM4710690",PDAC_3="GSM4710691",
  PDAC_4="GSM4710692",PDAC_5="GSM4710693",PDAC_6="GSM4710694",
  Norm_1="GSM4710706",Norm_2="GSM4710707",Norm_3="GSM4710708")
objs <- list()
for (s in names(gsm)) {
  d <- file.path("sc", s); dir.create(d, showWarnings=FALSE)
  if (!length(list.files(d, recursive=TRUE))) {
    f <- getGEOSuppFiles(gsm[[s]], baseDir="sc/dl", makeDirectory=FALSE)
    untar(rownames(f)[1], exdir=d)
  }
  ff <- list.files(d, recursive=TRUE, full.names=TRUE)
  mtx <- grep("matrix\\.mtx", ff, value=TRUE)[1]
  bc  <- grep("barcodes", ff, value=TRUE)[1]
  ft  <- grep("features|genes", ff, value=TRUE)[1]
  cat(s, ":", basename(mtx), basename(bc), basename(ft), "\n")
  m <- ReadMtx(mtx, cells=bc, features=ft, feature.column=2)
  o <- CreateSeuratObject(m, project=s, min.cells=3, min.features=200)
  o$sample <- s; o$group <- ifelse(grepl("Norm", s), "Normal", "PDAC")
  o[["pct.mt"]] <- PercentageFeatureSet(o, pattern="^MT-")
  objs[[s]] <- subset(o, nFeature_RNA > 300 & pct.mt < 20)
  cat("   cells kept:", ncol(objs[[s]]), "\n")
}
obj <- merge(objs[[1]], objs[-1], add.cell.ids = names(objs))
obj <- NormalizeData(obj, verbose=FALSE)
obj <- FindVariableFeatures(obj, nfeatures=2000, verbose=FALSE)
obj <- ScaleData(obj, verbose=FALSE); obj <- RunPCA(obj, verbose=FALSE)
obj <- IntegrateLayers(obj, method=CCAIntegration, orig.reduction="pca",
                       new.reduction="integ", verbose=FALSE)
obj <- FindNeighbors(obj, reduction="integ", dims=1:30, verbose=FALSE)
obj <- FindClusters(obj, resolution=0.6, verbose=FALSE)
obj <- RunUMAP(obj, reduction="integ", dims=1:30, verbose=FALSE)
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
lab <- setNames(names(mk)[max.col(cl[,-1])], cl$cluster)
obj$celltype <- unname(lab[as.character(obj$seurat_clusters)])
obj$ct_group <- paste(obj$celltype, obj$group, sep="|")
cat("\nCells by type and group:\n"); print(table(obj$celltype, obj$group))
tg <- c("ACACB","ARNTL2","TPX2","TFF1","AGR2","LGALS4","CTSE","ANXA10","GPX2",
  "TSPAN8","VSIG2","CAPN8","TFF2","ZG16B","ST6GALNAC1","TMEM45B","EPS8L3","VSIG1")
tg <- intersect(tg, rownames(obj))
dp <- DotPlot(obj, features=tg, group.by="ct_group")$data
out <- data.frame(gene=dp$features.plot, group=dp$id,
                  avg=round(dp$avg.exp,2), pct=round(dp$pct.exp,1))
write.csv(out, "results/sc/target_expression.csv", row.names=FALSE)
w <- reshape(out[,c("gene","group","pct")], idvar="group", timevar="gene", direction="wide")
colnames(w) <- sub("pct\\.", "", colnames(w)); print(w, row.names=FALSE)
ggsave("results/sc/dotplot_targets.png", DotPlot(obj, features=tg, group.by="ct_group")+RotatedAxis(), width=11, height=6)
ggsave("results/sc/umap_celltype.png", DimPlot(obj, group.by="celltype", label=TRUE), width=7, height=6)
saveRDS(obj, "sc/seurat_obj.rds")
