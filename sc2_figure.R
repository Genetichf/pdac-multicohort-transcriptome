suppressMessages({library(Seurat); library(ggplot2)})
dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)
o <- readRDS("sc2/GSE212966/seurat_obj.rds"); DefaultAssay(o) <- "RNA"; o <- JoinLayers(o)
md <- o@meta.data
cc <- names(md)[sapply(md, function(x) any(grepl("PDAC|ADJ", as.character(x))))][1]
cat("sample column:", cc, "\n")
o$tissue <- ifelse(grepl("PDAC", as.character(md[[cc]])), "PDAC", "Adjacent")
o$grp <- paste(o$celltype, o$tissue, sep = "|")
red <- grep("umap", Reductions(o), value = TRUE, ignore.case = TRUE)[1]
if (is.na(red)) {
  rd <- grep("integrated|rpca|cca", Reductions(o), value = TRUE, ignore.case = TRUE)[1]
  if (is.na(rd)) rd <- "pca"
  cat("computing UMAP from:", rd, "\n")
  o <- RunUMAP(o, reduction = rd, dims = 1:min(20, ncol(Embeddings(o, rd))), reduction.name = "umap", verbose = FALSE)
  red <- "umap"
}
p1 <- DimPlot(o, reduction = red, group.by = "celltype", label = TRUE, repel = TRUE) + ggtitle("Atlas 2 (GSE212966)")
genes <- intersect(c("ACACB","ARNTL2","TPX2","TFF1","AGR2","LGALS4","CTSE","ANXA10","GPX2","TSPAN8",
                     "VSIG2","CAPN8","TFF2","ZG16B","ST6GALNAC1","TMEM45B","EPS8L3","VSIG1"), rownames(o))
p2 <- DotPlot(o, features = genes, group.by = "grp") + RotatedAxis()
ggsave("results/figures/sc2_umap.png", p1, width = 7, height = 6, dpi = 300, bg = "white")
ggsave("results/figures/sc2_dotplot.png", p2, width = 11, height = 6, dpi = 300, bg = "white")
cat("done\n")
