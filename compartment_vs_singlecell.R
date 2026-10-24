suppressMessages({ library(Seurat); library(Matrix) })
D <- read.csv("results/meta/robust_funnel_table.csv", stringsAsFactors = FALSE)[, c("gene", "comp", "lfc_core")]
D$dir <- ifelse(D$lfc_core > 0, "up", "down")
cat("Genes:", nrow(D), "| marker-based compartments:\n"); print(table(D$comp))
map1 <- c(Ductal = "epithelial", Acinar = "acinar", Fibroblast = "stroma", Endothelial = "stroma",
          B_plasma = "immune", Myeloid = "immune", T_NK = "immune", Endocrine = "endocrine/other")

topcomp <- function(o, ct) {
  g <- intersect(D$gene, rownames(o))
  m <- tryCatch(GetAssayData(o, assay = "RNA", layer = "data"), error = function(e) GetAssayData(o, assay = "RNA", slot = "data"))[g, , drop = FALSE]
  lv <- sort(unique(ct))
  em <- sapply(lv, function(l) Matrix::rowMeans(m[, which(ct == l), drop = FALSE]))
  top <- lv[apply(em, 1, which.max)]
  srt <- t(apply(em, 1, sort, decreasing = TRUE))
  data.frame(gene = g, top_celltype = top, sc_comp = unname(map1[top]),
             clear = (srt[, 1] + 0.01) / (srt[, 2] + 0.01) >= 1.5, stringsAsFactors = FALSE)
}
agree <- function(df, label) {
  a <- df[df$comp %in% c("epithelial", "stroma", "immune", "acinar") & !is.na(df$sc_comp), ]
  cat("\n=====", label, "| genes found:", sum(!is.na(df$sc_comp)), "| marker-assigned genes compared:", nrow(a), "\n")
  print(table(marker = a$comp, singlecell = a$sc_comp))
  cat("Overall agreement:", round(100 * mean(a$comp == a$sc_comp), 1), "%",
      "| among genes with a clear top cell type (>=1.5x):", round(100 * mean(a$comp[a$clear] == a$sc_comp[a$clear]), 1), "% of", sum(a$clear), "\n")
  cat("Agreement by marker compartment (%):\n"); print(round(100 * tapply(a$comp == a$sc_comp, a$comp, mean), 1))
  k <- function(lab, i) cat(lab, ": n =", sum(i), "| top cell type matches:", round(100 * mean(a$comp[i] == a$sc_comp[i]), 1), "%\n")
  k("Tumor-down genes assigned acinar", a$dir == "down" & a$comp == "acinar")
  k("Tumor-up genes assigned stroma", a$dir == "up" & a$comp == "stroma")
  k("Tumor-up genes assigned epithelial", a$dir == "up" & a$comp == "epithelial")
  u <- df[df$comp == "unassigned" & !is.na(df$sc_comp), ]
  cat("Genes the markers left unassigned (n =", nrow(u), ") by single-cell top compartment:\n"); print(table(u$sc_comp))
}
o1 <- readRDS("sc/seurat_obj.rds"); DefaultAssay(o1) <- "RNA"
r1 <- topcomp(o1, sub("\\|.*$", "", as.character(o1$ct_group))); rm(o1); gc(verbose = FALSE)
o2 <- readRDS("sc2/GSE212966/seurat_obj.rds"); DefaultAssay(o2) <- "RNA"
r2 <- topcomp(o2, as.character(o2$celltype)); rm(o2); gc(verbose = FALSE)
names(r1)[-1] <- paste0(names(r1)[-1], "_atlas1"); names(r2)[-1] <- paste0(names(r2)[-1], "_atlas2")
M <- merge(merge(D, r1, by = "gene", all.x = TRUE), r2, by = "gene", all.x = TRUE)
agree(data.frame(M[, c("gene", "comp", "dir")], sc_comp = M$sc_comp_atlas1, clear = M$clear_atlas1), "ATLAS 1 (GSE155698)")
agree(data.frame(M[, c("gene", "comp", "dir")], sc_comp = M$sc_comp_atlas2, clear = M$clear_atlas2), "ATLAS 2 (GSE212966)")
both <- M$comp %in% c("epithelial", "stroma", "immune", "acinar") & !is.na(M$sc_comp_atlas1) & M$sc_comp_atlas1 == M$sc_comp_atlas2
cat("\nGenes where both atlases give the same top compartment:", sum(both, na.rm = TRUE),
    "| of these, marker-based agrees:", sum(both & M$comp == M$sc_comp_atlas1, na.rm = TRUE),
    "| disagrees:", sum(both & M$comp != M$sc_comp_atlas1, na.rm = TRUE), "\n")
dis <- M[both & M$comp != M$sc_comp_atlas1, c("gene", "comp", "dir", "top_celltype_atlas1", "top_celltype_atlas2")]
cat("Examples of disagreement (marker vs both atlases):\n"); print(head(dis, 25), row.names = FALSE)
write.csv(M, "results/meta/compartment_vs_singlecell.csv", row.names = FALSE)
