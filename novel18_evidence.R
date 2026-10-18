suppressMessages({library(Seurat); library(jsonlite)})
genes <- read.csv("results/meta/module3_not_in_published_lists.csv")$gene

## 1. bulk meta-analysis + pooled survival
cx <- read.csv("results/survival/cox_arrays_robust.csv")
sv <- read.csv("results/survival/survival_meta_fdr.csv")
cat("survival_meta_fdr columns:", names(sv), "\n")
bulk <- merge(data.frame(gene=genes),
              cx[, intersect(c("gene","lfc_core","I2_core","HR_T","p_pro"), names(cx))],
              by="gene", all.x=TRUE)
keep <- setdiff(unique(c("gene", grep("fdr|FDR|pool|n_sig|same", names(sv), value=TRUE))),
                setdiff(names(bulk), "gene"))
bulk <- merge(bulk, sv[, keep, drop=FALSE], by="gene", all.x=TRUE)

## 2. Europe PMC: papers with gene AND pancrea* in title/abstract
epmc <- function(g){
  q <- sprintf('TITLE_ABS:"%s" AND TITLE_ABS:"pancrea*"', g)
  u <- paste0("https://www.ebi.ac.uk/europepmc/webservices/rest/search?query=",
              URLencode(q, reserved=TRUE), "&format=json&pageSize=1&resultType=lite")
  Sys.sleep(0.3); tryCatch(fromJSON(u)$hitCount, error=function(e) NA)
}
bulk$epmc_pdac_papers <- sapply(bulk$gene, epmc)
bulk$epmc_note <- ifelse(bulk$gene %in% c("PI3","F5","NHS"), "ambiguous symbol", "")

## 3. single-cell localisation
sc_stats <- function(path, label){
  o <- readRDS(path); DefaultAssay(o) <- "RNA"; o <- JoinLayers(o)
  md <- o@meta.data
  cc <- names(md)[sapply(md, function(x) any(grepl("PDAC|Norm|ADJ", as.character(x))))][1]
  cat(label, ": sample column =", cc, "| celltypes:", paste(unique(md$celltype), collapse=", "), "\n")
  tum <- grepl("PDAC", as.character(md[[cc]]))
  duct <- grepl("ductal", md$celltype, ignore.case=TRUE)
  g <- intersect(genes, rownames(o))
  X <- GetAssayData(o, layer="data")[g, , drop=FALSE]
  pct <- function(m) if (sum(m) == 0) rep(NA, length(g)) else Matrix::rowMeans(X[, m, drop=FALSE] > 0) * 100
  res <- data.frame(gene=g,
    duct_tumor = pct(duct & tum), duct_normal = pct(duct & !tum),
    nonduct_tumor = pct(!duct & tum))
  nd <- setdiff(unique(md$celltype[tum & !duct]), NA)
  M <- sapply(nd, function(ct) pct(tum & md$celltype == ct)); rownames(M) <- g
  res$top_nonduct <- apply(M, 1, function(r) paste0(nd[which.max(r)], " ", round(max(r, na.rm=TRUE))))
  res[, -1] <- lapply(res[, -1, drop=FALSE], function(x) if (is.numeric(x)) round(x, 1) else x)
  names(res)[-1] <- paste0(label, "_", names(res)[-1]); res
}
a <- sc_stats("sc/seurat_obj.rds", "sc1")
b <- sc_stats("sc2/GSE212966/seurat_obj.rds", "sc2")
out <- Reduce(function(x, y) merge(x, y, by="gene", all.x=TRUE), list(bulk, a, b))
write.csv(out, "results/meta/module3_novel18_evidence.csv", row.names=FALSE)
options(width=250); print(out, row.names=FALSE)
