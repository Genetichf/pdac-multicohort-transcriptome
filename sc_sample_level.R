suppressMessages({ library(Seurat); library(Matrix) })
keyg <- c("FAR2", "PI3", "NHS", "ELOVL6", "TFF1", "CTSE", "VSIG2", "CAPN8", "ARNTL2", "MUC5AC", "MUC13")
ma <- read.csv("results/meta/module_assignments.csv"); mod3 <- ma$gene[as.character(ma$module) == "3"]
mucins <- c("MUC5AC", "MUC13", "MUC4", "MUC1", "MUC6", "MUC2")
MIN <- 20

analyse <- function(o, ct, cond, smp, label) {
  DefaultAssay(o) <- "RNA"
  keep <- ct == "Ductal"
  d <- o[, keep]; d$cond <- cond[keep]; d$smp <- smp[keep]
  g3 <- setdiff(intersect(mod3, rownames(d)), mucins)
  d <- AddModuleScore(d, features = list(g3), name = "gast", seed = 1)
  genes <- intersect(keyg, rownames(d))
  X <- as.matrix(tryCatch(GetAssayData(d, assay = "RNA", layer = "data"),
                          error = function(e) GetAssayData(d, assay = "RNA", slot = "data"))[genes, , drop = FALSE])
  cat("\n=====", label, "| ductal cells:", ncol(d), "| module genes used:", length(g3), "\n")
  cat("Ductal cells per sample and condition:\n"); print(table(sample = d$smp, condition = d$cond))
  idx <- split(seq_len(ncol(d)), d$smp)
  ps <- do.call(rbind, lapply(names(idx), function(s) {
    i <- idx[[s]]
    data.frame(sample = s, cond = d$cond[i[1]], n_ductal = length(i), module = mean(d$gast1[i]),
               t(setNames(rowMeans(X[, i, drop = FALSE]), paste0("mean_", genes))),
               t(setNames(100 * rowMeans(X[, i, drop = FALSE] > 0), paste0("pct_", genes))), check.names = FALSE)
  }))
  ps$used <- ps$n_ductal >= MIN
  cat("\nSamples excluded (<", MIN, "ductal cells):", paste(ps$sample[!ps$used], collapse = ", "), "\n")
  u <- ps[ps$used, ]
  tum <- u$cond == "PDAC"
  cat("Tumor samples used:", sum(tum), "| non-tumor samples used:", sum(!tum), "\n")
  feats <- c("module", paste0("mean_", genes))
  tests <- do.call(rbind, lapply(feats, function(f) {
    a <- u[[f]][tum]; b <- u[[f]][!tum]
    data.frame(feature = sub("^mean_", "", f), median_tumor = median(a), median_nontumor = median(b),
               n_tumor_above_all_nontumor = sum(a > max(b)), n_tumor = length(a),
               p_wilcoxon = suppressWarnings(wilcox.test(a, b, exact = TRUE)$p.value))
  }))
  tests[, 2:3] <- lapply(tests[, 2:3], round, 3); tests$p_wilcoxon <- signif(tests$p_wilcoxon, 2)
  cat("\nSample-level test (each sample = one observation):\n"); print(tests, row.names = FALSE)
  write.csv(ps, sprintf("results/sc/sample_level_%s.csv", label), row.names = FALSE)
  write.csv(tests, sprintf("results/sc/sample_level_tests_%s.csv", label), row.names = FALSE)
  invisible(NULL)
}

o1 <- readRDS("sc/seurat_obj.rds")
cand <- intersect(c("sample", "Sample", "orig.ident", "sample_id", "patient"), names(o1@meta.data))
if (length(cand) == 0) { print(names(o1@meta.data)); stop("No sample column found in atlas 1; send me the column names above.") }
cat("Atlas 1 sample column:", cand[1], "\n")
analyse(o1, sub("\\|.*$", "", as.character(o1$ct_group)), sub("^.*\\|", "", as.character(o1$ct_group)) |> (\(x) ifelse(x == "PDAC", "PDAC", "Normal"))(),
        as.character(o1@meta.data[[cand[1]]]), "atlas1")
rm(o1); invisible(gc())
o2 <- readRDS("sc2/GSE212966/seurat_obj.rds")
analyse(o2, as.character(o2$celltype), ifelse(o2$group == "PDAC", "PDAC", "Adjacent"), as.character(o2$sample), "atlas2")
