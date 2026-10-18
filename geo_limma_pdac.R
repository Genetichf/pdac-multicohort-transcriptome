#!/usr/bin/env Rscript
# geo_limma_pdac.R - GEO microarray Tumor vs Normal with limma
# Usage: Rscript geo_limma_pdac.R GSE28735 results

suppressPackageStartupMessages({
  pkgs_bioc <- c("GEOquery", "limma", "Biobase")
  for (p in pkgs_bioc) {
    if (!requireNamespace(p, quietly = TRUE)) {
      if (!requireNamespace("BiocManager", quietly = TRUE))
        install.packages("BiocManager", repos = "https://cloud.r-project.org")
      BiocManager::install(p, update = FALSE, ask = FALSE)
    }
  }
  library(GEOquery); library(limma); library(Biobase)
})

args   <- commandArgs(trailingOnly = TRUE)
gse_id <- if (length(args) >= 1) args[1] else "GSE28735"
outdir <- if (length(args) >= 2) args[2] else "results"
outdir <- file.path(outdir, gse_id)
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# CONFIG per dataset: mode = "paired", "dupcor" or "unpaired"
# patient_col = pData column with patient ID (NULL = auto-detect)
# patient_rx  = optional regex with one capture group to extract the ID
CONFIG <- list(
  GSE16515 = list(mode = "dupcor", patient_col = "title", patient_rx = "^Pancreatic Sample[[:space:]]*([0-9]*)-.*$", text_cols = NULL),
  GSE28735 = list(mode = "paired", patient_col = "title", patient_rx = "^.*patient sample ([0-9]+).*$", text_cols = NULL),
  GSE62452 = list(mode = "dupcor", patient_col = "title", patient_rx = "^E0*([0-9]+)-?.*$", text_cols = NULL, patient_remap = c("9728T"="9727", "17531T"="17530", "17798T"="17797", "19927T"="19926", "22506T"="22505", "26172T"="26171", "27881T"="27880", "29164T"="29163", "29818T"="29819", "29827T"="29826", "31263T"="31262", "31393T"="31392", "32018T"="32019", "34245T"="34244", "34435T"="34434")),
  GSE15471 = list(mode = "paired", patient_col = NULL, patient_rx = NULL, text_cols = NULL)
)
cfg <- if (gse_id %in% names(CONFIG)) CONFIG[[gse_id]] else
  list(mode = "unpaired", patient_col = NULL, patient_rx = NULL, text_cols = NULL)

rx_normal <- "non-?\\s?tumou?r|adjacent|normal|control|healthy|benign|non-?[[:space:]]?malignant"
rx_tumor  <- "tumou?r|pdac|adenocarcinoma|cancer|carcinoma"

log_msg <- function(...) cat(format(Sys.time(), "[%H:%M:%S] "), ..., "\n", sep = "")

# 1. Download
log_msg("Downloading ", gse_id, " from GEO ...")
options(timeout = max(600, getOption("timeout")))
gset_list <- getGEO(gse_id, GSEMatrix = TRUE, AnnotGPL = TRUE, destdir = outdir)
if (length(gset_list) > 1) {
  log_msg("Series has ", length(gset_list), " platforms: ",
          paste(names(gset_list), collapse = ", "), " -> using the first.")
}
gset <- gset_list[[1]]
log_msg("Platform: ", annotation(gset), " | features: ", nrow(gset),
        " | samples: ", ncol(gset))

pd <- pData(gset)
write.csv(pd, file.path(outdir, "metadata_preview.csv"))

# 2. Tissue annotation
text_cols <- cfg$text_cols
if (is.null(text_cols)) {
  text_cols <- names(pd)[vapply(pd, function(x) is.character(x) || is.factor(x), logical(1))]
  text_cols <- intersect(text_cols,
    c("title", grep("characteristics", names(pd), value = TRUE), "source_name_ch1"))
}
sample_text <- apply(pd[, text_cols, drop = FALSE], 1, function(r) tolower(paste(r, collapse = " | ")))

tissue <- ifelse(grepl(rx_normal, sample_text), "Normal",
          ifelse(grepl(rx_tumor,  sample_text), "Tumor", NA))
log_msg("Tissue classification:"); print(table(tissue, useNA = "ifany"))

if (any(is.na(tissue)) || length(unique(na.omit(tissue))) < 2) {
  cat("\nCould not classify all samples. Candidate labels per column:\n")
  for (cn in text_cols) { cat("--", cn, "\n"); print(head(sort(table(pd[[cn]]), decreasing = TRUE), 8)) }
  stop("Adjust rx_normal / rx_tumor or CONFIG$text_cols and re-run.")
}

# Patient ID (needed for paired / dupcor)
patient <- rep(NA_character_, ncol(gset))
if (cfg$mode %in% c("paired", "dupcor")) {
  if (!is.null(cfg$patient_col)) {
    patient <- as.character(pd[[cfg$patient_col]])
  } else {
    cand <- grep("patient|subject|individual|case|donor|pair", names(pd), ignore.case = TRUE, value = TRUE)
    ch <- grep("characteristics", names(pd), value = TRUE)
    key_hits <- ch[vapply(ch, function(cn) any(grepl("^(patient|subject|case|pair|id)[^:]*:",
                                                      tolower(as.character(pd[[cn]])))), logical(1))]
    cand <- unique(c(cand, key_hits))
    if (length(cand) > 0) {
      log_msg("Using patient column: ", cand[1])
      patient <- as.character(pd[[cand[1]]])
    }
  }
  if (!is.null(cfg$patient_rx) && !all(is.na(patient))) {
    patient <- sub(cfg$patient_rx, "\\1", patient)
  }
  patient <- sub("^[^:]*:\\s*", "", patient)
  if (!is.null(cfg$patient_remap)) {
    hit <- patient %in% names(cfg$patient_remap)
    patient[hit] <- unname(cfg$patient_remap[patient[hit]])
  }

  if (all(is.na(patient))) {
    stop("No patient ID found. Open metadata_preview.csv, find where the pairing ",
         "is encoded, then set CONFIG$", gse_id, "$patient_col and patient_rx. ",
         "Or set mode = 'unpaired'.")
  }
}

meta <- data.frame(gsm = colnames(gset), tissue = factor(tissue, levels = c("Normal", "Tumor")),
                   patient = patient, stringsAsFactors = FALSE)
write.csv(meta, file.path(outdir, "sample_annotation_used.csv"), row.names = FALSE)

keep <- rep(TRUE, nrow(meta))
if (cfg$mode == "paired") {
  tab  <- table(meta$patient, meta$tissue)
  good <- rownames(tab)[tab[, "Normal"] >= 1 & tab[, "Tumor"] >= 1]
  keep <- meta$patient %in% good
  log_msg("Paired mode: ", length(good), " complete pairs, ", sum(!keep), " samples dropped.")
  if (length(good) < 5) stop("Fewer than 5 complete pairs - check patient parsing.")
}
gset <- gset[, keep]
meta <- meta[keep, ]

# 3. Expression scale and QC
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0, 0.25, 0.5, 0.75, 0.99, 1), na.rm = TRUE))
needs_log <- (qx[5] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0)
if (needs_log) {
  log_msg("Values look un-logged -> applying log2.")
  ex[ex <= 0] <- NA
  ex <- log2(ex)
} else {
  log_msg("Values already look log2-scaled (99th pct = ", round(qx[5], 2), ").")
}
ex <- normalizeBetweenArrays(ex, method = "quantile")

grp_col <- c(Normal = "#3B7DD8", Tumor = "#D8503B")[as.character(meta$tissue)]

pdf(file.path(outdir, "QC_plots.pdf"), width = 9, height = 7)
boxplot(ex, col = grp_col, las = 2, outline = FALSE, names = NULL,
        main = paste(gse_id, "- expression per sample"), ylab = "log2 expression")
legend("topright", legend = c("Normal", "Tumor"), fill = c("#3B7DD8", "#D8503B"), bty = "n")
plotDensities(ex, group = meta$tissue, legend = "topright", main = "Density")

ex_pca <- ex[complete.cases(ex), ]
ex_pca <- ex_pca[apply(ex_pca, 1, var) > 0, ]
pca <- prcomp(t(ex_pca), scale. = TRUE)
ve  <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
plot(pca$x[, 1], pca$x[, 2], col = grp_col, pch = 19, cex = 1.3,
     xlab = paste0("PC1 (", ve[1], "%)"), ylab = paste0("PC2 (", ve[2], "%)"),
     main = "PCA (all probes)")
legend("topright", legend = c("Normal", "Tumor"), col = c("#3B7DD8", "#D8503B"), pch = 19, bty = "n")
dev.off()

cm <- cor(ex_pca); diag(cm) <- NA
mc <- rowMeans(cm, na.rm = TRUE)
z  <- (mc - mean(mc)) / sd(mc)
out_df <- data.frame(gsm = meta$gsm, tissue = meta$tissue, mean_cor = round(mc, 4), z = round(z, 2))
write.csv(out_df[order(out_df$z), ], file.path(outdir, "QC_sample_correlation.csv"), row.names = FALSE)
if (any(z < -3)) log_msg("WARNING: possible outliers (z < -3): ",
                         paste(meta$gsm[z < -3], collapse = ", "), " - not removed automatically.")

# 4. Probe -> gene symbol, collapse to one row per gene
fd <- fData(gset)
sym <- NULL
for (cn in c("Gene symbol", "Gene Symbol", "GENE_SYMBOL", "Symbol", "gene_symbol")) {
  if (cn %in% names(fd)) { sym <- as.character(fd[[cn]]); break }
}
if (is.null(sym) && "gene_assignment" %in% names(fd)) {
  sym <- vapply(strsplit(as.character(fd$gene_assignment), " // ", fixed = TRUE),
                function(x) if (length(x) >= 2) trimws(x[2]) else NA_character_, character(1))
}
if (is.null(sym)) stop("No gene symbol column found. Columns: ", paste(names(fd), collapse = ", "))

sym <- sub("[[:space:]]*///.*$", "", sym)
sym[sym %in% c("", "---", "NA")] <- NA
log_msg("Probes with a symbol: ", sum(!is.na(sym)), " / ", length(sym))

ok  <- !is.na(sym)
ex  <- ex[ok, ]; sym <- sym[ok]; probe_ids <- rownames(ex)

mu    <- rowMeans(ex, na.rm = TRUE)
best  <- tapply(seq_along(sym), sym, function(i) i[which.max(mu[i])])
ex_g  <- ex[unlist(best), , drop = FALSE]
rownames(ex_g) <- names(best)
probe_map <- data.frame(gene = names(best), probe = probe_ids[unlist(best)])
log_msg("Genes after collapse: ", nrow(ex_g))
write.csv(ex_g, file.path(outdir, "expression_gene_level.csv"))

# 5. limma
# average technical replicates (same patient + tissue) -> one column per patient x tissue
key <- paste(meta$patient, meta$tissue, sep = "__")
if (!all(is.na(meta$patient)) && anyDuplicated(key)) {
  log_msg("Averaging ", sum(duplicated(key)), " replicate arrays (same patient + tissue).")
  ex_g <- avearrays(ex_g, ID = key)
  keep1 <- !duplicated(key)
  meta <- meta[keep1, ]
  ex_g <- ex_g[, key[keep1], drop = FALSE]
  grp_col <- c(Normal = "#3B7DD8", Tumor = "#D8503B")[as.character(meta$tissue)]
}
tissue_f <- meta$tissue
if (cfg$mode == "paired") {
  patient_f <- factor(meta$patient)
  design <- model.matrix(~ patient_f + tissue_f)
  fit <- lmFit(ex_g, design)
} else if (cfg$mode == "dupcor") {
  design <- model.matrix(~ tissue_f)
  blk <- ifelse(is.na(meta$patient), meta$gsm, meta$patient)
  dc  <- duplicateCorrelation(ex_g, design, block = blk)
  log_msg("duplicateCorrelation consensus: ", round(dc$consensus.correlation, 3))
  fit <- lmFit(ex_g, design, block = blk, correlation = dc$consensus.correlation)
} else {
  design <- model.matrix(~ tissue_f)
  fit <- lmFit(ex_g, design)
}
coef_name <- "tissue_fTumor"
fit <- eBayes(fit, trend = TRUE, robust = TRUE)
tt  <- topTable(fit, coef = coef_name, number = Inf, sort.by = "P", adjust.method = "BH")
tt  <- cbind(gene = rownames(tt), tt)
tt$probe <- probe_map$probe[match(tt$gene, probe_map$gene)]

FDR_CUT <- 0.05; LFC_CUT <- 1
tt$direction <- ifelse(tt$adj.P.Val < FDR_CUT & tt$logFC >=  LFC_CUT, "Up_in_Tumor",
                ifelse(tt$adj.P.Val < FDR_CUT & tt$logFC <= -LFC_CUT, "Down_in_Tumor", "NS"))
write.csv(tt, file.path(outdir, paste0(gse_id, "_limma_Tumor_vs_Normal.csv")), row.names = FALSE)
log_msg("DEGs (FDR<", FDR_CUT, ", |log2FC|>=", LFC_CUT, "):")
print(table(tt$direction))

# 6. Plots
pdf(file.path(outdir, "DE_plots.pdf"), width = 8, height = 7)
hist(tt$P.Value, breaks = 50, col = "grey70", border = "white",
     main = "P-value distribution", xlab = "raw p-value")

cols <- c(Up_in_Tumor = "#D8503B", Down_in_Tumor = "#3B7DD8", NS = "grey75")
plot(tt$logFC, -log10(tt$P.Value), col = cols[tt$direction], pch = 19, cex = 0.5,
     xlab = "log2 fold change (Tumor vs Normal)", ylab = "-log10(p)",
     main = paste(gse_id, "- volcano"))
abline(v = c(-LFC_CUT, LFC_CUT), lty = 2)
top_lab <- head(tt[tt$direction != "NS", ], 15)
if (nrow(top_lab) > 0)
  text(top_lab$logFC, -log10(top_lab$P.Value), top_lab$gene, pos = 3, cex = 0.6)

top_genes <- head(tt$gene[tt$direction != "NS"], 50)
if (length(top_genes) >= 2) {
  mat <- t(scale(t(ex_g[top_genes, ])))
  ann <- data.frame(Tissue = meta$tissue, row.names = colnames(mat))
  if (requireNamespace("pheatmap", quietly = TRUE)) {
    print(pheatmap::pheatmap(mat, annotation_col = ann, show_colnames = FALSE,
                             fontsize_row = 6, main = "Top 50 DEGs (row z-score)"))
  } else {
    heatmap(mat, ColSideColors = grp_col, labCol = NA, scale = "none",
            main = "Top 50 DEGs (row z-score)", cexRow = 0.5)
  }
}
dev.off()

# 7. Reproducibility
writeLines(capture.output(sessionInfo()), file.path(outdir, "sessionInfo.txt"))
log_msg("Done. Outputs in: ", normalizePath(outdir))
