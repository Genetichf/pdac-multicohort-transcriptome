suppressMessages(library(ggplot2))
dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)
tum <- "#D55E00"; nor <- "#0072B2"           # colorblind-safe pair
th <- theme_minimal(base_size = 10) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.y = element_blank(),
        plot.title = element_text(face = "bold", size = 10), legend.position = "bottom",
        legend.title = element_blank())
save_fig <- function(plots, name, w, h, ncol = 1) {
  if (length(plots) == 1) { p <- plots[[1]]
  } else if (requireNamespace("patchwork", quietly = TRUE)) { p <- patchwork::wrap_plots(plots, ncol = ncol)
  } else { for (i in seq_along(plots)) {                       # no patchwork: save panels separately
      for (ext in c("pdf","png")) ggsave(sprintf("results/figures/%s_panel%d.%s", name, i, ext),
                                         plots[[i]], width = w, height = h / length(plots), dpi = 300)
    }; return(invisible()) }
  for (ext in c("pdf","png")) ggsave(sprintf("results/figures/%s.%s", name, ext), p, width = w, height = h, dpi = 300)
}
f <- read.csv("results/meta/robust_funnel_table.csv")

## Figure 1: signature
top <- data.frame(
  gene = c("CEACAM5","SLC6A14","POSTN","COL11A1","S100P","LAMC2","SERPINB5","INHBA","CTSE",
           "CUZD1","ALB","SERPINI2","PNLIPRP1","ERP27","PDIA2","PNLIPRP2","TMED6","AQP8"),
  lfc  = c(3.93,3.48,3.46,3.33,3.23,3.18,3.09,3.09,3.07,
           -3.72,-3.52,-3.46,-3.34,-3.11,-3.09,-3.07,-3.04,-3.00))
top$gene <- factor(top$gene, levels = top$gene[order(top$lfc)])
top$dir <- ifelse(top$lfc > 0, "Up in tumor", "Down in tumor")
p1a <- ggplot(top, aes(lfc, gene, colour = dir)) +
  geom_segment(aes(x = 0, xend = lfc, yend = gene), linewidth = 0.5) + geom_point(size = 2.2) +
  scale_colour_manual(values = c("Up in tumor" = tum, "Down in tumor" = nor)) +
  labs(x = "Pooled log2 fold change", y = NULL, title = "A. Strongest changes in the five-cohort core") + th
p1b <- ggplot(f, aes(lfc_core, I2_core, colour = dir)) + geom_point(size = 1, alpha = 0.6) +
  scale_colour_manual(values = c(up = tum, down = nor), labels = c(up = "Up", down = "Down")) +
  labs(x = "Pooled log2 fold change (core)", y = "Heterogeneity I2 (%)",
       title = "B. Effect size and heterogeneity") + th
save_fig(list(p1a, p1b), "Fig1_signature", 9, 5, ncol = 2)

## Figure 2: composition
d <- as.data.frame(table(comp = f$comp, dir = f$dir))
d$comp <- factor(d$comp, levels = c("epithelial","stroma","acinar","immune","unassigned"))
cat("Composition counts (compare with Table 3):\n"); print(table(f$comp, f$dir))
p2 <- ggplot(d, aes(comp, Freq, fill = dir)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.65) +
  geom_text(aes(label = Freq), position = position_dodge(width = 0.75), vjust = -0.3, size = 3) +
  scale_fill_manual(values = c(up = tum, down = nor), labels = c(up = "Up in tumor", down = "Down in tumor")) +
  labs(x = NULL, y = "Genes", title = "Compartment assignment of tumor-up and tumor-down genes") + th
save_fig(list(p2), "Fig2_composition", 6, 4)

## Figure 3: forest plot of the 16 pooled prognostic genes (CIs derived from p-values)
g16 <- c("TPX2","ITGA3","ARNTL2","LYZ","NT5E","MET","LAMC2","LAMA3","ACACB","ITGB6","ITGA2","F5","GPR87","RARRES3","LPCAT2","GNA15")
s <- read.csv("results/survival/survival_meta_fdr.csv"); s <- s[match(g16, s$gene), ]
mk <- function(hr, p, lab, floor) { p <- pmax(p, floor); se <- abs(log(hr)) / qnorm(1 - p / 2)
  data.frame(gene = s$gene, cohort = lab, HR = hr, lo = exp(log(hr) - 1.96 * se), hi = exp(log(hr) + 1.96 * se)) }
fo <- rbind(mk(s$HR_A, s$p_A, "GSE28735", 5e-5), mk(s$HR_B, s$p_B, "GSE62452", 5e-5),
            mk(s$HR_T, s$p_T, "TCGA", 5e-5), mk(s$HR_pooled, s$p_pooled, "Pooled", 1e-300))
fo$gene <- factor(fo$gene, levels = rev(s$gene[order(s$HR_pooled)]))
fo$cohort <- factor(fo$cohort, levels = c("GSE28735","GSE62452","TCGA","Pooled"))
p3 <- ggplot(fo, aes(HR, gene, colour = cohort, shape = cohort)) +
  geom_vline(xintercept = 1, linewidth = 0.3, colour = "grey50") +
  geom_errorbarh(aes(xmin = lo, xmax = hi), height = 0, linewidth = 0.4, position = position_dodge(width = 0.7)) +
  geom_point(size = 1.8, position = position_dodge(width = 0.7)) +
  scale_x_log10(breaks = c(0.5, 1, 2, 3)) +
  scale_colour_manual(values = c("#56B4E9","#009E73","#E69F00","black")) +
  scale_shape_manual(values = c(16, 16, 16, 18)) +
  labs(x = "Hazard ratio per SD of expression (log scale; intervals derived from p-values)", y = NULL,
       title = "Pooled survival analysis of the 16 genes surviving proliferation adjustment") + th
save_fig(list(p3), "Fig3_survival_forest", 7, 6.5)

## Figure 4: module overlap with published lists, and module correlation with classical/basal scores
o <- read.csv("results/meta/module3_overlap_stats.csv"); o <- o[o$overlap > 0, ]
o <- o[order(o$overlap), ]; o$signature <- factor(o$signature, levels = o$signature)
o$lab <- sprintf("%d / %d  (p = %s)", o$overlap, o$list_size, formatC(o$p_hyper, format = "e", digits = 1))
p4a <- ggplot(o, aes(overlap, signature)) + geom_col(fill = "#009E73", width = 0.6) +
  geom_text(aes(label = lab), hjust = -0.05, size = 2.7) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.6))) +
  labs(x = "Module genes in list (of 56)", y = NULL, title = "A. Overlap with published subtype lists") + th
mods <- data.frame(module = rep(paste0("M", 1:6), 2),
  r = c(0.45,0.01,0.55,0.21,0.12,0.13,  0.05,0.34,-0.15,0.31,0.10,0.21),
  score = rep(c("Classical","Basal"), each = 6))
p4b <- ggplot(mods, aes(module, r, fill = score)) + geom_col(position = position_dodge(width = 0.75), width = 0.65) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  scale_fill_manual(values = c(Classical = "#009E73", Basal = "#CC79A7")) +
  labs(x = "Co-expression module (M3 = gastric-like)", y = "Mean correlation with score",
       title = "B. Module correlation with subtype scores in TCGA") + th
save_fig(list(p4a, p4b), "Fig4_module_overlap", 9, 8, ncol = 1)

## Figure 6: ductal expression of the 18 module genes outside published lists
e <- read.csv("results/meta/module3_novel18_evidence.csv")
lg <- rbind(data.frame(gene = e$gene, atlas = "Atlas 1 (GSE155698)", tumor = e$sc1_duct_tumor, normal = e$sc1_duct_normal),
            data.frame(gene = e$gene, atlas = "Atlas 2 (GSE212966)", tumor = e$sc2_duct_tumor, normal = e$sc2_duct_normal))
ord <- aggregate(I(tumor - normal) ~ gene, lg, mean); ord <- ord[order(ord[[2]]), ]
lg$gene <- factor(lg$gene, levels = ord$gene)
hl <- c("FAR2","PI3","NHS","ELOVL6")
p6 <- ggplot(lg, aes(y = gene)) +
  geom_segment(aes(x = normal, xend = tumor, yend = gene), colour = "grey60", linewidth = 0.4) +
  geom_point(aes(x = normal, colour = "Normal / adjacent ducts"), size = 2) +
  geom_point(aes(x = tumor, colour = "Tumor ducts"), size = 2) +
  scale_colour_manual(values = c("Normal / adjacent ducts" = nor, "Tumor ducts" = tum)) +
  facet_wrap(~atlas) +
  labs(x = "Ductal cells expressing the gene (%)", y = NULL,
       title = "Tumor-selective expression among module genes absent from published lists") + th +
  theme(axis.text.y = element_text(face = ifelse(levels(lg$gene) %in% hl, "bold", "plain")))
save_fig(list(p6), "Fig6_candidates", 9, 5.5)

cat("\nFiles written:\n"); print(list.files("results/figures"))
