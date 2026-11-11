
suppressMessages(library(ggplot2))

dir.create("results/figures", showWarnings = FALSE, recursive = TRUE)

input <- "results/meta/module1_candidate_evidence_ward_corrected_epmc_fixed.csv"

e <- read.csv(input, stringsAsFactors = FALSE, check.names = FALSE)

required <- c(
  "gene",
  "sc1_duct_tumor",
  "sc1_duct_normal",
  "sc2_duct_tumor",
  "sc2_duct_normal"
)

missing <- setdiff(required, names(e))

if (length(missing) > 0) {
  stop("Missing required columns: ", paste(missing, collapse = ", "))
}

e <- e[!duplicated(e$gene), ]

if (nrow(e) != 22) {
  stop("Expected 22 candidate genes, but found ", nrow(e))
}

cat("Candidate genes:", nrow(e), "\n")

# Prepare expression percentages from both independent atlases
lg <- rbind(
  data.frame(
    gene = e$gene,
    atlas = "Atlas 1 (GSE155698)",
    tumor = as.numeric(e$sc1_duct_tumor),
    normal = as.numeric(e$sc1_duct_normal)
  ),
  data.frame(
    gene = e$gene,
    atlas = "Atlas 2 (GSE212966)",
    tumor = as.numeric(e$sc2_duct_tumor),
    normal = as.numeric(e$sc2_duct_normal)
  )
)

if (anyNA(lg$tumor) || anyNA(lg$normal)) {
  stop("Missing ductal expression percentages detected.")
}

# Order genes according to their mean tumor-normal difference
delta <- data.frame(
  gene = lg$gene,
  difference = lg$tumor - lg$normal
)

ord <- aggregate(difference ~ gene, delta, mean)
ord <- ord[order(ord$difference), ]

lg$gene <- factor(lg$gene, levels = ord$gene)

# Highlight exploratory candidates, not validated biomarkers
highlight <- c("NHS", "CORO2A", "MTMR11")

p <- ggplot(lg, aes(y = gene)) +

  geom_segment(
    aes(x = normal, xend = tumor, yend = gene),
    colour = "grey65",
    linewidth = 0.55
  ) +

  geom_point(
    aes(x = normal, colour = "Normal / adjacent ducts"),
    size = 2.4
  ) +

  geom_point(
    aes(x = tumor, colour = "Tumor ducts"),
    size = 2.4
  ) +

  scale_colour_manual(
    values = c(
      "Normal / adjacent ducts" = "#0072B2",
      "Tumor ducts" = "#D55E00"
    )
  ) +

  scale_y_discrete(
    labels = function(x) {
      ifelse(x %in% highlight, paste0(x, " *"), x)
    }
  ) +

  scale_x_continuous(
    limits = c(0, 100),
    breaks = seq(0, 100, 20),
    expand = expansion(mult = c(0.01, 0.03))
  ) +

  facet_wrap(~atlas, nrow = 1) +

  labs(
    x = "Ductal cells expressing the gene (%)",
    y = NULL,
    colour = NULL,
    title = "Ductal expression of 22 corrected Module 1 candidates",
    subtitle = "Comparison of tumor and non-tumor ductal cells across two independent PDAC atlases",
    caption = "* NHS, CORO2A and MTMR11 are highlighted as exploratory candidates for further investigation."
  ) +

  theme_minimal(base_size = 10) +

  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major.y = element_blank(),
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 9),
    legend.position = "bottom",
    strip.text = element_text(face = "bold"),
    axis.text.y = element_text(size = 8.5),
    plot.caption = element_text(size = 8, hjust = 0)
  )

ggsave(
  "results/figures/Fig6_candidates_ward_corrected.pdf",
  p,
  width = 11,
  height = 8.5
)

ggsave(
  "results/figures/Fig6_candidates_ward_corrected.png",
  p,
  width = 11,
  height = 8.5,
  dpi = 600,
  bg = "white"
)

cat("\nCorrected Figure 6 generated successfully.\n")
cat("Genes plotted:", length(unique(lg$gene)), "\n")
cat("Files saved:\n")
cat("Fig6_candidates_ward_corrected.pdf\n")
cat("Fig6_candidates_ward_corrected.png\n")
