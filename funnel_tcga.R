suppressPackageStartupMessages({ library(SummarizedExperiment); library(recount3); library(edgeR); library(survival) })
S <- read.csv("results/meta/meta_robust_DEGs.csv", stringsAsFactors = FALSE)
S <- S[S$k_core >= 3 & !grepl("^(LOC|LINC|MIR|SNOR)", S$gene), ]
mk <- list(stroma = c("COL1A1","COL1A2","COL3A1","DCN","LUM","FAP"), immune = c("PTPRC","CD53","LAPTM5","CORO1A","CD3E","CD68"),
           epithelial = c("KRT19","KRT7","EPCAM","KRT8","KRT18","CLDN4"), acinar = c("PRSS1","CPA1","PNLIP","CTRB1"),
           prolif = c("MKI67","CCNB1","CDK1","MCM2","PCNA","AURKA"), classical = c("GATA6","TFF1","AGR2","LGALS4"), basal = c("KRT5","KRT6A","TP63","S100A2"))
rt <- readRDS("results/rnaseq/rse_tcga_PAAD.rds")
assay(rt, "counts") <- compute_read_counts(rt)
rt <- rt[, rt$tcga.cgc_sample_sample_type %in% "Primary Tumor" & rt$tcga.cgc_case_histological_diagnosis %in% "Pancreas-Adenocarcinoma Ductal Type"]
num <- function(x) suppressWarnings(as.numeric(as.character(x)))
dead <- rt$tcga.gdc_cases.diagnoses.vital_status == "dead"
time <- ifelse(dead, num(rt$tcga.gdc_cases.diagnoses.days_to_death), num(rt$tcga.gdc_cases.diagnoses.days_to_last_follow_up)) / 30.44
stage <- tolower(as.character(rt$tcga.gdc_cases.diagnoses.tumor_stage))
adv <- ifelse(stage %in% c("stage iii", "stage iv"), 1, ifelse(is.na(stage) | stage == "not reported", NA, 0))
age <- num(rt$tcga.gdc_cases.diagnoses.age_at_diagnosis) / 365.25
ok <- is.finite(time) & time > 0 & !is.na(dead)
rt <- rt[, ok]; time <- time[ok]; event <- as.numeric(dead[ok]); adv <- adv[ok]; age <- age[ok]
cat("Tumors:", ncol(rt), "| deaths:", sum(event), "| genes in robust set:", nrow(S), "\n")
rd <- rowData(rt)
lcpm <- cpm(calcNormFactors(DGEList(assay(rt, "counts"))), log = TRUE, prior.count = 1)
pc <- rd$gene_type %in% "protein_coding" & rd$gene_name %in% unique(c(S$gene, unlist(mk)))
lcpm <- lcpm[pc, ]; sym <- rd$gene_name[pc]
mu <- rowMeans(lcpm); best <- tapply(seq_along(sym), sym, function(i) i[which.max(mu[i])])
X <- lcpm[unlist(best), ]; rownames(X) <- names(best)
sc <- sapply(mk, function(g) { g <- intersect(g, rownames(X)); colMeans(t(scale(t(X[g, , drop = FALSE])))) })
pz <- function(excl) { g <- setdiff(intersect(mk$prolif, rownames(X)), excl); colMeans(t(scale(t(X[g, , drop = FALSE])))) }
pro_all <- pz(character(0))
G <- intersect(S$gene, rownames(X)); cat("Genes found in TCGA data:", length(G), "\n")
r <- t(sapply(G, function(g) as.numeric(cor(X[g, ], sc, method = "spearman")))); colnames(r) <- paste0("r_", colnames(sc))
main <- c("stroma", "immune", "epithelial", "acinar"); rm4 <- r[, paste0("r_", main)]
comp <- ifelse(apply(rm4, 1, max) >= 0.4, main[apply(rm4, 1, which.max)], "unassigned")
cx <- t(sapply(G, function(g) {
  pro <- if (g %in% mk$prolif) pz(g) else pro_all
  d <- data.frame(time, event, z = as.numeric(scale(X[g, ])), pro = as.numeric(scale(pro)), adv, age)
  f <- function(fm) { m <- try(coxph(fm, data = d), silent = TRUE); if (inherits(m, "try-error")) return(c(NA, NA)); s <- summary(m)$coefficients; c(exp(s["z", "coef"]), s["z", "Pr(>|z|)"]) }
  a <- f(Surv(time, event) ~ z); b <- f(Surv(time, event) ~ z + adv + age); c3 <- f(Surv(time, event) ~ z + pro)
  c(HR = a[1], p = a[2], HR_adj = b[1], p_adj = b[2], HR_pro = c3[1], p_pro = c3[2])
}))
Tf <- data.frame(gene = G, comp = comp, round(r, 2), round(cx, 4), row.names = NULL)
write.csv(Tf, "results/meta/robust_tcga_features.csv", row.names = FALSE)
dirn <- ifelse(S$lfc_core[match(G, S$gene)] > 0, "up", "down")
cat("\nCompartment by direction in tumor:\n"); print(table(dirn, comp))
cat("\nGenes with TCGA p<0.05:", sum(Tf$p < 0.05, na.rm = TRUE), "| also p<0.05 after proliferation adjustment:", sum(Tf$p < 0.05 & Tf$p_pro < 0.05, na.rm = TRUE), "\n")
