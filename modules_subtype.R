suppressPackageStartupMessages({ library(SummarizedExperiment); library(recount3); library(edgeR); library(survival) })
D <- read.csv("results/meta/robust_funnel_table.csv", stringsAsFactors = FALSE)
cand <- D[D$dir == "up" & D$comp %in% c("epithelial", "unassigned") & D$r_prolif < 0.5 & D$rna_ok, ]
prog <- D[D$dir == "up" & !is.na(D$p) & D$p < 0.05 & D$p_pro < 0.05, ]
mk <- list(classical = c("GATA6","TFF1","AGR2","LGALS4"), basal = c("KRT5","KRT6A","TP63","S100A2"), prolif = c("MKI67","CCNB1","CDK1","MCM2","PCNA","AURKA"))
rt <- readRDS("results/rnaseq/rse_tcga_PAAD.rds"); assay(rt, "counts") <- compute_read_counts(rt)
rt <- rt[, rt$tcga.cgc_sample_sample_type %in% "Primary Tumor" & rt$tcga.cgc_case_histological_diagnosis %in% "Pancreas-Adenocarcinoma Ductal Type"]
num <- function(x) suppressWarnings(as.numeric(as.character(x)))
dead <- rt$tcga.gdc_cases.diagnoses.vital_status == "dead"
time <- ifelse(dead, num(rt$tcga.gdc_cases.diagnoses.days_to_death), num(rt$tcga.gdc_cases.diagnoses.days_to_last_follow_up)) / 30.44
stage <- tolower(as.character(rt$tcga.gdc_cases.diagnoses.tumor_stage))
adv <- ifelse(stage %in% c("stage iii", "stage iv"), 1, ifelse(is.na(stage) | stage == "not reported", NA, 0))
age <- num(rt$tcga.gdc_cases.diagnoses.age_at_diagnosis) / 365.25
ok <- is.finite(time) & time > 0 & !is.na(dead)
rt <- rt[, ok]; time <- time[ok]; event <- as.numeric(dead[ok]); adv <- adv[ok]; age <- age[ok]
rd <- rowData(rt); lcpm <- cpm(calcNormFactors(DGEList(assay(rt, "counts"))), log = TRUE, prior.count = 1)
pc <- rd$gene_type %in% "protein_coding" & rd$gene_name %in% unique(c(cand$gene, prog$gene, unlist(mk)))
lcpm <- lcpm[pc, ]; sym <- rd$gene_name[pc]; mu <- rowMeans(lcpm)
best <- tapply(seq_along(sym), sym, function(i) i[which.max(mu[i])]); X <- lcpm[unlist(best), ]; rownames(X) <- names(best)
sc <- sapply(mk, function(g) colMeans(t(scale(t(X[intersect(g, rownames(X)), , drop = FALSE])))))
cb <- as.numeric(scale(scale(sc[, "classical"]) - scale(sc[, "basal"]))); pro <- as.numeric(scale(sc[, "prolif"]))
cat("Tumors:", ncol(X), "| deaths:", sum(event), "| classical-basal axis built from", length(mk$classical), "+", length(mk$basal), "markers\n")
f <- function(fm, d) { m <- try(coxph(fm, data = d), silent = TRUE); if (inherits(m, "try-error")) return(c(NA, NA)); s <- summary(m)$coefficients; c(exp(s["z", "coef"]), s["z", "Pr(>|z|)"]) }
cox4 <- function(z) { d <- data.frame(time, event, z = as.numeric(scale(z)), cb, pro, adv, age)
  a <- f(Surv(time, event) ~ z, d); b <- f(Surv(time, event) ~ z + cb, d); c3 <- f(Surv(time, event) ~ z + cb + pro, d); e <- f(Surv(time, event) ~ z + cb + pro + adv + age, d)
  c(HR = a[1], p = a[2], HR_cb = b[1], p_cb = b[2], HR_cbpro = c3[1], p_cbpro = c3[2], HR_full = e[1], p_full = e[2]) }
Gc <- intersect(cand$gene, rownames(X)); Z <- t(scale(t(X[Gc, ])))
mod <- cutree(hclust(as.dist(1 - cor(t(Z), method = "spearman")), method = "ward.D2"), k = 6)
ms <- sapply(1:6, function(m) colMeans(Z[mod == m, , drop = FALSE]))
Dm <- D[match(Gc, D$gene), ]
tab <- do.call(rbind, lapply(1:6, function(m) { i <- mod == m
  data.frame(module = m, n = sum(i), n_new = sum(Dm$new[i]), r_classical = round(mean(Dm$r_classical[i]), 2), r_basal = round(mean(Dm$r_basal[i]), 2),
             r_prolif = round(mean(Dm$r_prolif[i]), 2), r_epith = round(mean(Dm$r_epithelial[i]), 2), cor_cb = round(cor(ms[, m], cb, method = "spearman"), 2), round(t(cox4(ms[, m])), 3)) }))
cat("\nModules among", length(Gc), "candidates (module score = mean z; HR per SD):\n"); print(tab, row.names = FALSE)
for (m in 1:6) { g <- Gc[mod == m]; rr <- sapply(g, function(x) cor(Z[x, ], ms[, m], method = "spearman")); top <- names(sort(rr, decreasing = TRUE))[1:min(14, length(g))]
  cat("\nModule", m, "top members (* = new vs earlier 402):", paste0(top, ifelse(Dm$new[match(top, Dm$gene)], "*", ""), collapse = ", "), "\n") }
write.csv(data.frame(gene = Gc, module = mod), "results/meta/module_assignments.csv", row.names = FALSE)
Gp <- intersect(prog$gene, rownames(X))
P <- data.frame(gene = Gp, new = D$new[match(Gp, D$gene)], r_classical = D$r_classical[match(Gp, D$gene)], r_basal = D$r_basal[match(Gp, D$gene)], round(t(sapply(Gp, function(g) cox4(X[g, ]))), 4), row.names = NULL)
P <- P[order(P$p_cbpro), ]; write.csv(P, "results/meta/prognostic_subtype_adjusted.csv", row.names = FALSE)
cat("\nPrognostic genes: HR alone vs adjusted for classical-basal axis (cb), + proliferation, + stage/age:\n"); print(P, row.names = FALSE, digits = 3)
cat("\nStill p<0.05 after cb:", sum(P$p_cb < 0.05, na.rm = TRUE), "| after cb+proliferation:", sum(P$p_cbpro < 0.05, na.rm = TRUE), "| fully adjusted:", sum(P$p_full < 0.05, na.rm = TRUE), "of", nrow(P), "\n")
