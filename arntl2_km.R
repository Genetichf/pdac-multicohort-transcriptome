suppressMessages({library(GEOquery); library(survival)})
gene <- "ARNTL2"
get_surv <- function(id) {
  g <- getGEO(id, GSEMatrix = TRUE, getGPL = FALSE)[[1]]; pd <- pData(g)
  pick <- function(rx) { cn <- grep(rx, colnames(pd), value = TRUE, ignore.case = TRUE)[1]; if (is.na(cn)) stop("no column ", rx, " in ", id); pd[[cn]] }
  tm <- as.numeric(as.character(pick("^survival[ _]month")))
  st <- as.character(pick("^(survival status|cancer_death)"))
  cat(id, "status values:", paste(unique(st), collapse = " | "), "\n")
  ev <- as.numeric(grepl("^(1|dead|deceased|death|yes)", trimws(tolower(st))))
  data.frame(gsm = pd$geo_accession, time = tm, event = ev)
}
ids <- c("GSE28735", "GSE62452"); plots <- list()
png("FigS3_KM_ARNTL2.png", width = 2400, height = 1100, res = 300)
par(mfrow = c(1, 2), mar = c(4.2, 4.2, 3, 1))
for (id in ids) {
  sv <- get_surv(id)
  ex <- read.csv(file.path("results", id, "expression_gene_level.csv"), check.names = FALSE)
  if (!is.numeric(ex[[1]])) { rownames(ex) <- ex[[1]]; ex <- ex[, -1] }
  v <- as.numeric(ex[gene, ]); names(v) <- colnames(ex)
  sv$expr <- v[sv$gsm]; sv <- sv[complete.cases(sv), ]
  cat(id, ": n =", nrow(sv), "events =", sum(sv$event), "\n")
  sv$grp <- factor(ifelse(sv$expr > median(sv$expr), "High", "Low"), levels = c("Low", "High"))
  lr <- survdiff(Surv(time, event) ~ grp, sv); p <- 1 - pchisq(lr$chisq, 1)
  cx <- summary(coxph(Surv(time, event) ~ scale(expr), sv))$coefficients
  cat(id, ": logrank p =", signif(p, 3), "| HR per SD =", round(exp(cx[1, 1]), 2), "p =", signif(cx[1, 5], 3), "\n")
  plot(survfit(Surv(time, event) ~ grp, sv), col = c("#2c7bb6", "#d7191c"), lwd = 2, xlab = "Months", ylab = "Survival probability",
       main = paste0(gene, " - ", id, " (n=", nrow(sv), ")"))
  legend("topright", c("Low (<= median)", "High (> median)"), col = c("#2c7bb6", "#d7191c"), lwd = 2, bty = "n")
  text(0, 0.05, sprintf("log-rank p = %.3g\nHR/SD = %.2f", p, exp(cx[1, 1])), adj = 0)
}
dev.off()
