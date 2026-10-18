pd <- read.csv("results/GSE62452/metadata_preview.csv", check.names = FALSE, stringsAsFactors = FALSE)
m  <- read.csv("results/GSE62452/sample_annotation_used.csv", stringsAsFactors = FALSE)
x  <- merge(m, pd[, c("geo_accession", "title", "characteristics_ch1.1", "characteristics_ch1.2")],
            by.x = "gsm", by.y = "geo_accession")
names(x)[names(x) == "characteristics_ch1.1"] <- "grade"
names(x)[names(x) == "characteristics_ch1.2"] <- "stage"
tum <- x[x$tissue == "Tumor" & !grepl("^[0-9]+$", x$patient), ]
tp  <- x$patient[x$tissue == "Tumor"]
nor <- x[x$tissue == "Normal" & grepl("^[0-9]+$", x$patient) & !x$patient %in% tp, ]
tum$num <- as.integer(sub("^E?([0-9]+).*$", "\\1", tum$patient))
res <- do.call(rbind, lapply(seq_len(nrow(tum)), function(i) {
  c <- nor[abs(as.integer(nor$patient) - tum$num[i]) <= 1, ]
  if (nrow(c) == 0) return(data.frame(tumor = tum$title[i], t_grade = tum$grade[i], t_stage = tum$stage[i],
                                      normal = NA, n_grade = NA, n_stage = NA))
  data.frame(tumor = tum$title[i], t_grade = tum$grade[i], t_stage = tum$stage[i],
             normal = c$title, n_grade = c$grade, n_stage = c$stage)
}))
options(width = 200); print(res, row.names = FALSE)
