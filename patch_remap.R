f <- "geo_limma_pdac.R"
x <- readLines(f)
i <- grep("GSE62452 = list", x)
x[i] <- '  GSE62452 = list(mode = "dupcor", patient_col = "title", patient_rx = "^E0*([0-9]+)-?.*$", text_cols = NULL, patient_remap = c("9728T"="9727", "17531T"="17530", "17798T"="17797", "19927T"="19926", "22506T"="22505", "26172T"="26171", "27881T"="27880", "29164T"="29163", "29818T"="29819", "29827T"="29826", "31263T"="31262", "31393T"="31392", "32018T"="32019", "34245T"="34244", "34435T"="34434")),'
j <- grep('patient <- sub("^[^:]*:', x, fixed = TRUE)
stopifnot(length(i) == 1, length(j) == 1)
add <- c('  if (!is.null(cfg$patient_remap)) {',
         '    hit <- patient %in% names(cfg$patient_remap)',
         '    patient[hit] <- unname(cfg$patient_remap[patient[hit]])',
         '  }')
x <- append(x, add, after = j)
writeLines(x, f)
