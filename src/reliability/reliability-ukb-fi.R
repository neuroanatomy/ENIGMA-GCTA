library(meta)

script.dir <- {
  initial.options <- commandArgs(trailingOnly = FALSE)
  file.arg.name <- "--file="
  script.name <- sub(file.arg.name, "", initial.options[grep(file.arg.name, initial.options)])
  sourceDir <- getSrcDirectory(function(dummy) {dummy})
  if (length(script.name)) { # called from command
    (dirname(script.name))
  } else if (nchar(sourceDir)) { # called with source
    sourceDir
  } else if (rstudioapi::isAvailable()) { # called from RStudio
    dirname(rstudioapi::getSourceEditorContext()$path)
  } else getwd()
}

annex_dir <- system(paste("cd", script.dir, "&& git rev-parse --show-toplevel"), intern=T)
ukb_all.pheno_dir <- file.path(annex_dir, "data/derived/ukb-all/00.phenotype_raw")

fi_files <- dir(ukb_all.pheno_dir, "f.20016.*.0.txt", full.names=T)

fi_tables <- sapply(fi_files, read.table, header=T, simplify=F)
fi_table <- Reduce(function(x,y) merge(x, y, all=T, by=c("FID", "IID")), fi_tables)

irr::icc(fi_table[, 3:5])
cor(fi_table[, 3:5], use="na.or.complete")
