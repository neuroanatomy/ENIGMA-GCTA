library(msm)

script.dir <- {function() {
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
}}()

ukb_annex <- system(paste("cd", script.dir, "&& git rev-parse --show-toplevel"), intern=T)
src_dir <- file.path(ukb_annex, "src")
dataset <- file.path(ukb_annex, "data/derived/ukb-22110_maf001")
hsq_dir <- file.path(dataset, "06.hsq_FreeSurfer")

freesurfer2pheno_bis <- function(pheno) {
  pheno <- sub("sum\\.", "", tolower(pheno))
  pheno <- sub("(\\.proper|\\.area)", "", pheno)
  pheno <- sub("brainseg", "brain", pheno)
  pheno <- sub("estimatedtotalintracranialvol", "icv", pheno)
}

plot_hsq_bis <- function(h2, outpdf = NULL, pattern_pheno_ignore = 'left|right', pdfwidth = 15, pdfheight = 8) {
  h2 <- h2[grep(pattern_pheno_ignore, pheno, invert = T),]
  setnames(h2, "V(G)/Vp", "vgvp")
  setnames(h2, "V(G)/Vp_se", "vgvp_se")
  
  h2 <- as.data.frame(h2, stringsAsFactors = F)
  h2$pheno <- sub('log10', '', h2$pheno)
  h2$pheno <- sub('_freesurfer', '', h2$pheno)
  row.names(h2) <- h2$pheno
  last_rows <- c("brain", "ICV", "intelligence", "height")
  h2 <- rbind(h2[!row.names(h2) %in% last_rows, ], h2[last_rows, ])
  h2$pheno <- factor(h2$pheno, levels=h2$pheno)
  
  h2 <- h2[row.names(h2),]
  h2$pheno <- factor(h2$pheno, levels=h2$pheno)
  
  h2$vgvp_lcl <- h2$vgvp - 1.96 * h2$vgvp_se
  h2$vgvp_ucl <- h2$vgvp + 1.96 * h2$vgvp_se
  
  gg <- ggplot(h2, aes(x=pheno, y=vgvp)) + #, group=variable)) + 
    geom_bar(position=position_dodge2(1), stat="identity", 
             width = 0.5, fill = "aquamarine4",
             colour="black",  # Use black outlines,
             size=.3) +
    ylab('VG/VP' ) +
    xlab('Phenotype') + 
    #scale_fill_brewer(palette = 'Set3') +
    #scale_fill_hue(name="Supplement type", # Legend label, use darker colors
    #               breaks=c("OJ", "VC"),
    #               labels=c("Orange juice", "Ascorbic acid")) +
    #facet_wrap(~margin, nrow = 3)+ 
    theme_bw( ) +
    theme( axis.text.x = element_text(angle = 90, size = 11, vjust = 0.5, hjust = 1), 
           axis.title.y = element_text(size = 10),
           axis.text.y = element_text(size = 10),
           axis.text.x.top = element_text(size = 10),
           strip.text.x = element_text(size=11), 
           legend.text = element_text(size = 10) )
  
  gg <- gg + geom_errorbar(aes(ymin=vgvp_lcl, ymax=vgvp_ucl),
                           size=.2,    # Thinner lines
                           width=0.5, 
                           position=position_dodge2(1))
  if (!is.null(outpdf)) {
    pdf(outpdf, width = pdfwidth, height = pdfheight)
    print(gg)
    dev.off()
  } else {
    print(gg)
  }
}

cor.file <- file.path(script.dir, "cor-fs-volumes.tsv")
cor.table <- read.csv(cor.file, sep="\t")
names(cor.table)[1] <- "phenotype"

# cor_cov.file <- file.path(script.dir, "cor-cov-fs-volumes.tsv")
# cor_cov.table <- read.csv(cor_cov.file, sep="\t")
# 
# cor_cov.table$pheno1 <- freesurfer2pheno_bis(cor_cov.table$Volume1)
# cor_cov.table$pheno2 <- freesurfer2pheno_bis(cor_cov.table$Volume2)
# 
# cor.table <- subset(cor_cov.table, pheno1==pheno2)
# # remove duplicate sum phenotypes
# cor.table <- subset(cor.table, !(Volume1 %in% c("Sum.CerebralWhiteMatter", "Sum.Cortex", "Sum.SurfaceHoles")))
# cor.table <- subset(cor.table, !(Volume2 %in% c("Sum.CerebralWhiteMatter", "Sum.Cortex", "Sum.SurfaceHoles")))
# cor.table <- data.frame(phenotype=cor.table$Volume1, pheno_bis=cor.table$pheno1, ICC=cor.table$ICC_cov)

# remove duplicate sum phenotypes
cor.table <- subset(cor.table, !(phenotype %in% c("Sum.CerebralWhiteMatter", "Sum.Cortex", "Sum.SurfaceHoles")))
cor.table$pheno_bis <- freesurfer2pheno_bis(cor.table$pheno)

row.names(cor.table) <- cor.table$pheno_bis
# cor.table$ICC.se <- (cor.table$ICC.ub - cor.table$ICC.lb) / (qnorm(0.975) * 2)

source(file.path(src_dir, "process_hsq_nopartition.R"))


groups <- c("hsq-all", "hsq-all_diff")

for (group in groups) {
  hsq.table <- readHsq(file.path(hsq_dir, group))
  pheno_bis <- sub("_freesurfer", "", tolower(hsq.table$pheno))
  pheno_bis <- sub("log10", "", pheno_bis)
  pheno_bis <- sub("(.*)left$", "left.\\1", pheno_bis)
  pheno_bis <- sub("(.*)right$", "right.\\1", pheno_bis)
  pheno_bis <- sub("(.*)_diff$", "diff.\\1", pheno_bis)
  
  # adjust estimated values for ICC
  ICC <- ifelse(pheno_bis == "height", 1, ifelse(pheno_bis == "intelligence", 0.63, cor.table[pheno_bis, "ICC"]))
  hsq.table.corrected <- hsq.table
  hsq.table.corrected$Vp <- hsq.table$Vp * ICC
  hsq.table.corrected$Vp_se <- hsq.table$Vp_se * ICC
  hsq.table.corrected$`V(G)/Vp` <- hsq.table$`V(G)/Vp` / ICC
  hsq.table.corrected$`V(G)/Vp_se` <- hsq.table$`V(G)/Vp_se` / ICC
  hsq.table.corrected$`V(e)` <- hsq.table$`V(e)` + (ICC - 1) * hsq.table$Vp
  # var(Vg) = var(Vp) + Var(Ve) - 2*Cov(Vp, Ve)
  hsq.table.corrected$`V(e)_se` <- sapply(1:nrow(hsq.table), function(i) {
    cor.row <- hsq.table[i,]
    cov_ve_vp <- (cor.row$Vp_se^2 + cor.row$`V(e)_se`^2 - cor.row$`V(G)_se`^2) / 2
    cov_matrix <- matrix(c(cor.row$`V(e)_se`^2, cov_ve_vp, cov_ve_vp, cor.row$Vp_se^2), nrow=2)
    cov_matrix <- cbind(rbind(cov_matrix, 0), 0)
    deltamethod(~x1 + (x3-1) * x2, c(cor.row$`V(e)`, cor.row$Vp, ICC[i]), cov_matrix)
  })

  corrected.file <- file.path(script.dir, paste0(group, "_corrected.tsv"))
  write.table(hsq.table.corrected, corrected.file, sep="\t", row.names=F)
  
  hsq.table.corrected <- hsq.table.corrected[grep('left|right', pheno, invert = T),]
  hsq.table.corrected[, pheno := gsub('log10|_freesurfer', '', pheno)]
  if (all(hsq.table.corrected$pheno %in% names(fsnames))) {
    fsnames <- c("accumbens" = 'Acc',
                 "amygdala" = 'Amy',
                 "putamen" = 'Pu',
                 "pallidum" = 'Pa',
                 "caudate" = 'Ca',
                 "thalamus" = 'Th',
                 "hippocampus" = 'Hip',
                 'brain' = 'BV',
                 'ICV' = 'ICV',
                 'height' = 'Height',
                 'intelligence' = 'FI')
    hsq.table.corrected[, pheno := fsnames[pheno]]
    hsq.table.corrected[, pheno := factor(pheno, levels = fsnames)]
  }
  pdf.file <- corrected.file <- file.path(script.dir, paste0(group, "_corrected.pdf"))
  plot_hsq(hsq.table.corrected, outpdf = pdf.file, pdfwidth = 7, pdfheight = 4.5)
  # plot_hsq_bis(hsq.table.corrected, outpdf=pdf.file)
}
