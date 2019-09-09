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

source(file.path(src_dir, "process_bivariate.R"))

correct_bivariate <- function(path_res, prefix_res, cor_table, cor_cov_table, output_file = NULL) {
  
  #### ==== Read original GCTA results  ==== ###
  res <- readHsqForBiv(path_res = path_res, prefix_res = prefix_res)
  
  res.zscores <- sapply(res, function(r) {
    
    #print(r$partition)
    nbsamples <- as.numeric(r$nbsamples)
    
    # find phenotypes in cor_table
    phenos_bis <- sub("_freesurfer", "", tolower(r$phenos))
    phenos_bis <- sub("log10", "", phenos_bis)
    phenos_bis <- sub("(.*)left$", "left.\\1", phenos_bis)
    phenos_bis <- sub("(.*)right$", "right.\\1", phenos_bis)
    
    cor.row1 <- cor_table[cor_table$pheno_bis == phenos_bis[1],]
    cor.row2 <- cor_table[cor_table$pheno_bis == phenos_bis[2],]
    cor_cov.row <- cor_cov_table[cor_cov_table$pheno1 == phenos_bis[1] & cor_cov_table$pheno2 == phenos_bis[2],]
    
    # give ICC of 1 for height and 0.63 for functional intelligence
    ICC1 <- ifelse(phenos_bis[1] == "height", 1, ifelse(phenos_bis[1] == "intelligence", 0.63, cor.row1$ICC))
    ICC2 <- ifelse(phenos_bis[2] == "height", 1, ifelse(phenos_bis[2] == "intelligence", 0.63, cor.row2$ICC))
    ICC_cov <- ifelse(any(phenos_bis %in% c("height", "intelligence")), 1, cor_cov.row$ICC_cov)
    
    # melted summary table
    
    # heritability values
    sum.melt <- melt(r$summary_table, id.vars=1)
    sum.melt <- structure(sum.melt$value, names=paste(sum.melt$Source, sum.melt$variable, sep="_"))
    names(sum.melt) <- sub("_Variance", "", names(sum.melt))
    # h <- h[grep("/Vp", names(h))]
    
    # variance and covariance explained by genetic and environment
    v <- c(sum.melt[1:6], ICC_tr1=ICC1, ICC_tr2=ICC2, ICC_cov=ICC_cov)
    w <- rbind(cbind(r$cov_matrix[1:6, 1:6], 0, 0, 0), 0, 0, 0)
    # w[7, 7] <- cor.row1$ICC.se^2
    # w[8, 8] <- cor.row2$ICC.se^2
    w[7, 7] <- 0
    w[8, 8] <- 0
    w[9, 9] <- 0
    
    w2 <- deltamethod(list(~x1, ~x2, ~x3, ~x7*x4+(x7-1)*x1, ~x8*x5+(x8-1)*x2, ~x6/sqrt(x7*x8)*x9), v, w, ses=F)
    v2 <- c(v[1], v[2], v[3], v[4] * v[7] + (v[7] - 1) * v[1], v[5] * v[8] + (v[8] - 1) * v[2], v[6] * v[9] + (v[9] - 1) * v[3])
    
    # compute correlations
    rg <- v2[3] / sqrt(v2[1] * v2[2])
    re <- v2[6] / sqrt(v2[4] * v2[5])
    rp <- (v2[6]+v2[3]) / sqrt((v2[1]+v2[4]) * (v2[2]+v2[5]))
    rs <- structure(c(rg, re, rp), names=c("rG", "rE", "rP"))
    rs_v <- deltamethod(list(~ x3 / sqrt(x1 * x2),
                             ~ x6 / sqrt(x4 * x5),
                             ~ (x6 + x3) / sqrt((x1+x4) * (x2+x5))),
                        v2, w2, ses=F)
    
    rs_se <- structure(sqrt(diag(rs_v)), names=paste0(names(rs), "_SE"))
    
    # difference between rE and rG correlations
    dr_rgre <- -diff(rs[1:2])
    dr_rgre_se <- deltamethod(~ x1 - x2, rs, rs_v)
    p.dr_rgre <- min(pnorm(dr_rgre, 0, dr_rgre_se, lower.tail = F),
                     pnorm(dr_rgre, 0, dr_rgre_se, lower.tail = T)) * 2
    
    # difference between rP and rG correlations
    dr_rgrp <- -diff(rs[c(1,3)])
    dr_rgrp_se <- deltamethod(~ x1 - x3, rs, rs_v)
    p.dr_rgrp <- min(pnorm(dr_rgrp, 0, dr_rgrp_se, lower.tail = F),
                     pnorm(dr_rgrp, 0, dr_rgrp_se, lower.tail = T)) * 2
    
    return(c(list(phenotype1=r$phenos[1], phenotype2=r$phenos[2]), rs, rs_se,
             dr_rgre=dr_rgre, dr_rgre_SE=dr_rgre_se, p.dr_rgre=p.dr_rgre,
             dr_rgrp=dr_rgrp, dr_rgrp_SE=dr_rgrp_se, p.dr_rgrp=p.dr_rgrp, 
             nind = nbsamples))
    
  }, simplify=F)
  
  pvals <- rbindlist(res.zscores)
  
  if (!is.null(output_file))
    fwrite(pvals, file = output_file, sep = '\t', row.names = F)
  
  return(pvals)
}

freesurfer2pheno_bis <- function(pheno) {
  pheno <- sub("sum\\.", "", tolower(pheno))
  pheno <- sub("(\\.proper|\\.area)", "", pheno)
  pheno <- sub("brainseg", "brain", pheno)
  pheno <- sub("estimatedtotalintracranialvol", "icv", pheno)
}

# cor.file <- file.path(script.dir, "cor-fs-volumes.tsv")
# cor.table <- read.csv(cor.file, sep="\t")
# names(cor.table)[1] <- "phenotype"
# # cor.table$ICC.se <- (cor.table$ICC.ub - cor.table$ICC.lb) / (qnorm(0.975) * 2)
# 
# cor.table$pheno_bis <- freesurfer2pheno_bis(cor.table$pheno)

cor_cov.file <- file.path(script.dir, "cor-cov-fs-volumes.tsv")
cor_cov.table <- read.csv(cor_cov.file, sep="\t")


cor_cov.table$pheno1 <- freesurfer2pheno_bis(cor_cov.table$Volume1)
cor_cov.table$pheno2 <- freesurfer2pheno_bis(cor_cov.table$Volume2)

cor.table <- subset(cor_cov.table, pheno1==pheno2)
cor.table <- data.frame(phenotype=cor.table$Volume1, pheno_bis=cor.table$pheno1, ICC=cor.table$ICC_cov)


res_cor <- correct_bivariate(path_res=file.path(hsq_dir, 'hsq-biv'),
                             prefix_res="all\\.rg\\\\?=0\\..*.hsq",
                             cor_table=cor.table,
                             cor_cov_table=cor_cov.table,
                             output_file = file.path(script.dir, paste0("summary_bivariate_corrected.tsv")))

# res_cor <- subset(res_cor, phenotype1 != "intelligence" & phenotype2 != "intelligence")
# res_cor[phenotype1 != "intelligence" & phenotype2 != "intelligence", 1:15] = 0

plot_bivariate(res_cor, script.dir)

