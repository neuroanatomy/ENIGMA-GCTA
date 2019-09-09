library(foreach)
library(ggplot2)
library(psychometric)

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

ukb_annex <- system(paste("cd", script.dir, "&& git rev-parse --show-toplevel"), intern=T)
dataset <- file.path(ukb_annex, "data/derived/ukb-24220_maf001")
result_dir <- file.path(dataset, "PRSice/all-sans.ukb-all_maf001")
# result_dir <- file.path(dataset, "PRSice/all.ukb-22110_maf001")
pheno_dir <- file.path(dataset, "00.phenotype")

quant_covar <- c(file.path(pheno_dir, "age.txt"), file.path(dataset, "04.grm/grm-all-0.025/all-0.025.pca.eigenvec"))
qual_covar <- c(file.path(pheno_dir, "sex.txt"), file.path(pheno_dir, "centre.txt"))

quant_tables <- sapply(quant_covar, read.table, header=T)
qual_tables <- sapply(qual_covar, function(x) {
  y <- read.table(x, header=T)
  for (i in 3:ncol(y)) {
    y[,i] <- as.factor(y[,i])
  }
  y
}, simplify=F)

cov_tables <- c(quant_tables, qual_tables)
cov_table <- Reduce(function(x,y) merge(x, y, all=T, by=c("FID", "IID")), cov_tables)
cov_table <- na.omit(cov_table)

phenos <- sapply(dir(pheno_dir), tools::file_path_sans_ext)
phenos <- intersect(phenos, dir(result_dir))

phenosel <- phenos
phenosel <- phenosel[grep('left|right',phenosel, invert = TRUE)]
phenosel <- phenosel[phenosel %in% c("intelligence", "height") | grepl('_freesurfer', phenosel)]

last_ind <- unlist(sapply(c("brain", "ICV", "intelligence", "height"), grep, phenosel))
phenosel <- c(phenosel[-last_ind], phenosel[last_ind])

pdf('gps_plots.pdf', 8.27*1.2, 11.69*1.2)
par(mfrow=c(6,4))
r2_dt <- foreach (pheno=phenosel, .combine=rbind) %do% {
  
  print(pheno)
  
  pheno_file <- file.path(pheno_dir, paste0(pheno, ".txt"))
  prs_file <- file.path(result_dir, pheno, paste0(pheno, ".best"))
  
  prs_table <- read.table(prs_file, header=T)
  # fam_table <- read.table(fam_file, col.names=c("FID", "IID", "PID", "MID", "Sex", "Status"), na.strings=-9)
  pheno_table <- read.table(pheno_file, header=T)
  pheno_var <- names(pheno_table)[3]
  # mds_table <- read.table("covariates.tsv", header=T)
  
  if (grepl('log10', pheno_var)) {
    pheno_table[,3] = 10^pheno_table[,3] / 1000
  }
  
  table <- merge(prs_table, pheno_table, by=c("FID", "IID"))
  table <- merge(table, cov_table, by=c("FID", "IID"))
  
  table <- subset(table, In_Regression=="Yes")
  
  # model 1: compute r^2 on residual~PRS
  
  cov_cols <- c(sapply(cov_tables, function(x) colnames(x)[3:ncol(x)]), gps="PRS")
  r2.cov <- sapply(cov_cols, function(x) {
    mod_cov <- lm(as.formula(paste(pheno_var, "~", paste(x, collapse='+'))), data=table)
    summary(mod_cov)$r.squared
  })
  
  names(r2.cov) <- sapply(sapply(names(r2.cov), basename), tools::file_path_sans_ext)
  names(r2.cov) <- sub(".*\\.", "", names(r2.cov))
  names(r2.cov) <- paste0("Rsq_", names(r2.cov), "_only")
  
  mod0 <- lm(as.formula(paste(pheno_var, "~", paste(colnames(cov_table)[3:ncol(cov_table)], collapse=' + '))), data=table)
  table$resid <- mod0$residuals
  
  mod1 <- lm(as.formula(paste("resid", "~ PRS")), data=table)
  
  summary(mod1)
  rsq <- summary(mod1)$r.squared
  rsq.adj <- summary(mod1)$adj.r.squared
  rsq.p <- do.call(pf, c(unname(as.list(summary(mod1)$fstatistic)), lower.tail=F))
  
  mean.resid <- mean(table$resid)
  sd.resid <- sd(table$resid)
  
  mean.prs <- mean(table$PRS)
  sd.prs <- sd(table$PRS)
  rsqbis = signif(rsq,2)
  phenobis <- gsub('log10', '',pheno_var)
  rsq.pbis <- signif(rsq.p, 2)
  if (grepl("height", pheno_var)) {
    unit = "~(cm)"
  } else if (grepl("intelligence", pheno_var)) {
    unit = ""
  } else {
    unit = "~(cm^3)"
  }
  plot(table$PRS, table$resid,
       xlab=NA, ylab=NA,
       col=rgb(0,0,.5,alpha=.2), pch=16,
       main = phenobis)
  title(xlab = "GPS", ylab=parse(text = paste(phenobis, "~residual", unit)), line = 2.5)
  title(bquote( list(R^2 == .(rsqbis), p == .(rsq.pbis))), line = 1)
  
  abline(mod1$coefficients, col='red', lwd=1)
  abline(v=mean.prs + -2:2 * sd.prs, lty=2, lwd=1)
  abline(h=mean.resid + -2:2 * sd.resid, lty=2, lwd=1)
  

  low.resid <- which(table$resid < mean.resid - 2 * sd.resid)
  high.resid <- which(table$resid > mean.resid + 2 * sd.resid)
   
  low.prs <- which(table$PRS < mean.prs - 2 * sd.prs)
  high.prs <- which(table$PRS > mean.prs + 2 * sd.prs)
    
  if (length(low.prs)>0 & length(high.prs)>0) {
      boxplot(list(low_gps = table$resid[low.prs], high_gps = table$resid[high.prs]), 
              ylab=NA,
              main=gsub('log10', '', pheno_var))
    title(ylab=parse(text = paste(phenobis, "~residual", unit)), line=2.5)
    pval=signif(t.test(table$resid[low.prs], table$resid[high.prs])$p.value,2)
    st=signif(t.test(table$resid[low.prs], table$resid[high.prs])$statistic,2)
    title(parse(text = paste("list(t ==", st, ", p ==", pval, ")")), line = 1)
  }
  
  # model 2: compute partial epsilon squared for effect of PRS on full model
  
  mod2 <- lm(as.formula(paste(pheno_var, "~ PRS +", paste(colnames(cov_table)[3:ncol(cov_table)], collapse=' + '))), data=table)
  
  summary(mod2)
  
  drop1.mod2 <- drop1(mod2, test='F')
  SSB <- drop1.mod2["PRS", "Sum of Sq"]
  SST <- sum((mod2$model[[1]] - mean(mod2$model[[1]]))^2)
  SSE <- drop1.mod2[1, "RSS"]
  MSE <- SSE / mod2$df.residual
  DFB <- drop1.mod2["PRS", "Df"]
  
  etasq <- SSB / SST
  petasq <- SSB / (SSE + SSB)
  
  epsilonsq <- (SSB - DFB * MSE) / SST
  pepsilonsq <- (SSB - DFB * MSE) / (SSE + SSB)
  
  
  c(pheno=pheno, mean_pheno=mean(table[[pheno_var]]), sd_pheno=sd(table[[pheno_var]]),
    r2.cov, sd_resid=sd.resid, unlist(CI.Rsqlm(mod1)[1,,drop=TRUE]))
  
}
par(mfrow=c(1,1))
dev.off()

r2_dt <- as.data.frame(r2_dt, stringsAsFactors = F)
r2_dt$pheno <- sub('log10', '', r2_dt$pheno)
r2_dt$pheno <- sub('_freesurfer', '', r2_dt$pheno)
row.names(r2_dt) <- r2_dt$pheno
last_rows <- c("brain", "ICV", "intelligence", "height")
r2_dt <- rbind(r2_dt[!row.names(r2_dt) %in% last_rows, ], r2_dt[last_rows, ])
r2_dt$pheno <- factor(r2_dt$pheno, levels=r2_dt$pheno)

gg <- ggplot(r2_dt, 
       aes(x = pheno, y = as.numeric(Rsq))) + 
  geom_bar(stat='identity', 
           position=position_dodge2(1),
           size=.3,
           colour="black", 
           width = 0.5, fill = "cornflowerblue") +
    ylab( expression('R'^2)) +
  xlab('Phenotype') + 
    coord_cartesian(ylim = c(0, max(as.numeric(r2_dt$UCL)))) + 
  theme_bw( ) +
  theme( axis.text.x = element_text(angle = 90, size = 11, vjust = 0.5, hjust = 1), 
         axis.title.y = element_text(size = 10),
         axis.text.y = element_text(size = 10),
         axis.text.x.top = element_text(size = 10),
         strip.text.x = element_text(size=11), 
         legend.text = element_text(size = 10))


gg <- gg + geom_errorbar(aes(ymin=as.numeric(r2_dt$Rsq)-as.numeric(r2_dt$SErsq), ymax=as.numeric(r2_dt$Rsq)+as.numeric(r2_dt$SErsq)),
                         size=.2,    # Thinner lines
                         width=0.5, 
                         position=position_dodge2(1))

print(gg)

source(file.path(ukb_annex, 'src/process_hsq_nopartition.R'))
h2 = readHsq(file.path(ukb_annex, "data/derived/ukb-22110_maf001/06.hsq_FreeSurfer/hsq-all/"))
plot_hsq(h2)

setnames(h2, "V(G)/Vp", "vgvp")
setnames(h2, "V(G)/Vp_se", "vgvp_se")

h2 <- as.data.frame(h2, stringsAsFactors = F)
h2$pheno <- sub('log10', '', h2$pheno)
h2$pheno <- sub('_freesurfer', '', h2$pheno)
row.names(h2) <- h2$pheno
h2 <- h2[row.names(r2_dt),]
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

print(gg)

WriteXLS::WriteXLS(x = r2_dt,
                   ExcelFileName = 'gps_ukb.xls', row.names = F)

