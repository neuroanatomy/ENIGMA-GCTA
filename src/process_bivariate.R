# This script plots and writes in a table the genetic correlation estimates
# for a given 

library(tools)
library(data.table)
# library(foreach)
library(msm)
library(ggplot2)
library(psych)
library(corrplot)

readHsqForBiv <- function(path_res, prefix_res = '.hsq') {

  fres <- list.files(path_res, pattern = prefix_res , full.names = T)

  # load files
  res <- sapply(fres, function(f) {
    print(f)
    # read log file instead of hsq file
    f2 <- paste0(file_path_sans_ext(f), '.log')
    lines <- readLines(f2)
    line_summary <- grep("Summary result of REML analysis:", lines)
    line_variance <- grep("Sampling variance/covariance", lines)
    line_hsq <- grep("Summary result of REML analysis has been saved", lines)
    line_n <- grep("individuals are in common in these files", lines)
    
    summary_table <- read.delim(text=lines[c((line_summary+1):(line_variance-1))], stringsAsFactors=F)
    nbsamples <-  sub(' .*', '', lines[line_n])

    pval <- sub("Pval\t(\\S+)\\s.*", "\\1", grep("^Pval", readLines(paste0(file_path_sans_ext(f), '.hsq')), value=T))
    
    cov_matrix <- read.delim(text=lines[(line_variance+1):(line_hsq-1)], header=F)
    cov_matrix <- as.matrix(cov_matrix[,1:nrow(cov_matrix)])
    
    phenos=c(sub("^.*\\.([^.]+)(\\.[^.]+){2}$", "\\1", f), sub("^.*\\.([^.]+)\\.[^.]+$", "\\1", f))
    
    return(list(summary_table=summary_table, cov_matrix=cov_matrix, phenos=phenos, nbsamples = nbsamples, Pval = pval))
    
  }, simplify = F)
  return(res)
}

process_bivariate <- function(path_res, prefix_res, output_file = NULL) {

  #### ==== Read original GCTA results  ==== ###
  res <- readHsqForBiv(path_res = path_res, prefix_res = prefix_res)

  res.zscores <- sapply(res, function(r) {
    
    #print(r$partition)
    nbsamples <- r$nbsamples
    pvallrt <- r$Pval
    
    # melted summary table
    
    # heritability values
    sum.melt=melt(r$summary_table, id.vars=1)
    sum.melt <- structure(sum.melt$value, names=paste(sum.melt$Source, sum.melt$variable, sep="_"))
    names(sum.melt) <- sub("_Variance", "", names(sum.melt))
    # h <- h[grep("/Vp", names(h))]
    
    # variance and covariance explained by genetic and environment
    v <- c(sum.melt[1:6])
    w <- r$cov_matrix[1:6, 1:6]
    
    # compute correlations
    rg <- v[3] / sqrt(v[1] * v[2])
    re <- v[6] / sqrt(v[4] * v[5])
    rp <- (v[6]+v[3]) / sqrt((v[1]+v[4]) * (v[2]+v[5]))
    rs <- structure(c(rg, re, rp), names=c("rG", "rE", "rP"))
    rs_v <- deltamethod(list(~ x3 / sqrt(x1 * x2), 
                             ~ x6 / sqrt(x4 * x5),
                             ~ (x6 + x3) / sqrt((x1+x4) * (x2+x5))), 
                        v, w, ses=F)
    
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
             nind = nbsamples, Plrt = pvallrt))
    
  }, simplify=F)
  
  pvals <- rbindlist(res.zscores)

  if (!is.null(output_file))
    fwrite(pvals, file = output_file, sep = '\t', row.names = F)

  return(pvals)
}

plot_bivariate <- function(res, dirout = '') {
  
  if (is.character(res)) {
    res <- fread(res, sep = '\t', stringsAsFactors = F)
  }
  
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
  
  lrfsnames <- c("accumbensleft" = 'Acc_L', "accumbensright" = 'Acc_R', 
                 "amygdalaleft" = 'Amy_L', "amygdalaright" = 'Amy_R', 
                 "putamenleft" = 'Pu_L', "putamenright" = 'Pu_R', 
                 "pallidumleft" = 'Pa_L', "pallidumright"  = 'Pa_R', 
                 "caudateleft" = 'Ca_L',  "caudateright"  = 'Ca_R', 
                 "thalamusleft" = 'Th_L', "thalamusright" = 'Th_R',
                 "hippocampusleft" = 'Hip_L', "hippocampusright" = 'Hip_R', 
                 'brain' = 'BV',
                 'ICV' = 'ICV',
                 'height' = 'Height',
                 'intelligence' = 'FI')
  
  
  ## add p-value for rP 
  res$p.rp = signif(psych::r.test(res$nind, res$rP, twotailed = T)$p, 4)
  
  ## ==== Scatter plots rP vs rG ====
  
  lm_eqn <- 
    function(df){ #source http://goo.gl/K4yh
      m <- lm(y ~ x, df);
      eq <- substitute(italic(y) == a + b %.% italic(x)*","~~italic(r)^2~"="~r2, 
                       list(a = format(coef(m)[1], digits = 2), 
                            b = format(coef(m)[2], digits = 2), 
                            r2 = format(summary(m)$r.squared, digits = 3)))
      as.character(as.expression(eq));
    }
  # select pairs of phenotypes only once
  res_un <- res[-(which(duplicated((mapply(x=res$phenotype1, y=res$phenotype2, function(x,y) paste(sort(c(x,y))), SIMPLIFY = F))))),]
  res_un <- res_un[grep('left|right', phenotype1, invert = T),]
  res_un <- res_un[grep('left|right', phenotype2, invert = T),]
  res_un <- res_un[grep('_noCovar', phenotype2, invert = T),]
  res_un <- res_un[grep('_noCovar', phenotype1, invert = T),]
  
  gg <- ggplot(res_un, aes(x = rG, y  = rP)) + 
    geom_point(alpha = 0.4, size = 5) + 
    theme_bw() + 
    coord_cartesian(ylim = c(0,0.9), xlim = c(0, 0.9)) +
    geom_smooth(method=lm) +
    geom_text(x = 0.4, y = 0.75,
              label = lm_eqn(data.frame(x = res$rG,y = res$rP)),
              parse = TRUE, size = 4)
  pdf(file.path(dirout,'rP_vs_rG.pdf'), 5, 5)
  print(gg)
  dev.off()
  
  ## ==== Scatter plots rE vs rG ====
  
  # select pairs of phenotypes only once
  res_un <- res[-(which(duplicated((mapply(x=res$phenotype1, y=res$phenotype2, function(x,y) paste(sort(c(x,y))), SIMPLIFY = F))))),]
  res_un <- res_un[grep('left|right', phenotype1, invert = T),]
  res_un <- res_un[grep('left|right', phenotype2, invert = T),]
  res_un <- res_un[grep('_noCovar', phenotype2, invert = T),]
  res_un <- res_un[grep('_noCovar', phenotype1, invert = T),]
  
  gg <- ggplot(res_un, aes(x = rG, y  = rE)) + 
    geom_point(alpha = 0.4, size = 5) + 
    theme_bw() + 
    coord_cartesian(ylim = c(0,0.9), xlim = c(0, 0.9)) +
    geom_smooth(method=lm) +
    geom_text(x = 0.4, y = 0.75,
              label = lm_eqn(data.frame(x = res$rG,y = res$rE)),
              parse = TRUE, size = 4)
  pdf(file.path(dirout,'rE_vs_rG.pdf'), 5, 5)
  print(gg)
  dev.off()
  
  
  
  ## ==== corrplot rE ===
  
  res[, phenotype1 := gsub('log10|_freesurfer', '', phenotype1)]
  res[, phenotype2 := gsub('log10|_freesurfer', '', phenotype2)]
  res <- res[!(phenotype1 %in% 'height_noCovariates' | phenotype2 %in% 'height_noCovariates' ), ]
  
  
  res2 <- copy(res)
  res2[, phenotype1 := res$phenotype2]
  res2[, phenotype2 := res$phenotype1]
  res2 <- rbind(res,res2)
  res2 <- unique(res2)
  #res2 %>% reshape2::dcast(phenotype1 ~ phenotype2, value.var = 'rE') 
  
  res2 <- res2[grep('right|left', phenotype1, invert = T),]
  res2 <- res2[grep('right|left', phenotype2, invert = T),]
  res2 <- res2[grep('_noCovar', phenotype1, invert = T),]
  res2 <- res2[grep('_noCovar', phenotype2, invert = T),]
  
  res2[, fdr.dr_rgre := p.adjust(p.dr_rgre, method = 'BY')]
  res2[, fdr.dr_rgrp := p.adjust(p.dr_rgrp, method = 'BY')]
  
  res2[, phenotype1 := fsnames[phenotype1]]
  res2[, phenotype2 := fsnames[phenotype2]]
  
  m <- matrix(nrow = length(unique(res2$phenotype1)), ncol = length(unique(res2$phenotype2)))
  rownames(m) <- unique(res2$phenotype1)
  colnames(m) <- unique(res2$phenotype2)
  apply(res2, 1, function(x) m[x[1],x[2]] <<- m[x[2],x[1]] <<- as.numeric(x['rE']))
  mcor <- m[unique(res2$phenotype1), unique(res2$phenotype1)]
  
  m <- matrix(nrow = length(unique(res2$phenotype1)), ncol = length(unique(res2$phenotype2)))
  rownames(m) <- unique(res2$phenotype1)
  colnames(m) <- unique(res2$phenotype2)
  apply(res2, 1, function(x) m[x[1],x[2]] <<- m[x[2],x[1]] <<- as.numeric(x['rE_SE']))
  mse <- m[unique(res2$phenotype1), unique(res2$phenotype1)]
  
  # m <- matrix(nrow = length(unique(res2$phenotype1)), ncol = length(unique(res2$phenotype2)))
  # rownames(m) <- unique(res2$phenotype1)
  # colnames(m) <- unique(res2$phenotype2)
  # apply(res2, 1, function(x) m[x[1],x[2]] <<- m[x[2],x[1]] <<- as.numeric(x[11]))
  # mpval <- m[unique(res2$phenotype1), unique(res2$phenotype1)]
  
  # mpval <- mpval[fsnames,fsnames]
  mcor <- mcor[fsnames,fsnames]
  mse <- mse[fsnames,fsnames]
  mlci <- mcor - mse
  muci <- mcor + mse
  
  ## round correlation values to stay in [-1, 1] range
  mcor[mcor > 1] <- 1
  mcor[mcor < -1] <- -1
  mlci[mlci > 1] <- 1
  mlci[mlci < -1] <- -1
  muci[muci > 1] <- 1
  muci[muci < -1] <- -1
  
  pdf(file.path(dirout,'rE.pdf'), 8, 8)
  col3 <- colorRampPalette(c("red", "white", "blue")) 
  cp1 <- corrplot(mcor, col =  col3(200), order = 'original',
                  type = "lower", method = 'shade',#'number',
                  diag = TRUE, tl.pos = "n", #cl.pos = "n",
                  mar = c(0,0,4,0),
                  addCoef.col = "black", number.cex = .8)
  cp2 <- corrplot(mcor, order = 'original', # p.mat = mpval,
                  type = "upper", method = 'circle', diag = T, col =  col3(200),
                  tl.pos = 'd', tl.cex = 0.8, plotCI = 'circle',  add = TRUE, cl.pos = 'n',
                  #addCoef.col = "black", number.cex = .7, #col = upper.col,
                  lowCI.mat = mlci, uppCI.mat = muci, number.cex = .8,
                  insig = 'label_sig', pch = 1, pch.cex = 1, sig.level = c(.001, .01, .05))
  dev.off()
  
  ## ==== corrplot rP ===
  
  m <- matrix(nrow = length(unique(res2$phenotype1)), ncol = length(unique(res2$phenotype2)))
  rownames(m) <- unique(res2$phenotype1)
  colnames(m) <- unique(res2$phenotype2)
  apply(res2, 1, function(x) m[x[1],x[2]] <<- m[x[2],x[1]] <<- as.numeric(x['rP']))
  mcor <- m[unique(res2$phenotype1), unique(res2$phenotype1)]
  
  m <- matrix(nrow = length(unique(res2$phenotype1)), ncol = length(unique(res2$phenotype2)))
  rownames(m) <- unique(res2$phenotype1)
  colnames(m) <- unique(res2$phenotype2)
  apply(res2, 1, function(x) m[x[1],x[2]] <<- m[x[2],x[1]] <<- as.numeric(x['rP_SE']))
  mse <- m[unique(res2$phenotype1), unique(res2$phenotype1)]
  
  m <- matrix(nrow = length(unique(res2$phenotype1)), ncol = length(unique(res2$phenotype2)))
  rownames(m) <- unique(res2$phenotype1)
  colnames(m) <- unique(res2$phenotype2)
  apply(res2, 1, function(x) m[x[1],x[2]] <<- m[x[2],x[1]] <<- as.numeric(x['p.rp']))
  mpval <- m[unique(res2$phenotype1), unique(res2$phenotype1)]
  
  mpval <- mpval[fsnames,fsnames]
  mcor <- mcor[fsnames,fsnames]
  mse <- mse[fsnames,fsnames]
  mlci <- mcor - mse
  muci <- mcor + mse
  
  ## round correlation values to stay in [-1, 1] range
  mcor[mcor > 1] <- 1
  mcor[mcor < -1] <- -1
  mlci[mlci > 1] <- 1
  mlci[mlci < -1] <- -1
  muci[muci > 1] <- 1
  muci[muci < -1] <- -1
  
  pdf(file.path(dirout,'rP.pdf'), 8, 8)
  col3 <- colorRampPalette(c("red", "white", "blue")) 
  cp1 <- corrplot(mcor, col =  col3(200), order = 'original',
                  type = "lower", method = 'shade',#'number',
                  diag = TRUE, tl.pos = "n", #cl.pos = "n",
                  mar = c(0,0,4,0),
                  addCoef.col = "black", number.cex = .8)
  cp2 <- corrplot(mcor, order = 'original',  p.mat = mpval,
                  type = "upper", method = 'circle', diag = T, col =  col3(200),
                  tl.pos = 'd', tl.cex = 0.8, plotCI = 'circle',  add = TRUE, cl.pos = 'n',
                  #addCoef.col = "black", number.cex = .7, #col = upper.col,
                  lowCI.mat = mlci, uppCI.mat = muci, number.cex = .8,
                  insig = 'label_sig', pch = 1, pch.cex = 1, sig.level = c(.001, .01, .05))
  dev.off()
  
  
  ## ==== corrplot rG ===
  
  m <- matrix(nrow = length(unique(res2$phenotype1)), ncol = length(unique(res2$phenotype2)))
  rownames(m) <- unique(res2$phenotype1)
  colnames(m) <- unique(res2$phenotype2)
  apply(res2, 1, function(x) m[x[1],x[2]] <<- m[x[2],x[1]] <<- as.numeric(x['rG']))
  mcor <- m[unique(res2$phenotype1), unique(res2$phenotype1)]
  
  m <- matrix(nrow = length(unique(res2$phenotype1)), ncol = length(unique(res2$phenotype2)))
  rownames(m) <- unique(res2$phenotype1)
  colnames(m) <- unique(res2$phenotype2)
  apply(res2, 1, function(x) m[x[1],x[2]] <<- m[x[2],x[1]] <<- as.numeric(x['rG_SE']))
  mse <- m[unique(res2$phenotype1), unique(res2$phenotype1)]
  
  m <- matrix(nrow = length(unique(res2$phenotype1)), ncol = length(unique(res2$phenotype2)))
  rownames(m) <- unique(res2$phenotype1)
  colnames(m) <- unique(res2$phenotype2)
  apply(res2, 1, function(x) m[x[1],x[2]] <<- m[x[2],x[1]] <<- as.numeric(x['Plrt']))
  mpval <- m[unique(res2$phenotype1), unique(res2$phenotype1)]
  
  mpval <- mpval[fsnames,fsnames]
  mcor <- mcor[fsnames,fsnames]
  mse <- mse[fsnames,fsnames]
  mlci <- mcor - mse
  muci <- mcor + mse
  
  ## round correlation values to stay in [-1, 1] range
  mcor[mcor > 1] <- 1
  mcor[mcor < -1] <- -1
  mlci[mlci > 1] <- 1
  mlci[mlci < -1] <- -1
  muci[muci > 1] <- 1
  muci[muci < -1] <- -1
  
  pdf(file.path(dirout,'rG.pdf'), 8, 8)
  col3 <- colorRampPalette(c("red", "white", "blue")) 
  cp1 <- corrplot(mcor, col =  col3(200), order = 'original',
                  type = "lower", method = 'shade',#'number',
                  diag = TRUE, tl.pos = "n", #cl.pos = "n",
                  mar = c(0,0,4,0),
                  addCoef.col = "black", number.cex = .8)
  cp2 <- corrplot(mcor, order = 'original',  p.mat = mpval,
                  type = "upper", method = 'circle', diag = T, col =  col3(200),
                  tl.pos = 'd', tl.cex = 0.8, plotCI = 'circle',  add = TRUE, cl.pos = 'n',
                  #addCoef.col = "black", number.cex = .7, #col = upper.col,
                  lowCI.mat = mlci, uppCI.mat = muci, number.cex = .8,
                  insig = 'label_sig', pch = 1, pch.cex = 1, sig.level = c(.001, .01, .05))
  dev.off()
  
  
  
  ## ==== corrplot rG-rE ===
  m <- matrix(nrow = length(unique(res2$phenotype1)), ncol = length(unique(res2$phenotype2)))
  rownames(m) <- unique(res2$phenotype1)
  colnames(m) <- unique(res2$phenotype2)
  apply(res2, 1, function(x) m[x[1],x[2]] <<- m[x[2],x[1]] <<- as.numeric(x["dr_rgre.rE"]))
  mcor <- m[unique(res2$phenotype1),unique(res2$phenotype1)]
  
  m <- matrix(nrow = length(unique(res2$phenotype1)), ncol = length(unique(res2$phenotype2)))
  rownames(m) <- unique(res2$phenotype1)
  colnames(m) <- unique(res2$phenotype2)
  apply(res2, 1, function(x) m[x[1],x[2]] <<- m[x[2],x[1]] <<- as.numeric(x["dr_rgre_SE"]))
  mse <- m[unique(res2$phenotype1),unique(res2$phenotype1)]
  
  m <- matrix(nrow = length(unique(res2$phenotype1)), ncol = length(unique(res2$phenotype2)))
  rownames(m) <- unique(res2$phenotype1)
  colnames(m) <- unique(res2$phenotype2)
  apply(res2, 1, function(x) m[x[1],x[2]] <<- m[x[2],x[1]] <<- as.numeric(x["fdr.dr_gre"]))
  mpval <- m[unique(res2$phenotype1),unique(res2$phenotype1)]
  
  mlci <- mcor - mse
  muci <- mcor + mse
  
  ## round correlation values to stay in [-1, 1] range
  mcor[mcor > 1] <- 1
  mcor[mcor < -1] <- -1
  mlci[mlci > 1] <- 1
  mlci[mlci < -1] <- -1
  muci[muci > 1] <- 1
  muci[muci < -1] <- -1
  
  pdf(file.path(dirout,'rG-rE.pdf'), 8, 8)
  col3 <- colorRampPalette(c("red", "white", "blue")) 
  cp1 <- corrplot(mcor, col =  col3(200),order = 'original',
                  type = "lower", method = 'shade',#'number',
                  diag = TRUE, tl.pos = "n", #cl.pos = "n",
                  mar = c(0,0,4,0),
                  addCoef.col = "black", number.cex = .8)
  cp2 <- corrplot(mcor, order = 'original',  p.mat = mpval,
                  type = "upper", method = 'circle', diag = T, col =  col3(200),
                  tl.pos = 'd', tl.cex = 0.8, plotCI = 'circle',  add = TRUE, cl.pos = 'n',
                  #addCoef.col = "black", number.cex = .7, #col = upper.col,
                  lowCI.mat = mlci, uppCI.mat = muci, number.cex = .8,
                  insig = 'label_sig', pch = 1, pch.cex = 1, sig.level = c(.001, .01, .05))
  dev.off()
  
  res3 <- copy(res)
  res3[, phenotype1 := res$phenotype2]
  res3[, phenotype2 := res$phenotype1]
  
  res3 <- rbind(res,res3)
  res3 <- unique(res3)
  res3[, fdr.dr_rgre := p.adjust(p.dr_rgre, method = 'BY')]
  res3[, fdr.dr_rgrp := p.adjust(p.dr_rgrp, method = 'BY')]
  
  #res3 %>% reshape2::dcast(phenotype1 ~ phenotype2, value.var = 'rE') 
  
  res3 <- res3[grep('right|left|height|ICV|brain|intelligence', phenotype1),]
  res3 <- res3[grep('right|left|height|ICV|brain|intelligence', phenotype2),]
  
  res3[, phenotype1 := lrfsnames[phenotype1]]
  res3[, phenotype2 := lrfsnames[phenotype2]]
  
  m <- matrix(nrow = length(unique(res3$phenotype1)), ncol = length(unique(res3$phenotype1)))
  rownames(m) <- unique(res3$phenotype1)
  colnames(m) <- unique(res3$phenotype1)
  apply(res3, 1, function(x) m[x[1],x[2]] <<- m[x[2],x[1]] <<- as.numeric(x["dr_rgre.rE"]))
  mcor <- m[unique(res3$phenotype1),unique(res3$phenotype1)]
  
  m <- matrix(nrow = length(unique(res3$phenotype1)), ncol = length(unique(res3$phenotype2)))
  rownames(m) <- unique(res3$phenotype1)
  colnames(m) <- unique(res3$phenotype1)
  apply(res3, 1, function(x) m[x[1],x[2]] <<- m[x[2],x[1]] <<- as.numeric(x["dr_rgre_SE"]))
  mse <- m[unique(res3$phenotype1),unique(res3$phenotype1)]
  
  m <- matrix(nrow = length(unique(res3$phenotype1)), ncol = length(unique(res3$phenotype2)))
  rownames(m) <- unique(res3$phenotype1)
  colnames(m) <- unique(res3$phenotype1)
  apply(res3, 1, function(x) m[x[1],x[2]] <<- m[x[2],x[1]] <<- as.numeric(x["fdr.dr_rgre"]))
  mpval <- m[unique(res3$phenotype1),unique(res3$phenotype1)]
  
  mpval <- mpval[rownames(mpval) %in% lrfsnames,rownames(mpval) %in% lrfsnames]
  mcor <- mcor[rownames(mcor) %in% lrfsnames,rownames(mcor) %in% lrfsnames]
  mse <- mse[rownames(mse) %in% lrfsnames,rownames(mse) %in% lrfsnames]
  mlci <- mcor - mse
  muci <- mcor + mse
  
  ## round correlation values to stay in [-1, 1] range
  mcor[mcor > 1] <- 1
  mcor[mcor < -1] <- -1
  mlci[mlci > 1] <- 1
  mlci[mlci < -1] <- -1
  muci[muci > 1] <- 1
  muci[muci < -1] <- -1
  
  pdf(file.path(dirout, 'rG-rE_leftright.pdf'), 9, 9)
  col3 <- colorRampPalette(c("red", "white", "blue")) 
  cp1 <- corrplot(mcor, col =  col3(200),order = 'original',
                  type = "lower", method = 'shade',#'number',
                  diag = TRUE, tl.pos = "n", #cl.pos = "n",
                  mar = c(0,0,4,0), title = 'rG - rE', tl.col = 'black',
                  addCoef.col = "black", number.cex = .8)
  cp2 <- corrplot(mcor, order = 'original',  p.mat = mpval, tl.col = 'black',
                  type = "upper", method = 'circle', diag = T, col =  col3(200),
                  tl.pos = 'd', tl.cex = 0.8, plotCI = 'circle',  add = TRUE, cl.pos = 'n',
                  #addCoef.col = "black", number.cex = .7, #col = upper.col,
                  lowCI.mat = mlci, uppCI.mat = muci, number.cex = .8,
                  insig = 'label_sig', pch = 1, pch.cex = 1, sig.level = c(.001, .01, .05))
  dev.off()
  
  if (sum(!is.na(mcor)) > ncol(mcor)) {
    gplots::heatmap.2(mcor, 
                      col=(col3(20)), scale = 'none',  
                      hclustfun = function(x) hclust(x, method = 'ward.D2'), margins = c(10,10),
                      main = 'rE', trace  = 'none')
  }
  
  
  ## ==== corrplot rG-rP ===
  m <- matrix(nrow = length(unique(res2$phenotype1)), ncol = length(unique(res2$phenotype2)))
  rownames(m) <- unique(res2$phenotype1)
  colnames(m) <- unique(res2$phenotype2)
  apply(res2, 1, function(x) m[x[1],x[2]] <<- m[x[2],x[1]] <<- as.numeric(x["dr_rgrp.rP"]))
  mcor <- m[unique(res2$phenotype1),unique(res2$phenotype1)]
  
  m <- matrix(nrow = length(unique(res2$phenotype1)), ncol = length(unique(res2$phenotype2)))
  rownames(m) <- unique(res2$phenotype1)
  colnames(m) <- unique(res2$phenotype2)
  apply(res2, 1, function(x) m[x[1],x[2]] <<- m[x[2],x[1]] <<- as.numeric(x["dr_rgrp_SE"]))
  mse <- m[unique(res2$phenotype1),unique(res2$phenotype1)]
  
  m <- matrix(nrow = length(unique(res2$phenotype1)), ncol = length(unique(res2$phenotype2)))
  rownames(m) <- unique(res2$phenotype1)
  colnames(m) <- unique(res2$phenotype2)
  apply(res2, 1, function(x) m[x[1],x[2]] <<- m[x[2],x[1]] <<- as.numeric(x["fdr.dr_rgrp"]))
  mpval <- m[unique(res2$phenotype1),unique(res2$phenotype1)]
  
  mpval <- mpval[fsnames,fsnames]
  mcor <- mcor[fsnames,fsnames]
  mse <- mse[fsnames,fsnames]
  mlci <- mcor - mse
  muci <- mcor + mse
  
  ## round correlation values to stay in [-1, 1] range
  mcor[mcor > 1] <- 1
  mcor[mcor < -1] <- -1
  mlci[mlci > 1] <- 1
  mlci[mlci < -1] <- -1
  muci[muci > 1] <- 1
  muci[muci < -1] <- -1
  
  pdf(file.path(dirout,'rG-rP.pdf'), 8, 8)
  col3 <- colorRampPalette(c("red", "white", "blue")) 
  cp1 <- corrplot(mcor, col =  col3(200),order = 'original',
                  type = "lower", method = 'shade',#'number',
                  diag = TRUE, tl.pos = "n", #cl.pos = "n",
                  mar = c(0,0,4,0),
                  addCoef.col = "black", number.cex = .8)
  cp2 <- corrplot(mcor, order = 'original',  p.mat = mpval,
                  type = "upper", method = 'circle', diag = T, col =  col3(200),
                  tl.pos = 'd', tl.cex = 0.8, plotCI = 'circle',  add = TRUE, cl.pos = 'n',
                  #addCoef.col = "black", number.cex = .7, #col = upper.col,
                  lowCI.mat = mlci, uppCI.mat = muci, number.cex = .8,
                  insig = 'label_sig', pch = 1, pch.cex = 1, sig.level = c(.001, .01, .05))
  dev.off()
  
  res3 <- copy(res)
  res3[, phenotype1 := res$phenotype2]
  res3[, phenotype2 := res$phenotype1]
  res3 <- rbind(res,res3)
  res3 <- unique(res3)
  #res3 %>% reshape2::dcast(phenotype1 ~ phenotype2, value.var = 'rE') 
  
  res3[, fdr.dr_rgre := p.adjust(p.dr_rgre, method = 'BY')]
  res3[, fdr.dr_rgrp := p.adjust(p.dr_rgrp, method = 'BY')]
  
  res3 <- res3[grep('right|left|height|ICV|brain|intelligence', phenotype1),]
  res3 <- res3[grep('right|left|height|ICV|brain|intelligence', phenotype2),]
  
  res3[, phenotype1 := lrfsnames[phenotype1]]
  res3[, phenotype2 := lrfsnames[phenotype2]]
  
  m <- matrix(nrow = length(unique(res3$phenotype1)), ncol = length(unique(res3$phenotype2)))
  rownames(m) <- unique(res3$phenotype1)
  colnames(m) <- unique(res3$phenotype2)
  apply(res3, 1, function(x) m[x[1],x[2]] <<- m[x[2],x[1]] <<- as.numeric(x["dr_rgrp.rP"]))
  mcor <- m[unique(res3$phenotype1),unique(res3$phenotype1)]
  
  m <- matrix(nrow = length(unique(res3$phenotype1)), ncol = length(unique(res3$phenotype2)))
  rownames(m) <- unique(res3$phenotype1)
  colnames(m) <- unique(res3$phenotype2)
  apply(res3, 1, function(x) m[x[1],x[2]] <<- m[x[2],x[1]] <<- as.numeric(x["dr_rgrp_SE"]))
  mse <- m[unique(res3$phenotype1),unique(res3$phenotype1)]
  
  m <- matrix(nrow = length(unique(res3$phenotype1)), ncol = length(unique(res3$phenotype2)))
  rownames(m) <- unique(res3$phenotype1)
  colnames(m) <- unique(res3$phenotype2)
  apply(res3, 1, function(x) m[x[1],x[2]] <<- m[x[2],x[1]] <<- as.numeric(x["fdr.dr_rgrp"]))
  mpval <- m[unique(res3$phenotype1),unique(res3$phenotype1)]
  
  mpval <- mpval[rownames(mpval) %in% lrfsnames,rownames(mpval) %in% lrfsnames]
  mcor <- mcor[rownames(mcor) %in% lrfsnames,rownames(mcor) %in% lrfsnames]
  mse <- mse[rownames(mse) %in% lrfsnames,rownames(mse) %in% lrfsnames]
  mlci <- mcor - mse
  muci <- mcor + mse
  
  ## round correlation values to stay in [-1, 1] range
  mcor[mcor > 1] <- 1
  mcor[mcor < -1] <- -1
  mlci[mlci > 1] <- 1
  mlci[mlci < -1] <- -1
  muci[muci > 1] <- 1
  muci[muci < -1] <- -1
  
  pdf(file.path(dirout,'rG-rP_leftright.pdf'), 9, 9)
  col3 <- colorRampPalette(c("red", "white", "blue")) 
  cp1 <- corrplot(mcor, col =  col3(200),order = 'original',
                  type = "lower", method = 'shade',#'number',
                  diag = TRUE, tl.pos = "n", #cl.pos = "n",
                  mar = c(0,0,4,0), tl.col = 'black',title = 'rG - rP',
                  addCoef.col = "black", number.cex = .8)
  cp2 <- corrplot(mcor, order = 'original',  p.mat = mpval,
                  type = "upper", method = 'circle', diag = T, col =  col3(200),
                  tl.pos = 'd', tl.cex = 0.8, plotCI = 'circle',  add = TRUE, cl.pos = 'n',
                  #addCoef.col = "black", number.cex = .7, #col = upper.col,
                  lowCI.mat = mlci, uppCI.mat = muci, number.cex = .8, tl.col = 'black',
                  insig = 'label_sig', pch = 1, pch.cex = 1, sig.level = c(.001, .01, .05))
  dev.off()
  
  if (sum(!is.na(mcor)) > ncol(mcor)) {
    gplots::heatmap.2(mcor, 
                      col=(col3(20)), scale = 'none',  
                      hclustfun = function(x) hclust(x, method = 'ward.D2'), margins = c(10,10),
                      main = 'rE', trace  = 'none')
  }
  
  
  # res[, Phenotype:= gsub('_freesurfer|log10', '', Phenotype)]
  # 
  # pdf(paste0('hsq-', partition, '-', gsub('\\.', '', group), '.pdf'), 11, 5)
  # gg1 <- plot_hsq_partition(hsqsummary = res[grep('left|right', Phenotype, invert = T),],
  #                                  normByPropSNPs = FALSE,
  #                                  width = 0.7)
  # gg2 <- plot_hsq_partition(hsqsummary = res[grep('left|right', Phenotype, invert = T),],
  #                                  normByPropSNPs = TRUE,
  #                                  width = 0.7)
  # gg3 <- plot_hsq_partition(hsqsummary = res[grep('left|right', Phenotype, invert = F),],
  #                                  normByPropSNPs = FALSE,
  #                                  width = 0.7)
  # gg4 <- plot_hsq_partition(hsqsummary = res[grep('left|right', Phenotype, invert = F),],
  #                                  normByPropSNPs = TRUE,
  #                                  width = 0.7)
  # dev.off()
  # 
  # pdf(paste0('hsq-', partition, '-', gsub('\\.', '', group), '_fig.pdf'), 11, 8)
  #   print(plot_grid(gg1[[2]] + xlab('') +  theme(axis.text.x =element_blank()) ,
  #             gg2[[2]],
  #             labels = c("A", "B"), align = "v",  ncol = 1, rel_heights = c(1,1.3)))
  # dev.off()
  
  # names(allres) <- c("updown-margin20-50", #"updown-margin20-40", 
  #                    'cnsexpression', 'neurodev', 'maf', 
  #                    "genic-exon1intron1",  
  #                    'genic-margin10', 'genic-margin20', 'genic-margin30', 'genic-margin40', 'genic-margin50',
  #                    "updown-margin20", "updown-margin30", "updown-margin40", "updown-margin50")
  
}  
