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
cor_dir <- file.path(annex_dir, "data", "raw", "CoRR")
volume.file <- file.path(cor_dir, "volumes-corr.tsv")
volume.table <- read.csv(volume.file, sep="\t")

names(volume.table)[1] <- "Folder"
volume.table$Subject <- as.numeric(sub("^sub-([^_]+)_.*$", "\\1", volume.table$Folder))
volume.table$Session <- as.numeric(sub("^.*_ses-([^_]+)_.*$", "\\1", volume.table$Folder))
volume.table$Acquisition <- NA
acquisition_lines <- grep("^.*_acq-([^_]+)_.*$", volume.table$Folder)
volume.table$Acquisition[acquisition_lines] <- sub("^.*_acq-([^_]+)_.*", "\\1$", volume.table$Folder[acquisition_lines])
volume.table$Run <- as.numeric(sub("^.*_run-([^_]+)$", "\\1", volume.table$Folder))

phenotypic.file <- file.path(cor_dir, "CoRR_AggregatedPhenotypicData.csv")
phenotypic.table <- read.csv(phenotypic.file, na.strings="#") # , colClasses=c(SUBID="character"))

phenotypic.table$Subject <- phenotypic.table$SUBID
phenotypic.table$Session <- NA
phenotypic.table$Session[grep("Baseline", phenotypic.table$SESSION)] <- 1
retest_lines <- grep("Retest", phenotypic.table$SESSION)
phenotypic.table$Session[retest_lines] <- as.numeric(sub("Retest_([0-9]+)", "\\1", phenotypic.table$SESSION[retest_lines])) + 1

phenotypic.table$Age <- as.numeric(phenotypic.table$AGE_AT_SCAN_1)
phenotypic.table$Age[is.na(phenotypic.table$Age)] <- phenotypic.table$Age[match(phenotypic.table$SUBID[is.na(phenotypic.table$Age)], phenotypic.table$SUBID)]
which(is.na(phenotypic.table$Age))

phenotypic.table$Sex <- factor(phenotypic.table$SEX, labels = c("F", "M"))
phenotypic.table$Sex[is.na(phenotypic.table$Sex)] <- phenotypic.table$Sex[match(phenotypic.table$SUBID[is.na(phenotypic.table$Sex)], phenotypic.table$SUBID)]
which(is.na(phenotypic.table$Sex))

# dsubj = c()
# for (subj in unique(phenotypic.table$SUBID)) {
#   ages <- phenotypic.table$AGE_AT_SCAN_1[phenotypic.table$SUBID == subj]
#   if (!all(ages == ages[1], na.rm=T)) {
#     dsubj = c(dsubj, subj)
#   }
# }
# dphe = subset(phenotypic.table, SUBID %in% dsubj)
# View(dphe)

outlier.col <- c("EstimatedTotalIntraCranialVol", "BrainSeg")
outlier.table <- volume.table[, outlier.col]

complete.rows <- complete.cases(outlier.table)
volume.table <- volume.table[complete.rows,]
outlier.table <- outlier.table[complete.rows,]

pca <- prcomp(outlier.table)

median <- apply(pca$x, 2, median)
mad <- apply(pca$x, 2, function(x) median(abs(median(x)-x)))

density <- apply(dnorm(t(t(pca$x)/mad*qnorm(3/4))), 1, prod)
density0 <- dnorm(0) ** ncol(pca$x)
outliers <- density < density0/100

plot(outlier.table, asp=1, col=ifelse(outliers, "red", "blue"))
b1 <- (pca$rotation %*% c(1, 0))[2] / (pca$rotation %*% c(1, 0))[1]
a1 <- pca$center[2] - b1 * pca$center[1]
abline(a1, b1, col="black")
b2 <- (pca$rotation %*% c(0, 1))[2] / (pca$rotation %*% c(0, 1))[1]
a2 <- pca$center[2] - b2 * pca$center[1]
abline(a2, b2, col="black")

table(outliers)

# remove outliers from volume table
volume.table <- volume.table[!(outliers),-1]
volume.table.log10 <- log10(Filter(is.numeric, volume.table))
volume.table.used <- volume.table

merged.table <- merge(phenotypic.table, volume.table.used, by=c("Subject", "Session"))

# Remove NYU_1 volumes as sessions are all with the same anatomical MRI
merged.table <- subset(merged.table, SITE != "NYU_1")

# Remove MPG_1 volumes as they have only one session with different acquisitions parameters
merged.table <- subset(merged.table, SITE != "MPG_1")

# Remove SWU_3 volumes as they seem to use same anatomical MRI for each session
merged.table <- subset(merged.table, SITE != "SWU_3")


fs.volumes <- grep("hypointensities", head(names(volume.table.used), ncol(volume.table.used)-4), value=T, inv=T)
left.volumes <- grep("^(Left\\.|lh)", fs.volumes, value=T)
right.volumes <- grep("^(Right\\.|rh)", fs.volumes, value=T)

diff.vol <- merged.table[, left.volumes] - merged.table[, right.volumes]
names(diff.vol) <- sub("^Left\\.|lh", "Diff.", names(diff.vol))

sum.vol <- merged.table[, left.volumes] + merged.table[, right.volumes]
names(sum.vol) <- sub("^Left\\.|lh", "Sum.", names(sum.vol))

merged.table <- cbind(merged.table, diff.vol, sum.vol)

first.table <- merged.table[!duplicated(merged.table$Subject),]

vol.names <- c(fs.volumes, names(diff.vol), names(sum.vol))

for (vol in vol.names) {
    form <- as.formula(paste(vol, "~Age+Sex+SITE"))
    mod <- lm(form, data=first.table)
    merged.table[[paste0(vol, ".resid")]] <- merged.table[[vol]] - predict(mod, merged.table)
}
  
ses1.table <- subset(merged.table, Session==1)

if (any(duplicated(ses1.table[, c("Subject", "Run")]))) {
  print(paste("Warning: duplicates with same run found"))
}
run.tables <- split(ses1.table, ses1.table$Run)
table.runs <- sapply(vol.names, function(vol) {
  vol.runs <- lapply(run.tables, '[', c("Subject", "Session", paste0(vol, ".resid")))
  # ignore duplicated name warnings
  table <- suppressWarnings(Reduce(function(x,y) merge(x, y, all=T, by=c("Subject", "Session")), vol.runs))
  # as.matrix(table[,-c(1,2)])
  table
}, simplify=F)

run1.table <- subset(merged.table, Run==1)

if (any(duplicated(run1.table[, c("Subject", "Session")]))) {
  print(paste("Warning: duplicates with same session found"))
}
ses.tables <- split(run1.table, run1.table$Session)
table.sessions <- sapply(vol.names, function(vol) {
  vol.sessions <- lapply(ses.tables, '[', c("Subject", "Run", paste0(vol, ".resid")))
  # ignore duplicated name warnings
  table <- suppressWarnings(Reduce(function(x,y) merge(x, y, all=T, by=c("Subject", "Run")), vol.sessions))
  # as.matrix(table[,-c(1,2)])
  table
}, simplify=F)

icc.runs <- sapply(vol.names, function(vol) {
  table <- na.omit(table.runs[[vol]][, c(3,4)])
  irr::icc(table)
}, simplify=F)

icc.sessions <- sapply(vol.names, function(vol) {
  table <- na.omit(table.sessions[[vol]][, c(3,4)])
  irr::icc(table)
}, simplify = F)

cor.left.right <- diag(cor(ses.tables[[1]][left.volumes], ses.tables[[1]][right.volumes]))
names(cor.left.right) <- sub("^(Left\\.|lh)", "", left.volumes)

cor.sessions <- sapply(vol.names, function(vol1) {
  sapply(vol.names, function(vol2) {
    table1 <- as.matrix(na.omit(table.sessions[[vol1]][, c(3,4)]))
    table1 <- table1 - mean(table1)
    table2 <- as.matrix(na.omit(table.sessions[[vol2]][, c(3,4)]))
    table2 <- table2 - mean(table2)
    cov.table <- t(table1) %*% table2
    mean(cov.table[lower.tri(cov.table)|upper.tri(cov.table)]) / mean(diag(cov.table))
  })
})

icc.sessions.values <- sapply(icc.sessions, function(x) if(is.na(x[1])) NA else x$value)
icc.sessions.n <- sapply(icc.sessions, function(x) if(is.na(x[1])) NA else x$subjects)
icc.sessions.k <- sapply(icc.sessions, function(x) if(is.na(x[1])) NA else x$raters)
icc.sessions.lb <- sapply(icc.sessions, function(x) if(is.na(x[1])) NA else x$lbound)
icc.sessions.ub <- sapply(icc.sessions, function(x) if(is.na(x[1])) NA else x$ubound)
icc.sessions.se <- (icc.sessions.ub - icc.sessions.lb) / (qnorm(0.975) * 2)

cor.table = data.frame(ICC = icc.sessions.values, ICC.lb = icc.sessions.lb, ICC.ub = icc.sessions.ub)
cor.file <- file.path(script.dir, "cor-fs-volumes.tsv")
write.table(cor.table, cor.file, col.names = NA, sep = "\t")

# cor.left.right.table <- sapply(res.sites, "[[", "cor.left.right")
# cor.sessions.table <- sapply(sapply(res.sites, "[[", "cor.sessions", simplify=F), "[", cbind(left.volumes, right.volumes))
# rownames(cor.sessions.table) <- sub("^(Left\\.|lh)", "", left.volumes)
# icc.sessions.left <- icc.sessions.values[left.volumes,]
# icc.sessions.right <- icc.sessions.values[right.volumes,]

# cor.left.right.table.cor1 <- cor.left.right.table / sqrt(icc.sessions.left * icc.sessions.right)
# cor.left.right.table.cor2 <- cor.left.right.table / sqrt(icc.sessions.left * icc.sessions.right) * cor.sessions.table

# sel <- icc.sessions.n[1, ] * icc.sessions.k[1, ] >= 50 & !is.na(icc.sessions.k[1,])
# cor.left.right.table.cor1.sel <- cor.left.right.table.cor1[,sel]
# cor.left.right.table.cor2.sel <- cor.left.right.table.cor2[,sel]


cor.sessions.table <- sapply(vol.names, function(vol1) {
    t(sapply(vol.names, function(vol2) {
      c(vol1, vol2, cor.sessions[vol1, vol2])
    }, simplify = T))
  }, simplify=F)
cor.sessions.table <- do.call(rbind, cor.sessions.table)
cor.sessions.file <- file.path(script.dir, "cor-cov-fs-volumes.tsv")
write.table(cor.sessions.table, cor.sessions.file, col.names=c("Volume1", "Volume2", "ICC_cov"), row.names=FALSE, sep='\t')

sites.sessions <- unique(subset(merged.table, Subject %in% table.sessions[[vol[1]]][cc, "Subject"])$SITE)
length(sites.sessions)
