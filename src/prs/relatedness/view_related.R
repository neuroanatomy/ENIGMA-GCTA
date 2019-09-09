# R script to read the GRM binary file
ReadGRMBin=function(prefix, AllN=F, size=4){
  sum_i=function(i){
    return(sum(1:i))
  }
  BinFileName=paste(prefix,".grm.bin",sep="")
  NFileName=paste(prefix,".grm.N.bin",sep="")
  IDFileName=paste(prefix,".grm.id",sep="")
  id = read.table(IDFileName)
  n=dim(id)[1]
  BinFile=file(BinFileName, "rb");
  grm=readBin(BinFile, n=n*(n+1)/2, what=numeric(0), size=size)
  NFile=file(NFileName, "rb");
  if(AllN==T){
    N=readBin(NFile, n=n*(n+1)/2, what=numeric(0), size=size)
  }
  else N=readBin(NFile, n=1, what=numeric(0), size=size)
  i=sapply(1:n, sum_i)
  return(list(diag=grm[i], off=grm[-i], id=id, N=N))
}

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
grm_file <- file.path(ukb_annex, "data/derived/test/04.grm/grm-all/all")
rel_file <- "test.dat"
threshold <- 0.0125


rel_table <- read.table(rel_file, header=T)
rel_table <- subset(rel_table, ID1>0 & ID2>0)

grm <- ReadGRMBin(grm_file)
GRM=matrix(0,nrow(grm$id),nrow(grm$id))
diag(GRM)=grm$diag
GRM[upper.tri(GRM)]=grm$off
GRM[lower.tri(GRM)]=t(GRM)[lower.tri(t(GRM))]
rownames(GRM)=grm$id[,2]
colnames(GRM)=rownames(GRM)

indexes <- which(GRM >= threshold, arr.ind = T)
indexes <- indexes[indexes[,1] < indexes[,2],]
related <- data.frame(ID1=indexes[,2],
                      ID2=indexes[,1],
                      relatedness=GRM[indexes])

related[related$ID1 > related$ID2, c("ID1", "ID2")] <- related[related$ID1 > related$ID2, c("ID2", "ID1")]
which(related$ID1 > related$ID2)

rel_table[rel_table$ID1 > rel_table$ID2, c("ID1", "ID2")] <- rel_table[rel_table$ID1 > rel_table$ID2, c("ID2", "ID1")]
which(rel_table$ID1 > rel_table$ID2)

merged_table <- merge(related, rel_table, all=T)
merged_table1 <- subset(merged_table, !is.na(relatedness))
