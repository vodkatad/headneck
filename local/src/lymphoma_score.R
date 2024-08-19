#!/usr/bin/env Rscript

library(getopt)
library(tidyverse)

set.seed(42)

opts <- matrix(c(
  'help', 'h', 0, 'logical',
  'expr_in', 'e', 1, 'character',
  'genes_in', 'g', 1, 'character',
  'output', 'o', 1, 'character'), ncol=4, byrow=TRUE)
opt <- getopt(opts)


expr <- read.table(gzfile(opt$expr_in), sep='\t', quote="", header=TRUE, row.names=1)
colnames(expr) <- gsub('.', '-', colnames(expr), fixed = TRUE) # solito problema del . nei replicati che in altri file è un -
# rownames(expr) <- str_remove(rownames(expr),"H_") # rimuovo H_ dai geni

df <- read.table(gzfile(opt$genes_in), sep='\t', quote="", header=TRUE, row.names=1)

scores <- data.frame(row.names=names(expr),stringsAsFactors = FALSE,
                     Endothelial=rep("",ncol(expr)),
                     CAF=rep("",ncol(expr)),
                     Leucocyte=rep("",ncol(expr)))
                     
levels <- levels(df$Signature)

for ( i in 1:length(levels) ) {
    for ( j in 1:length(names(expr)) ) {
      tmp_mark <- df[df[,1]==levels[i],, drop=FALSE]
      tmp_expr <- expr[rownames(expr) %in% rownames(tmp_mark),]
      scores[j,i] <- mean(tmp_expr[,j], na.rm=TRUE)
  }
}

###TODO : aggiungi la lista dei geni trovati per ogni signature

write.table(scores, gzfile(opt$output), sep='\t', quote=FALSE, col.names=TRUE, row.names=TRUE)
