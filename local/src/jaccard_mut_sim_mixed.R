library(ggplot2)
library(proxy)
library(pheatmap)

pheat_f <- snakemake@output[['pheat']]
jac_f <- snakemake@output[['jac_f']]
mw_f <- snakemake@output[['mw']]
xthr <- as.numeric(snakemake@wildcards[['XAF']])
hthr <- as.numeric(snakemake@wildcards[['HAF']])

load(snakemake@input[['Rimage']])

# we need to fill with 0 missing muts
pdo[pdo < hthr] <- 0 
pdx[pdx < xthr] <- 0 
m1 <- merge(pdx, pdo, all.x = TRUE, all.y = TRUE, by='row.names')
rownames(m1) <- m1$Row.names
m1$Row.names <- NULL
m1[is.na(m1)] <- 0

mx <- m1[,grepl('.x', colnames(m1))]
mo <- m1[,grepl('.y', colnames(m1))]
colnames(mx) <- substr(colnames(mx), 0, 7)
colnames(mo) <- substr(colnames(mo), 0, 7)


jac <- proxy::simil(mx, mo, by_rows = FALSE, method = "Jaccard")

pdf(pheat_f, family="sans")# height=2.36, width=2.36)
pheatmap(jac, cluster_rows=F , cluster_cols=F, labels_col="Human", labels_row="Xeno", fontsize.number=1.5, color=colorRampPalette(c("white", "red"))(50), legend=FALSE)
graphics.off()

write.table(jac, jac_f, sep="\t", quote=FALSE)
pearson <- jac
diag <- diag(pearson)
pearson2 <- pearson
diag(pearson2) <- rep(NA, length(diag))
all <- as.numeric(unlist(pearson2))
all <- all[!is.na(all)]
#all <- upper.tri(pearson, diag = FALSE) # this is not a simmetric matrix!
pdata <- data.frame(pearson=c(all, diag), type=c(rep('unmatched', length(all)),rep('matched', length(diag))))



sink(mw_f)
wilcox.test(formula=as.formula("pearson~type"), data=pdata)$p.value
sink()