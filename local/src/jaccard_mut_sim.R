library(ggplot2)
library(proxy)
library(pheatmap)

pheat_f <- snakemake@output[['pheat']]
density_f <- snakemake@output[['density']]
violin_f <- snakemake@output[['violin']]
violin2_f <- snakemake@output[['violin2']]
jac_f <- snakemake@output[['jac_f']]
mw_f <- snakemake@output[['mw']]
thr <- as.numeric(snakemake@wildcards[['AF']])

load(snakemake@input[['Rimage']])

# we need to fill with 0 missing muts
m1 <- merge(pdx, pdo, all.x = TRUE, all.y = TRUE, by='row.names')
rownames(m1) <- m1$Row.names
m1$Row.names <- NULL
m1[is.na(m1)] <- 0
m1[m1 < thr] <- 0

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

ggplot(data=pdata, aes(x=pearson, color=type))+geom_density()+geom_density(size=1.5)+scale_color_manual(values=c('#004D40','#FFC107'))+theme_bw()+xlab('jaccard')
ggsave(density_f)


ggplot(data=pdata[pdata$type=="matched",], aes(y=pearson,x=""))+geom_violin()+geom_jitter(height = 0, width = 0.1)+ylim(-0.001,1.0001)+theme_bw()+ylab('Mutational Similarity?')
ggsave(violin_f)

ggplot(data=pdata, aes(y=pearson,x=type, fill=type))+geom_violin()+geom_jitter(height = 0, width = 0.1)+ylim(-0.001,1.0001)+scale_fill_manual(values=c('#004D40','#FFC107'))+theme_bw()+ylab('jaccard')
ggsave(violin2_f)

sink(mw_f)
wilcox.test(formula=as.formula("pearson~type"), data=pdata)$p.value
sink()