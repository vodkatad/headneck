library(ggplot2)

xeno_af <- snakemake@input[['matX']]
human_af <- snakemake@input[['matH']]
log_f <- snakemake@log[['log']]

savedata <- snakemake@output[['preprocAF']]

xeno <- read.table(xeno_af, sep="\t", header=T, row.names=1)
human <- read.table(human_af, sep="\t", header=T, row.names=1)

pdo <- human
pdx <- xeno

xenos <- substr(colnames(xeno), 0, 7)
human <- substr(colnames(human), 0, 7)
common <- intersect(human, xenos)

sink(log_f)
print('human')
length(human)
print('xenos')
length(xenos)
print('common')
length(common)
sink()

colnames(pdo) <- substr(colnames(pdo), 0, 7)
colnames(pdx) <- substr(colnames(pdx), 0, 7)

pdo <- pdo[, common]
pdx <- pdx[, common]

save.image(savedata)