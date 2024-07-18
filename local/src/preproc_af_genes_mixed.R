library(ggplot2)

xeno_af <- snakemake@input[['matX']]
human_af <- snakemake@input[['matH']]
log_f <- snakemake@log[['log']]

savedata <- snakemake@output[['preprocAF']]

xthr <- as.numeric(snakemake@wildcards[['XAF']])
hthr <- as.numeric(snakemake@wildcards[['HAF']])

xeno <- read.table(xeno_af, sep="\t", header=T, row.names=1)
human <- read.table(human_af, sep="\t", header=T, row.names=1)

# o takes the h place
pdobing <- ifelse(human > hthr, TRUE, FALSE)
pdxbing <- ifelse(xeno > xthr, TRUE, FALSE)

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

colnames(pdobing) <- substr(colnames(pdobing), 0, 7)
colnames(pdxbing) <- substr(colnames(pdxbing), 0, 7)

pdobing <- pdobing[, common]
pdxbing <- pdxbing[, common]

save.image(savedata)