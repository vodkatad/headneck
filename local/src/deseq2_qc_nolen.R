#!/usr/bin/env Rscript
#options(error=traceback

args <- commandArgs(trailingOnly=TRUE)
counts <- args[[1]]
metadataf <- args[[2]]
design <- args[[3]]
prefix <- args[[4]] # if all do not filter, otherwise keep ^prefix_ only genes. 
#useful for xeno mice-man counts matrices made by Ivan's pipeline.
minc <- as.numeric(args[[5]])
minsamples <- as.numeric(args[[6]])
image <- args[[7]]
tmmf <- args[[8]]
vsdf <- args[[9]]
cores <- as.numeric(args[[10]])

outprefix <- 'qc'

save.image(image)
library("BiocParallel")
register(MulticoreParam(cores))
library(DESeq2)
library(RColorBrewer)
library(pheatmap)
library(ggplot2)
library(edgeR)

data <- read.table(gzfile(counts), header=T, sep="\t", row=1)
if (prefix != "all") {
  data <- data[grep(paste0("^", prefix, "_"), rownames(data)),]
}
metadata <- read.table(metadataf, sep="\t", header=T, row=1)

if ('batch' %in% colnames(metadata)) {
    metadata[,'batch'] <- as.factor(metadata[,'batch'])
}

fdesign <- as.formula(design)
print(terms(fdesign)[[2]])

# this was needed for ad hoc biodiversa stuff/strunz
#rownames(metadata) <- gsub("-", ".", rownames(metadata), fixed=TRUE)
###  colnames(data) <- gsub(".2", "", colnames(data), fixed=TRUE)
if (length(intersect(rownames(metadata), colnames(data))) != nrow(metadata)) {
    stop('No correspondence between metadata and counts!')
}

new_data <- data[,match(rownames(metadata), colnames(data))]
if (!all(rownames(metadata)==colnames(new_data))) {
    stop('match issues...')
}
## TMM from edgeR without filtering expression ###
#group <- as.factor(metadataf$terms(fdesign)[[2]])
#manca filtro idiota!
# TMM Should be done after selectoin of samples though, move!
efilterGenes <- rowSums(new_data > minc) < minsamples
edata <- new_data[!efilterGenes,]
e_new_data <- edata[,match(rownames(metadata), colnames(edata))]

###esegui fino a qui e da qui faccio limma su un confronto
y <- DGEList(counts=e_new_data, samples=metadata)
#https://www.biostars.org/p/317701/
y <- calcNormFactors(y)
tmm <- cpm(y)
write.table(tmm, gzfile(tmmf), quote=F, row.names=T, col.names=T, sep="\t")
####

dds <- DESeqDataSetFromMatrix(countData = new_data, colData = metadata, design = fdesign)
# filterGenes are the genes that will be removed cause they have 'noise reads' in less than minsamples
filterGenes <- rowSums(counts(dds) > minc) < minsamples
dds <- dds[!filterGenes]
dds <- DESeq(dds, parallel=TRUE, betaPrior=TRUE)
vsd <- vst(dds, blind=FALSE)

save.image(image)

write.table(assay(vsd), gzfile(vsdf), quote=F, row.names=T, col.names=T, sep="\t")


sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- colnames(vsd)
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatf <- paste0(outprefix, "_vsd_dists.pdf")
pheatmap(sampleDistMatrix,clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, col=colors, file=pheatf)
if (design != "~1") {
    out <- tryCatch(
        {
           ints <- attr(terms(fdesign),"term.labels")
           for (i in seq(1, length(ints))) {
               plotPCA(vsd, intgroup=ints[i])
               ggsave(paste0(outprefix, "_", ints[i], "_pca.pdf"))
           }
        },
        error=function(cond) {
            message(cond)
        },
        warning=function(cond) {
            message(cond)
        },
        finally=function(cond) { dev.off() }
    )
}

cc <- counts(dds,normalized=TRUE)
select <- order(rowMeans(counts(dds,normalized=TRUE)),
                decreasing=TRUE)[1:20]
if (design != "~1") {
    df <- as.data.frame(colData(dds)[,attr(terms(fdesign),"term.labels"), drop=F])
} else {
    df <- as.data.frame(colData(dds)[,c(terms(fdesign)[[2]]),drop=F])
}
pheatf2 <- paste0(outprefix, "_cc_high.pdf")
pheatmap(cc[select,], cluster_rows=FALSE, show_rownames=FALSE, cluster_cols=FALSE, annotation_col=df, file=pheatf2)
highgenes <- rownames(cc[select,])
sds <- apply(cc, 1, sd)
highsd <- cc[order(sds, decreasing=TRUE)[1:20],]
pheatf3 <- paste0(outprefix, "_cc_highsd.pdf")
pheatmap(highsd, cluster_rows=FALSE, show_rownames=FALSE, cluster_cols=FALSE, annotation_col=df, file=pheatf3)
highsdgenes <- rownames(highsd)

res <- data.frame(highg=highgenes, highsdgenes=highsdgenes)
write.table(res, file=paste0(outprefix, "_high.tsv"), sep="\t", quote=FALSE, row.names=FALSE)

save.image(image)
