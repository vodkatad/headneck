library(DESeq2)
library(edgeR)
library(RColorBrewer)
library(pheatmap)
library(ggplot2)
setwd('/mnt/cold1/snaketree/prj/hn/local/share/data')
counts <- read.table('Merged_raw_counts_xeno.txt', sep="\t", header=TRUE, row.names=1)
annot <- read.table('DNA_RNA_HNC_filtersFPlympho_qc_defResponses_shared.csv', sep="\t", header=TRUE)

annot1 <- annot[, c('Case','Average.3.6WKS', 'Definitive.resp')]
df_annot <- data.frame(sample=colnames(counts), smodel= substr(colnames(counts),0, 7))

mannot <- merge(df_annot, annot1, by.x="smodel", by.y= 'Case')
mannot <- mannot[mannot$Definitive.resp != "",]
mannot$Average.3.6WKS <- gsub('%','', mannot$Average.3.6WKS, fixed=TRUE)
mannot$Definitive.resp <- gsub(' ',"_", mannot$Definitive.resp, fixed=TRUE)
mannot$Average.3.6WKS <- as.numeric(mannot$Average.3.6WKS)
rownames(mannot) <- mannot$sample
mannot$sample <- NULL
metadata <- mannot

design <- '~ Definitive.resp'

counts_fil <- counts[, colnames(counts) %in% rownames(metadata)]

minc <- 5
minsamples <- 3

new_data <- counts_fil[,match(rownames(metadata), colnames(counts_fil))]
if (!all(rownames(metadata)==colnames(new_data))) {
  stop('match issues...')
}


efilterGenes <- rowSums(new_data > minc) < minsamples
edata <- new_data[!efilterGenes,]
e_new_data <- edata[,match(rownames(metadata), colnames(edata))]

###esegui fino a qui e da qui faccio limma su un confronto
y <- DGEList(counts=e_new_data, samples=metadata)
#https://www.biostars.org/p/317701/
y <- calcNormFactors(y)
tmm <- cpm(y)
####

dds <- DESeqDataSetFromMatrix(countData = new_data, colData = metadata, design = as.formula(design))
filterGenes <- rowSums(counts(dds) > minc) < minsamples
dds <- dds[!filterGenes]
dds <- DESeq(dds, parallel=TRUE, betaPrior=TRUE)
vsd <- vst(dds, blind=FALSE)


sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- colnames(vsd)
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,clustering_distance_rows=sampleDists, clustering_distance_cols=sampleDists, col=colors)
plotPCA(vsd, intgroup='Definitive.resp')


dim(new_data)
dim(metadata)
table(filterGenes)
##
metadata <- all_metadata
tmm <- all_tmm
vsd <- all_vsd

genes <- c('KLF5', 'SERPINB8', 'CD24', 'ATP10B','FAXDC2', 'SASH1', 'KRT13', 'SP1')

if (!all(colnames(tmm)==rownames(metadata))) {
  stop('Piciucina')
}

if (!all(colnames(vsd)==rownames(metadata))) {
  stop('Piciucina')
}

cor.test(tmm[rownames(tmm)=="KLF5",],tmm[rownames(tmm)=="SP1",] )
avsd <- assay(vsd)
cor.test(avsd[rownames(avsd)=="KLF5",],avsd[rownames(avsd)=="SP1",] )

compare_two <- function(gene, expr, metadata, dolog=TRUE) {
  if (gene %in% rownames(expr)) {
    data <- data.frame(expr=expr[rownames(expr)==gene,], cetuxi=metadata$Average.3.6WKS)
    if (dolog){
      data$expr <- log(data$expr+1)
    }
    gg1 <- ggplot(data=data, aes(x=expr, y=cetuxi))+geom_point()+geom_smooth(method="lm")+theme_bw()+xlab(gene)
    ci <- cor.test(data$expr, data$cetuxi)
    return(list(pearson=c(ci$estimate, ci$p.value),plot=gg1))
  }
    return(list(pearson=c(NA, NA), plot=NULL))
}


tmm_cors <- lapply(genes, compare_two, tmm, metadata, TRUE)

pearson <- t(sapply(tmm_cors, function(x){x$pearson}))
rownames(pearson) <- genes
colnames(pearson) <- c('pearson', 'pvalue')
pearson
tmm_cors[[2]][['plot']]

tmm_cors <- lapply(genes, compare_two, assay(vsd), metadata, TRUE)
pearson <- t(sapply(tmm_cors, function(x){x$pearson}))
rownames(pearson) <- genes
colnames(pearson) <- c('pearson', 'pvalue')
pearson
tmm_cors[[1]][['plot']]

# 1 rappresentante per modello:
#egrassi@ulisse:/mnt/cold1/snaketree/prj/hn/dataset/V1/WES$ head -n1 xeno_gene_matrix.tsv  | tr "\t" "\n"| sed 1d | bawk '{print substr($1,0,7)}' | sort | uniq -d

# repeat on TP53 mutated only
muts <- read.table('/mnt/cold1/snaketree/prj/hn/dataset/V1/WES/xeno_long_fx.tsv', sep="\t", header=T)
tp53 <- muts[muts$Hugo_Symbol=="TP53",]

mutated <- substr(tp53$Tumor_Sample_Barcode, 0,7)

all_metadata <- metadata
all_tmm <- tmm
all_vsd <- vsd

metadata <- metadata[substr(rownames(metadata),0,7) %in% mutated,]
tmm <- tmm[,substr(colnames(tmm),0,7) %in% mutated]
vsd <- vsd[,substr(colnames(vsd),0,7) %in% mutated]
dim(tmm)
dim(metadata)
dim(vsd)

if (!all(colnames(tmm)==rownames(metadata))) {
  stop('Piciucina')
}

if (!all(colnames(vsd)==rownames(metadata))) {
  stop('Piciucina')
}

cor.test(tmm[rownames(tmm)=="KLF5",],tmm[rownames(tmm)=="SP1",] )
avsd <- assay(vsd)
cor.test(avsd[rownames(avsd)=="KLF5",],avsd[rownames(avsd)=="SP1",] )

data <- data.frame(KLF5=tmm[rownames(tmm)=="KLF5",],SP1=tmm[rownames(tmm)=="SP1",])
data$KLF5 <- log(data$KLF5+1)
data$SP1 <- log(data$SP1+1)
ggplot(data=data, aes(x=KLF5, y=SP1))+geom_point()+geom_smooth(method="lm")+theme_bw()
cor.test(data$KLF5, data$SP1)

tmm_cors <- lapply(genes, compare_two, tmm, metadata, TRUE)
pearson <- t(sapply(tmm_cors, function(x){x$pearson}))
rownames(pearson) <- genes
colnames(pearson) <- c('pearson', 'pvalue')
pearson
tmm_cors[[2]][['plot']]

tmm_cors <- lapply(genes, compare_two, assay(vsd), metadata, TRUE)
pearson <- t(sapply(tmm_cors, function(x){x$pearson}))
rownames(pearson) <- genes
colnames(pearson) <- c('pearson', 'pvalue')
pearson
tmm_cors[[1]][['plot']]

# repeat on TP53 wt only
metadata <- all_metadata
tmm <- all_tmm
vsd <- all_vsd

metadata <- metadata[!substr(rownames(metadata),0,7) %in% mutated,]
tmm <- tmm[,!substr(colnames(tmm),0,7) %in% mutated]
vsd <- vsd[,!substr(colnames(vsd),0,7) %in% mutated]
dim(tmm)
dim(metadata)
dim(vsd)


if (!all(colnames(tmm)==rownames(metadata))) {
  stop('Piciucina')
}

if (!all(colnames(vsd)==rownames(metadata))) {
  stop('Piciucina')
}

cor.test(tmm[rownames(tmm)=="KLF5",],tmm[rownames(tmm)=="SP1",] )
avsd <- assay(vsd)
cor.test(avsd[rownames(avsd)=="KLF5",],avsd[rownames(avsd)=="SP1",] )

data <- data.frame(KLF5=tmm[rownames(tmm)=="KLF5",],SP1=tmm[rownames(tmm)=="SP1",])
data$KLF5 <- log(data$KLF5+1)
data$SP1 <- log(data$SP1+1)
ggplot(data=data, aes(x=KLF5, y=SP1))+geom_point()+geom_smooth(method="lm")+theme_bw()
cor.test(data$KLF5, data$SP1)


tmm_cors <- lapply(genes, compare_two, tmm, metadata, TRUE)
pearson <- t(sapply(tmm_cors, function(x){x$pearson}))
rownames(pearson) <- genes
colnames(pearson) <- c('pearson', 'pvalue')
pearson
tmm_cors[[1]][['plot']]
tmm_cors[[2]][['plot']]
tmm_cors[[4]][['plot']]
tmm_cors[[8]][['plot']]


tmm_cors <- lapply(genes, compare_two, assay(vsd), metadata, TRUE)
pearson <- t(sapply(tmm_cors, function(x){x$pearson}))
rownames(pearson) <- genes
colnames(pearson) <- c('pearson', 'pvalue')
pearson
tmm_cors[[1]][['plot']]


# check on humans:
#egrassi@ulisse:/mnt/cold1/snaketree/prj/hn/dataset/V1/WES$ grep TP53 human_long_fx.tsv  | grep HNC0089
#chr17_7675997_G:T	TP53	HNC0089PRH00	p.C124*	stop_gained	0.695652173913043
#egrassi@ulisse:/mnt/cold1/snaketree/prj/hn/dataset/V1/WES$ grep TP53 human_long_fx.tsv  | grep HNC0098
#egrassi@ulisse:/mnt/cold1/snaketree/prj/hn/dataset/V1/WES$ grep TP53 human_long_fx.tsv  | grep HNC0107
#egrassi@ulisse:/mnt/cold1/snaketree/prj/hn/dataset/V1/WES$ grep TP53 human_long_fx.tsv  | grep HNC0148
#egrassi@ulisse:/mnt/cold1/snaketree/prj/hn/dataset/V1/WES$ grep TP53 human_long_fx.tsv  | grep HNC0151
#egrassi@ulisse:/mnt/cold1/snaketree/prj/hn/dataset/V1/WES$ grep TP53 human_long_fx.tsv  | grep HNC0164
