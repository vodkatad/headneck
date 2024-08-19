### GSEA enrichment analysis

library(clusterProfiler)
library(tidyverse)
library(dplyr)
library(msigdbr)
library(enrichplot)
library(DOSE)
library(ggplot2)

gene_res_f <- snakemake@input[["gene_res_freq"]]
GSEA_r <- snakemake@output[["GSEA_r"]]
GSEA_ridgeplot <- snakemake@output[["GSEA_ridgeplot"]]
type <- snakemake@wildcards[["msign"]]
subclass <- snakemake@wildcards[["sub"]]
gene_res_df <- read.table(gene_res_f, quote = "", sep = "\t", header = TRUE)
###order
geneList <- gene_res_df$Freq
names(geneList) <- as.character(gene_res_df$gene)
geneList <- sort(geneList, decreasing = TRUE)

m_t2g <- msigdbr(species = "Homo sapiens", category = type)  %>% filter (gs_subcat != "CGP") %>% 
  dplyr::select(gs_name, human_gene_symbol) ### altrimenti chiede gli id numerici


#m_t2g <- m_t2g[grep(subclass, m_t2g$gs_name),]
em <- GSEA(geneList, TERM2GENE = m_t2g, pvalueCutoff = 1)

save.image("gsea_results_CP.Rdata")

#GSEA_r <- write.table(em@result, quote = FALSE, row.names = TRUE, col.names = TRUE)

write.table(em@result, file = GSEA_r, quote = FALSE, sep = "\t", row.names = TRUE,
            col.names = TRUE)

ridgeplot(em, showCategory = 20)
ggsave(GSEA_ridgeplot, width = 300, height = 107, useDingbats=FALSE, units = "mm")

#save.image('GSEA.Rdata')
