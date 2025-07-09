ba_ims <- "/mnt/cold1/snaketree/prj/hn/dataset/V1/differenziali_clustering/classification_cutoff0.05-BA.vs.IMS.deseq2.tsv"
ba_ims <- read.table(ba_ims, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
ba_ims <- ba_ims %>% filter(padj < 0.05)
ba_ims <- ba_ims %>% filter(log2FoldChange > 0.5849625)
ba_ims$genes <- rownames(ba_ims)
ba_up <- ba_ims

ba_cl <- "/mnt/cold1/snaketree/prj/hn/dataset/V1/differenziali_clustering/classification_cutoff0.05-BA.vs.CL.deseq2.tsv"
ba_cl <- read.table(ba_cl, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
ba_cl <- ba_cl %>% filter(padj < 0.05)
ba_cl <- ba_cl %>% filter(log2FoldChange > 0.5849625)
ba_cl$genes <- rownames(ba_cl)
ba_up <- rbind(ba_up,ba_cl)

cl_ims <- "/mnt/cold1/snaketree/prj/hn/dataset/V1/differenziali_clustering/classification_cutoff0.05-CL.vs.IMS.deseq2.tsv"
cl_ims <- read.table(cl_ims, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
cl_ims <- cl_ims %>% filter(padj < 0.05)
cl_ims <- cl_ims %>% filter(log2FoldChange > 0.5849625)
cl_ims$genes <- rownames(cl_ims)
cl_up <- cl_ims

ims_ba <- "/mnt/cold1/snaketree/prj/hn/dataset/V1/differenziali_clustering/classification_cutoff0.05-BA.vs.IMS.deseq2.tsv"
ims_ba <- read.table(ims_ba, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
ims_ba <- ims_ba %>% filter(padj < 0.05)
ims_ba <- ims_ba %>% filter(log2FoldChange < - 0.5849625)
ims_ba$log2FoldChange <- ims_ba$log2FoldChange * -1
ims_ba$genes <- rownames(ims_ba)
ims_up <- ims_ba

ims_cl <- "/mnt/cold1/snaketree/prj/hn/dataset/V1/differenziali_clustering/classification_cutoff0.05-CL.vs.IMS.deseq2.tsv"
ims_cl <- read.table(ims_cl, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
ims_cl <- ims_cl %>% filter(padj < 0.05)
ims_cl <- ims_cl %>% filter(log2FoldChange < - 0.5849625)
ims_cl$genes <- rownames(ims_cl)
ims_cl$log2FoldChange <- ims_cl$log2FoldChange * -1
ims_up <- rbind(ims_up, ims_cl)

cl_ba <- "/mnt/cold1/snaketree/prj/hn/dataset/V1/differenziali_clustering/classification_cutoff0.05-BA.vs.CL.deseq2.tsv"
cl_ba <- read.table(cl_ba, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
cl_ba <- cl_ba %>% filter(padj < 0.05)
cl_ba <- cl_ba %>% filter(log2FoldChange < - 0.5849625)
cl_ba$genes <- rownames(cl_ba)
cl_ba$log2FoldChange <- cl_ba$log2FoldChange * -1
cl_up <- rbind(cl_up, cl_ba)

ba_h <-"/mnt/cold1/snaketree/prj/hn/dataset/V1/differenziali_clustering/GSEA_H_results_classification_cutoff0.05-BA.vs.CL.tsv"
ba_h <- read.table(ba_h, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

ba_h <-"/mnt/cold1/snaketree/prj/hn/dataset/V1/differenziali_clustering/GSEA_H_results_classification_cutoff0.05-BA.vs.IMS.tsv"
ba_h <- read.table(ba_h, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

cl <- "/mnt/cold1/snaketree/prj/hn/dataset/V1/differenziali_clustering/GSEA_H_results_classification_cutoff0.05-CL.vs.IMS.tsv"
cl <- read.table(cl, quote = "", sep = "\t",header = TRUE, stringsAsFactors = FALSE)
