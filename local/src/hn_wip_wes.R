library(ComplexHeatmap)
library(ggplot2)
library(ggrepel)
library(ggsignif)
setwd('/mnt/cold1/snaketree/prj/hn/dataset/V1/WES')

################## oncoprint aggiuntivi
load('oncoprint_0.2_0.05.RData')

oncoPrint(mat_list, alter_fun = alter_fun, col = col, column_order = names(su),
          remove_empty_columns = TRUE, remove_empty_rows = TRUE, pct_gp=gpar(fontsize=3), column_names_gp=gpar(fontsize=3),
          top_annotation = HeatmapAnnotation(cbar = anno_oncoprint_barplot(), annotation_name_gp=gpar(fontsize=3), cn=anno_text(colnames(mat_list[[1]]))))

oncoPrint(mat_list, alter_fun = alter_fun, col = col, column_order = names(su),
          remove_empty_columns = FALSE, remove_empty_rows = TRUE, pct_gp=gpar(fontsize=3), column_names_gp=gpar(fontsize=3),
          top_annotation = HeatmapAnnotation(cbar = anno_oncoprint_barplot(), annotation_name_gp=gpar(fontsize=3), cn=anno_text(colnames(mat_list[[1]]))))

################### correlazioni freq
tcga <- read.table('TCGA_freqs.tsv', sep="\t", header=F)
human <- read.table('human_freqs_0.05.tsv', sep="\t", header=F)

m <- merge(tcga, human, by="V1")
colnames(m) <- c('gene','TCGA', 'us')

cor.test(m$TCGA, m$us, method="spearman")
cor.test(m$TCGA, m$us, method="pearson")

# limit to mutsig significative in TCGA?
stcga <- read.table('TCGA_sign.tsv', sep="\t", header=F)
m2 <- m[m$gene %in% stcga$V1,]

cor.test(m2$TCGA, m2$us, method="spearman")
cor.test(m2$TCGA, m2$us, method="pearson")
m2 <- m2[order(m2$TCGA),]
lab <- tail(m2, n=7)
ggplot(data=m2, aes(x=TCGA, y=us))+ geom_point()+ geom_smooth(method=lm, se=TRUE, size=0.3)+theme_bw()+
  xlab('Frequency TCGA') + ylab('Our Frequency')+geom_text_repel(data=lab, aes(label=gene))


ggplot(data=m, aes(x=TCGA, y=us))+ geom_point()+ geom_smooth(method=lm, se=TRUE, size=0.3)+theme_bw()+
  xlab('Frequency TCGA') + ylab('Our Frequency')+geom_text_repel(data=lab, aes(label=gene))


tcga$is_mutsigcv <- ifelse(tcga$V1 %in% stcga$V1, 'yes', 'no')
tcga$we_have <- ifelse(tcga$V1 %in% human$V1, 'yes', 'no')

msig <- tcga[tcga$is_mutsigcv=="yes",]
ggplot(data=msig, aes(y=V2, x=we_have))+geom_boxplot(outlier.shape=NA)+geom_jitter(height=0)+theme_bw()
