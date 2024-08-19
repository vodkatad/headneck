library(ComplexHeatmap)
library(ggplot2)
library(ggrepel)
library(ggsignif)
setwd('/mnt/cold1/snaketree/prj/hn/dataset/V1/WES')
textSize <- 15
largerSize <- 20
slidetheme <- theme(
  text = element_text(size = textSize, family='sans'),
  axis.title = element_text(size = largerSize),
  axis.text.x = element_text(size = textSize, color="black"),#, angle = 90, vjust = 0.5, hjust=1)
  axis.text.y = element_text(size = textSize, color="black"),
  plot.title = element_text(size = largerSize, hjust = 0.5),
  legend.title = element_text(size=largerSize, hjust = 0.5),
  legend.text = element_text(size=textSize),
  panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  axis.line = element_line(colour = "black"),
  axis.ticks = element_line(color = "black"),
  panel.background = element_blank()
)


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

ggplot(data=m, aes(x=TCGA, y=us))+ geom_point()+ geom_smooth(method=lm, se=TRUE, size=0.3)+theme_bw()+
  xlab('Frequency TCGA') + ylab('Our Frequency')+geom_text_repel(data=lab, aes(label=gene))

# limit to mutsig significative in TCGA?
#q < 0.1 TODO
#https://www.nature.com/articles/nature14129
stcga <- read.table('TCGA_sign.tsv', sep="\t", header=F)
m2 <- m[m$gene %in% stcga$V1,]

cor.test(m2$TCGA, m2$us, method="spearman")
cor.test(m2$TCGA, m2$us, method="pearson")
m2 <- m2[order(m2$TCGA),]
lab <- tail(m2, n=7)

guess_ticks <- function(values, nticks=5, fixed_max=NULL, fixed_min=0) {
  vmax <- max(values)
  if (is.null(fixed_max)) {
    round_max <- ceiling(vmax)
  } else {
    round_max <- fixed_max
  }
  my_breaks <- seq(fixed_min, round_max, length.out=nticks)
  return(my_breaks)
}


x_breaks <- guess_ticks(m2$TCGA, fixed_max=0.8)
y_breaks <- guess_ticks(m2$us, fixed_max=0.8)
  
  ggplot(data=m2, aes(x=TCGA, y=us))+ geom_point(size=0.8)+ geom_smooth(method=lm, se=TRUE, size=0.3)+slidetheme+
  xlab('Frequency (TCGA, n=530)') + ylab('Frequency (HNC, n=70)')+geom_text_repel(data=lab, aes(label=gene))+
  scale_y_continuous(breaks=y_breaks,limits=c(-0.01, max(y_breaks)),expand = c(0, 0))+ #ylim(min(y_breaks),max(y_breaks))+
    scale_x_continuous(breaks=x_breaks,limits=c(-0.01, max(x_breaks)),expand = c(0, 0)) #+xlim(min(x_breaks),max(x_breaks))
  
  ggsave('~/hn.svg', height=9, width=9)

  write.table(m2, '~/hn.tsv', sep="\t", quote=F)
  
##########
ggplot(data=m, aes(x=TCGA, y=us))+ geom_point()+ geom_smooth(method=lm, se=TRUE, size=0.3)+theme_bw()+
  xlab('Frequency TCGA') + ylab('Our Frequency')+geom_text_repel(data=lab, aes(label=gene))


m4 <- m[m$gene != 'TP53',]
cor.test(m4$TCGA, m4$us, method="spearman")
cor.test(m4$TCGA, m4$us, method="pearson")

ggplot(data=m4, aes(x=TCGA, y=us))+ geom_point()+ geom_smooth(method=lm, se=TRUE, size=0.3)+theme_bw()+
  xlab('Frequency TCGA') + ylab('Our Frequency')+geom_text_repel(data=lab, aes(label=gene))

tcga$is_mutsigcv <- ifelse(tcga$V1 %in% stcga$V1, 'yes', 'no')
tcga$we_have <- ifelse(tcga$V1 %in% human$V1, 'yes', 'no')

msig <- tcga[tcga$is_mutsigcv=="yes",]
ggplot(data=msig, aes(y=V2, x=we_have))+geom_boxplot(outlier.shape=NA)+geom_jitter(height=0)+theme_bw()



###
m5 <- m4[m4$TCGA < 0.2,]
cor.test(m5$TCGA, m5$us, method="spearman")
cor.test(m5$TCGA, m5$us, method="pearson")

ggplot(data=m4, aes(x=log2(TCGA), y=log2(us)))+ geom_point()+ geom_smooth(method=lm, se=TRUE, size=0.3)+theme_bw()

ggplot(data=m5, aes(x=TCGA, y=us))+ geom_point()+ geom_smooth(method=lm, se=TRUE, size=0.3)+theme_bw()
