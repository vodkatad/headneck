set.seed(42)
library(ggplot2)
library(precrec)
library(pheatmap)

d <- read.table('/mnt/cold1/snaketree/prj/hn/local/share/data/cks_per_ROC.tsv', sep="\t", header=T, stringsAsFactors = F) 

d$cetuxi <- factor(d$cetuxi, levels=c('s', 'r'))
ggplot(data=d, aes(x=cetuxi, y=score_pos, fill=cetuxi))+geom_boxplot(outlier.shape=NA)+geom_jitter(aes(color=oral), height=0)+theme_bw(base_size=18)+
  scale_fill_manual(values=c('blue', 'red'))+scale_color_manual(values=c('black', 'orange'))


ggplot(data=d, aes(x=cetuxi, y=score_neg, fill=cetuxi))+geom_boxplot(outlier.shape=NA)+geom_jitter(aes(color=oral), height=0)+theme_bw(base_size=18)+
  scale_fill_manual(values=c('blue', 'red'))+scale_color_manual(values=c('black', 'orange'))

precrec_obj <- evalmod(scores = d$score_pos, labels = d$cetuxi, posclass='s', mode="rocprc", ties='equiv')
#pdf(auc_all_f, height=2.5, width=2.5)
autoplot(precrec_obj, curvetype = c("ROC"))
#graphics.off()

precrec_obj <- evalmod(scores = d$score_neg, labels = d$cetuxi, posclass='s', mode="rocprc", ties='equiv')
#pdf(auc_all_f, height=2.5, width=2.5)
autoplot(precrec_obj, curvetype = c("ROC"))
#graphics.off()



# define optimal threshold and try to apply on patients?
get_sensitivity_specificity <- function(thr, df) {
  p <- nrow(df[df$labels == 1,])
  tp <- nrow(df[df$scores > thr & df$labels==1,])
  n <- nrow(df[df$labels == 0,])
  tn <- nrow(df[df$scores <= thr & df$labels==0,])
  fn <- nrow(df[df$scores <= thr & df$labels==1,])
  fp <- nrow(df[df$scores > thr & df$labels==0,])
  return(c(tp/p, tn/n, tp, tn, fp, fn, fp/(fp+tp)))
}

compute_thr <- function(scores, labels) {
  df <- data.frame(scores=scores, labels=labels)
  intervals <- cut(scores, breaks=10)
  max_i <- levels(intervals)
  max_i <- gsub('(','', max_i, fixed=TRUE)
  max_i <- gsub(']','', max_i, fixed=TRUE)
  upper_bs <- sapply(strsplit(max_i, ','), '[[', 2)
  upper_bs <- as.numeric(upper_bs)
  upper_bs <- upper_bs[order(-upper_bs)]
  res <- as.data.frame(t(sapply(upper_bs, get_sensitivity_specificity, df)))
  colnames(res) <- c('sensitivity', 'specificity', 'tp', 'tn', 'fp', 'fn', 'fdr')
  res$thr <- upper_bs
  return(res)
}

sens <- compute_thr(d$score_neg, ifelse(d$cetuxi=='s', 1, 0))


library(pROC)
my_curve <- roc(predictor=d$score_pos, response=ifelse(d$cetuxi=='s', 1, 0))
plot(my_curve, print.thres=TRUE)

coords(my_curve, "best", best.method="closest.topleft")

cposd <- coords(my_curve, "best")

cpos <- cposd[1,1]

my_curve <- roc(predictor=d$score_neg, response=ifelse(d$cetuxi=='s', 1, 0))
plot(my_curve, print.thres=TRUE)

cnegd <- coords(my_curve, "best")

cneg <- cnegd[1,1]

coords(my_curve, "best", best.method="closest.topleft")

dp <- read.table('/mnt/cold1/snaketree/prj/hn/local/share/data/cks_per_ROC_patients.tsv', sep="\t", header=T, stringsAsFactors = F) 

dp$labels <- ifelse(dp$cetuxi=='s', 1, 0)
dp$scores <- dp$score_pos
get_sensitivity_specificity(cpos, dp)
dp$scores <- dp$score_neg
get_sensitivity_specificity(cneg, dp)

### def version only score_pos 
set.seed(42)
library(ggplot2)
library(precrec)
library(pheatmap)
library(RColorBrewer)

d <- read.table('/mnt/cold1/snaketree/prj/hn/local/share/data/cks_per_ROC_3.tsv', sep="\t", header=T, stringsAsFactors = F) 

cet <- read.table('/mnt/cold1/snaketree/prj/hn/local/share/data/def_cohort/MaryKate.tsv', sep="\t", header=T, stringsAsFactors = F)
cet$smodel <- substr(cet$X, 0, 7)
cet$X <- NULL
cet <- cet[!duplicated(cet),]
length(unique(cet$smodel))
dim(cet)

d$smodel <- d$Genealogy
dim(d)

# sanity check of RNaseq responses used by mk and new colors by Fra
m <-  merge(d, cet, by="smodel")
dim(m)
table(m$cet, m$Definitive.resp)

d$cetuxi <- factor(d$cet, levels=c('s', 'r'))
ggplot(data=d, aes(x=cetuxi, y=score_pos, fill=cetuxi))+geom_boxplot(outlier.shape=NA)+geom_jitter(height=0)+theme_bw(base_size=18)+
  scale_fill_manual(values=c('blue', 'red'))
dim(d)

precrec_obj <- evalmod(scores = d$score_pos, labels = d$cetuxi, posclass='s', mode="rocprc", ties='equiv')
#pdf(auc_all_f, height=2.5, width=2.5)
autoplot(precrec_obj, curvetype = c("ROC"))

library(pROC)
my_curve <- roc(predictor=d$score_pos, response=ifelse(d$cetuxi=='s', 1, 0))
pdf('~/ck_roc_fig5.pdf')
plot(my_curve, print.thres=TRUE)
graphics.off()


sd <- data.frame(my_curve$thresholds, my_curve$sensitivities, my_curve$specificities)
write.table(sd, '~/hn_roc_sd_figure5.tsv', sep="\t", row.names = F)

ci(my_curve)

coords(my_curve, "best", best.method="closest.topleft")

cposd <- coords(my_curve, "best")

cpos <- cposd[1,1]

d$labels <- ifelse(d$cetuxi=='s', 1, 0)
d$scores <- d$score_pos
get_sensitivity_specificity(cpos, d)


dp <- read.table('/mnt/cold1/snaketree/prj/hn/local/share/data/cks_per_ROC_patients_3.tsv', sep="\t", header=T, stringsAsFactors = F) 

dp$labels <- ifelse(dp$cetuxi=='PR', 1, 0)
dp$cet <- ifelse(dp$cetuxi=='PR', 's', 'r')

dp$scores <- dp$score_pos
get_sensitivity_specificity(cpos, dp)

dp$cet <- factor(dp$cet, levels=c('s', 'r'))
ggplot(data=dp, aes(x=cet, y=score_pos, fill=cet))+geom_boxplot(outlier.shape=NA)+geom_jitter(height=0)+theme_bw(base_size=18)+
  scale_fill_manual(values=c('blue', 'red'))+geom_hline(yintercept=4.5)

dp$predict <- ifelse(dp$score_pos > cpos, 's', 'r')

cfm <- as.matrix(table(dp$cet, dp$predict)) # on rows we have cet reality on columns the prediction

pheatmap(cfm, cluster_cols = F, cluster_rows = F, display_numbers = T, fontsize_number=20, 
         color = colorRampPalette(brewer.pal(n = 3, name ="YlGnBu"))(100), angle_col=0,
         labels_row=c('Non Responder', 'Responder'), labels_col=c('Predicted Non Responder', 'Predicted Responder'),
         filename='~/ck_confusionpatients_fig5G.pdf')#, labels_row = 'Response', labels_col='Prediction')

### sept 2026 ###############################################################

library(pROC)
set.seed(42)
library(ggplot2)
library(precrec)
library(pheatmap)

dp <- read.table('/mnt/cold1/snaketree/prj/hn/local/share/data/ck_sept2026.txt', sep="\t", header=T, stringsAsFactors = F) 
dp$labels <- ifelse(dp$`Response.to.CETUXIMAB`=='RESPONDER', 1, 0)
dp$cetuxi <- ifelse(dp$`Response.to.CETUXIMAB`=='RESPONDER', 's', 'r')

### sanity check vs MK
cet <- read.table('/mnt/cold1/snaketree/prj/hn/local/share/data/def_cohort/MaryKate.tsv', sep="\t", header=T, stringsAsFactors = F)
cet$smodel <- substr(cet$X, 0, 7)
cet$X <- NULL
cet <- cet[!duplicated(cet),]
length(unique(cet$smodel))
dim(cet)

dp$smodel <- dp$ID
dim(dp)

# sanity check of RNaseq responses used by mk and new colors by Fra
m <-  merge(dp, cet, by="smodel")
dim(m)
table(m$cet, m$Definitive.resp)

### sanity check vs old ck
ck <- read.table('/mnt/cold1/snaketree/prj/hn/local/share/data/cks_per_ROC_3.tsv', sep="\t", header=T, stringsAsFactors = F)
mm <- merge(ck,dp, by.y='smodel', by.x='Genealogy')
nrow(ck)
nrow(dp)
nrow(mm)
all(mm$CK1.IHC.SCORE==mm$CK1)
all(mm$CK5.IHC.SCORE==mm$CK5)
all(mm$CK10.IHC.SCORE==mm$CK10)

###

ggplot(data=dp, aes(x=cetuxi, y=overall.score, fill=cetuxi))+geom_boxplot(outlier.shape=NA)+geom_jitter(height=0)+theme_bw(base_size=18)+
  scale_fill_manual(values=c('blue', 'red'))


precrec_obj <- evalmod(scores = dp$overall.score, labels = dp$cetuxi, posclass='s', mode="rocprc", ties='equiv')
#pdf(auc_all_f, height=2.5, width=2.5)
autoplot(precrec_obj, curvetype = c("ROC"))

my_curve <- roc(predictor=dp$overall.score, response=ifelse(dp$cetuxi=='s', 1, 0))
plot(my_curve, print.thres=TRUE)
ci(my_curve)
my_curve$auc
coords(my_curve, "best", best.method="closest.topleft")
coords(my_curve, "best", best.method="youden")
plot(my_curve, print.thres='all')
plot(my_curve, print.thres='local maximas')

ggplot(data=dp, aes(x=cetuxi, y=pos.score, fill=cetuxi))+geom_boxplot(outlier.shape=NA)+geom_jitter(height=0)+theme_bw(base_size=18)+
  scale_fill_manual(values=c('blue', 'red'))
ggplot(data=dp, aes(x=cetuxi, y=-neg.score, fill=cetuxi))+geom_boxplot(outlier.shape=NA)+geom_jitter(height=0)+theme_bw(base_size=18)+
  scale_fill_manual(values=c('blue', 'red'))
precrec_obj <- evalmod(scores = dp$pos.score, labels = dp$cetuxi, posclass='s', mode="rocprc", ties='equiv')
#pdf(auc_all_f, height=2.5, width=2.5)
autoplot(precrec_obj, curvetype = c("ROC"))

my_curve <- roc(predictor=dp$pos.score, response=ifelse(dp$cetuxi=='s', 1, 0))
plot(my_curve, print.thres=TRUE)
ci(my_curve)
my_curve$auc

dp$score_old <- dp$pos.score - dp$CK8.IHC.SCORE
my_curve <- roc(predictor=dp$score_old, response=ifelse(dp$cetuxi=='s', 1, 0))
plot(my_curve, print.thres=TRUE)
ci(my_curve)
my_curve$auc

grade <- read.table('/mnt/cold1/snaketree/prj/hn/local/share/data/Supplementary_Data_1.txt', sep="\t", header=T, stringsAsFactors = F) 

mg <- merge(dp, grade, by.x='smodel', by.y='Case.ID')
mgg <- mg[mg$HISTOLOGICAL.GRADE!='',]
nrow(mgg)

fisher.test(table(mgg$HISTOLOGICAL.GRADE, mgg$Response.to.CETUXIMAB))
ggplot(data=mgg, aes(x=HISTOLOGICAL.GRADE, y=overall.score, fill=cetuxi))+geom_boxplot(outlier.shape=NA)+geom_jitter(height=0)+theme_bw(base_size=18)+
  scale_fill_manual(values=c('blue', 'red'))


ggplot(data=mgg, aes(x=HISTOLOGICAL.GRADE, y=overall.score))+geom_boxplot(outlier.shape=NA)+geom_jitter(height=0)+theme_bw(base_size=18)


wilcox.test(mgg[mgg$HISTOLOGICAL.GRADE=='g3', 'overall.score'], mgg[mgg$HISTOLOGICAL.GRADE=='g2' , 'overall.score'])
wilcox.test(mgg[mgg$HISTOLOGICAL.GRADE=='g3', 'overall.score'], mgg[mgg$HISTOLOGICAL.GRADE!='g3' , 'overall.score'])

d <- mgg[,c('CK1.IHC.SCORE','CK5.IHC.SCORE','CK8.IHC.SCORE','CK10.IHC.SCORE','CK18.IHC.SCORE', 'pos.score', 'overall.score')]

annot_rows <- mgg[, c('Response.to.CETUXIMAB', 'HISTOLOGICAL.GRADE', 'Site.of.Primary')]

library(pheatmap)
minv <- min(d)
maxv <- max(d)
neutral_value <- mean(c(minv, maxv))

bk1 <- c(seq(minv-0.1,neutral_value-0.1,by=0.2),neutral_value-0.0999)
bk2 <- c(neutral_value+0.001, seq(neutral_value+0.1,maxv+0.1,by=0.2))
bk <- c(bk1, bk2)
my_palette <- c(colorRampPalette(colors = c("darkblue",
                                            "lightblue"))(n = length(bk1)-1),
                "#FFFFFF", #"snow1",
                c(colorRampPalette(colors = c("tomato1", "darkred"))(n
                                                                     = length(bk2)-1)))

pheatmap(d, annotation_row=annot_rows, cluster_rows = T, cluster_cols=T,
         breaks = bk, color=my_palette)


mo <- glm(formula="labels~overall.score+HISTOLOGICAL.GRADE", data=mg, family='binomial')
summary(mo)
compute_thr(mg$overall.score, mg$labels)
table(mg[mg$overall.score > 1.5 & mg$Response.to.CETUXIMAB=="RESPONDER",c('HISTOLOGICAL.GRADE')])
table(mg[mg$overall.score <= 1.5 & mg$Response.to.CETUXIMAB=="RESPONDER",c('HISTOLOGICAL.GRADE')])
table(mg[mg$overall.score <= 1.5 & mg$Response.to.CETUXIMAB=="NON RESPONDER",c('HISTOLOGICAL.GRADE')])

