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


