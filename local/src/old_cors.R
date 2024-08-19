library(ggplot2)
library(ggsignif)
setwd('/mnt/cold1/snaketree/prj/hn/local/share/data')

file <- 'DESeq_normalised_counts_Shauna3009.txt'

set.seed(42)

## specific input management
d <- read.table(file, sep="\t", header=TRUE, row.names = 1)
names(d) <- gsub('PN0039B.', '', names(d), fixed=T) # uffa per far cosi` ogni volta mi dimentico se siano le col o le row...
#> unique(sapply(strsplit(colnames(d),"_"), '[[', 2))
#[1] "graft"
names(d) <- gsub('_graft', '', names(d), fixed=T)


## filtering expression data: we want high sd genes but not clear outliers / not expressed genes
# filter not expressed genes
means <- apply(d, 1, mean)
med <- median(means)

de <- d[means > med,]

sds <- apply(de, 1, sd)
# there are some very high sd?
#    Min.  1st Qu.   Median     Mean  3rd Qu.     Max. 
#0.9      5.6     28.9    354.0    149.4 860110.9 
#noisy_genes_thr <- quantile(sds, 0.9)

#de <- de[sds < noisy_genes_thr,]
#sds <- sds[sds < noisy_genes_thr]

# now we keep the top 10% variable genes 
sds <- sds[order(-sds)]

n <- length(sds)
keep <- head(sds, round(0.10*n))
keep_genes <- names(keep)
desd <- de[rownames(de) %in% keep_genes,]


## build correspondences of human/pdx
# ids here are "HNC0002PRH00" "HNC0002PRX0A" 
# we work separately on 0A and 0B or together? Start with together
samples <- names(desd)
# remove the  technical replicates, ends in 2
replicates <- samples[grepl('2$', samples)]
samples <- samples[!samples %in% replicates]


### from now on it could be a generalized script that gets the human/xeno model suffixes, length of common model name and computes correlations
# in an ideal pipeline the previous steps should be in rules with the different steps (in order to be able to check the results)
# and the correlations computed with this last 'generic' script.
# but...for CRC the picture is more complex therefore we should think about the structure before starting (LMX/LMH, more than only XA/XB double matches)
# and evaluate if we will ever obtain a deeply general solution
human <- samples[grepl('PRH', samples, fixed=TRUE)]

#pairs <- data.frame(human=human, xeno=rep('', length(human))
find_xeno <- function(h, s) {
  model <- substr(h, 0 , 7)
  model <- paste0(model, 'PRX')
  x <- s[grepl(model,s,fixed=TRUE)]
  res <- data.frame(human=rep(h, length(x)), xeno=x)
  return(res)
}

pairs <- lapply(human, find_xeno, samples)
pairsdf <- do.call(rbind, pairs)

# all 40 human samples have a paired pdx
#> length(human); length(unique(pairsdf$human))
#[1] 40
#[1] 40
#
#> table(substr(pairsdf$xeno,nchar(pairsdf$xeno), nchar(pairsdf$xeno)+1))
#A  B 
#30 27
#

# 17 xenos have both XA and XB, 23 only 1
# > tt <- table(substr(pairsdf$xeno,0, 7))
# > tt[tt==2] 
# > length(tt[tt==2])
# [1] 17
# > length(tt[tt==1])
# [1] 23

method <- 'pearson'

get_cor <- function(pair, data, method="spearman") {
  # we are sure that we will always find a match because we have build samples coming from our data itself
  expr_one <- data[, names(data) == pair[1]]
  expr_two <- data[, names(data) == pair[2]]
  corrs <- cor.test(expr_one, expr_two, method=method)
  res <- c(corrs$estimate, corrs$p.value)
  attributes(res) <- NULL
  return(res)
}

right_pairs <- as.data.frame(t(apply(pairsdf, 1, get_cor, desd, method)))

pairsdf_wrong <- pairsdf # we shuffle also A and B in this way
#pairsdf_wrong$human <- sample(pairsdf$human,nrow(pairsdf)) # this ends up in assigning somethimes the right pairs anyhow
for (i in seq(1, nrow(pairsdf_wrong))) {
  model <- substr(pairsdf_wrong[i, 'xeno'],0,7)
  random_model <- model
  while (random_model == model) {
    random <- sample(pairsdf$human[-i],1)
    random_model <- substr(random,0,7)
  }
  pairsdf_wrong[i,'human'] <- random
}

wrong_pairs <- as.data.frame(t(apply(pairsdf_wrong, 1, get_cor, desd, method)))

colnames(right_pairs) <- c(method, 'pvalue')
colnames(wrong_pairs) <- c(method, 'pvalue')
wrong_pairs$pairing <- 'wrong'
right_pairs$pairing <- 'right'

m <- rbind(right_pairs, wrong_pairs)

ggplot(m, aes_string(x="pairing", y=method, fill="pairing"))+geom_violin()+theme_bw()+scale_fill_manual(values=c("darkgoldenrod","darkgreen"))+theme(text=element_text(size=15))+geom_signif(comparisons=list(c('right','wrong')))

ggplot(m, aes_string(x="pairing", y=method, fill="pairing"))+geom_boxplot(outlier.shape = NA)+geom_jitter(height=NULL)+theme_bw()+scale_fill_manual(values=c("darkgoldenrod","darkgreen"))+theme(text=element_text(size=15))+geom_signif(comparisons=list(c('right','wrong')))


#> check <- cbind(right_pairs, pairsdf)
#> check[check$spearman < 0.5,]
#spearman        pvalue pairing        human         xeno
#22 0.4071286 2.380094e-109   right HNC0050PRH00 HNC0050PRX0A
#36 0.4147276 8.390562e-114   right HNC0098PRH00 HNC0098PRX0A
#38 0.4911699 1.544005e-165   right HNC0107PRH00 HNC0107PRX0A
#> median(check$spearman)
#[1] 0.8180192

#Moreover, in 88% of the cases, the best correlate for each PDX was its matched original counterpart 

# filter only all zero genes and keep all top sds:


sums <- rowSums(d)
de <- d[sums != 0,]

sds <- apply(de, 1, sd)

# now we keep thet top 10% variable genes 
sds <- sds[order(-sds)]

n <- length(sds)
keep <- head(sds, round(0.10*n))
keep_genes <- names(keep)

desd <- de[rownames(de) %in% keep_genes,]

method <- 'pearson'

right_pairs <- as.data.frame(t(apply(pairsdf, 1, get_cor, desd, method)))
wrong_pairs <- as.data.frame(t(apply(pairsdf_wrong, 1, get_cor, desd, method)))


colnames(right_pairs) <- c(method, 'pvalue')
colnames(wrong_pairs) <- c(method, 'pvalue')
wrong_pairs$pairing <- 'wrong'
right_pairs$pairing <- 'right'

m <- rbind(right_pairs, wrong_pairs)

ggplot(m, aes_string(x="pairing", y=method, fill="pairing"))+geom_boxplot(outlier.shape = NA)+geom_jitter()+theme_bw()+scale_fill_manual(values=c("darkgoldenrod","darkgreen"))+theme(text=element_text(size=15))+geom_signif(comparisons=list(c('right','wrong')))


#
p <- cbind(pairsdf, right_pairs)

