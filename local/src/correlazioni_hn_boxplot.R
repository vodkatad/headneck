vsd <- "/mnt/cold1/snaketree/prj/hn/dataset/V1/RNAseq/vsd.tsv.gz"
vsd <- read.table(vsd, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
vsd <- as.data.frame(t(vsd))
vsd$type <- substr(rownames(vsd), 8, 10)
prh <- vsd %>% filter(type == "PRH")
prh$type <- NULL
prh <- as.data.frame(t(prh))
vsd <- vsd %>% filter(!type == "PRH")
vsd$type <- NULL
vsd$model <- substr(rownames(vsd), 1, 7)

modelli <- as.data.frame(table(vsd$model))
modelli <- modelli %>% filter(Freq == 1)
modellinorepliche <- as.character(modelli$Var1)

vsdnorepliche <- vsd %>% filter(model %in% modellinorepliche)
rownames(vsdnorepliche) <- paste0(vsdnorepliche$model, "PRX")
vsdnorepliche$model <- NULL
vsdnorepliche <- as.data.frame(t(vsdnorepliche))

vsd <- vsd %>% filter(!model %in% modellinorepliche)

unique_pairs <- unique(vsd$model)

# Split the dataframe based on unique pairs
df_list <- lapply(unique_pairs, function(combination) {
  subset_df <-  vsd[grep(combination, vsd$model), ]
  subset_df$model <- NULL
  return(subset_df)
})



t_df_list <- lapply(seq_along(df_list), function(i) {
  transposed_df <- as.data.frame(t(df_list[[i]]))
  return(transposed_df)
})

named_df_list <- list()

for (i in seq_along(t_df_list)) {
  first_col_name <- names(t_df_list[[i]])[1]
  name <- substr(first_col_name, 1, 7)
  named_df_list[[name]] <- t_df_list[[i]]
}


#prova <- t_df_list[[1]]

t_df_list <- named_df_list

calculate_correlation <- function(df) {
  result <- cor.test(df[[1]], df[[2]])
  estimate <- result$estimate
  pvalue <- result$p.value
  return(c(estimate = estimate, pvalue = pvalue))
}

correlations <- sapply(t_df_list, calculate_correlation)

df_list <- lapply(t_df_list, function(df) {
  nuovo_nome_colonna <- substr(names(df)[1], 1, 10)
  if (ncol(df) == 3) {
    df <- mutate(df, !!nuovo_nome_colonna := rowMeans(df[, 1:3]))
    df <- subset(df, select = -c(1, 2,3))
  } else if (ncol(df) == 2) {
    df <- mutate(df, !!nuovo_nome_colonna := rowMeans(df[, 1:2]))
    df <- subset(df, select = -c(1, 2))
  }
  return(df)
})

nuovo_vsd <- do.call(cbind, df_list)
nuovo_vsd <- cbind(nuovo_vsd, vsdnorepliche)

## ma tutti i prx hanno il prh? mi sale il dubbio visti i numeri

model_prx <- substr(colnames(nuovo_vsd), 1, 7)
model_prh <- substr(colnames(prh), 1, 7)

both <- intersect(model_prx, model_prh)
soloprx <- setdiff(model_prx, model_prh)

## manca il 98 che dovrebbe esserci ma non c'è neanche nei counts di shauna -> eliminato per QC

nuovo_vsd2 <- nuovo_vsd
colnames(nuovo_vsd2) <- substr(colnames(nuovo_vsd2), 1, 7)
nuovo_vsd2 <- as.data.frame(t(nuovo_vsd2))
nuovo_vsd2$model <- rownames(nuovo_vsd2)
nuovo_vsd2 <- nuovo_vsd2 %>% filter(model %in% both)
nuovo_vsd2$model <- NULL
nuovo_vsd2 <- as.data.frame(t(nuovo_vsd2))
colnames(nuovo_vsd2) <- paste0(colnames(nuovo_vsd2), "PRX")

## tolgo ai prx tutti i geni espressi nei prh che non sono espressi

nuovo_vsd2$mean_lmx <- rowMeans(nuovo_vsd2)
prh$mean_prh <- rowMeans(prh)

mean_lmx <- nuovo_vsd2[,c(1,37)]
mean_lmx$HNC0002PRX <- NULL
mean_lmx$gene <- rownames(mean_lmx)

mean_prh <- prh[,c(1,37)]
mean_prh$HNC0002PRH00 <- NULL
mean_prh$gene <- rownames(mean_prh)

mean <- merge(mean_lmx, mean_prh, by="gene")

ggplot(mean, aes(mean_lmx, mean_prh))+geom_point()

med_lmx <- median(mean$mean_lmx)
med_prh <- median(mean$mean_prh)

mean <- mean %>% filter(!mean_prh < med_prh)
mean <- mean %>% filter(!mean_lmx > med_lmx)

ggplot(mean, aes(mean_lmx, mean_prh))+geom_point()

rownames(mean) <- mean$gene
mean$gene <- NULL
mean$mean_tot <- rowMeans(mean)

#mean$diff <- abs(mean$mean_lmx - mean$mean_prh)
#mean$diffsegn <- mean$mean_lmx-mean$mean_prh

mean$fc <- mean$mean_prh/mean$mean_lmx
mean <- mean %>% filter(fc >= 1.5)

#mean <- mean %>% filter(diff >=3)

#q <- quantile(mean$diffsegn[mean$diffsegn<0], 0.1)
# mean <- mean %>% filter(diffsegn <= q)
# 
# ggplot(mean, aes(mean_lmx, mean_prh))+geom_point()
# 
# ## prendere il valore del primo quartile degli lmx (3.44)
# 
# for (i in rownames(mean)) {
#   if (mean[i, "mean_lmx"] > 3.44) {
#     mean[i, "lmx_alto"] <- "sì"
#   } else {
#     mean[i, "lmx_alto"] <- "no"
#   }
# }
# 
# #mean$diff <- abs(mean$mean_lmx - mean$mean_prh)
# 
# mean <- mean %>% filter(lmx_alto == "no")
# 
# ggplot(mean, aes(mean_lmx, mean_prh))+geom_point()

genidaeliminare <- rownames(mean)

#mean <- mean %>% filter(diffsegn > 0)

d <- cbind(nuovo_vsd2, prh)

d$geni <- rownames(nuovo_vsd2)

d <- d %>% filter(!geni %in% genidaeliminare)

d$geni <- NULL
d$mean_lmx <- NULL
d$mean_prh <- NULL

de <- d
sds <- apply(de, 1, sd)
# ### there are some very high sd?
# ### Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# ### 0.1086  0.3540  0.4702  0.5313  0.6411  3.7110
# 
# #noisy_genes_thr <- quantile(sds, 0.9) # tenere o no? non � una differenza cos� enorme come prima...toglilo
# 
# #de <- de[sds < noisy_genes_thr,]
# #sds <- sds[sds < noisy_genes_thr]
# 

### now we keep thet top 10% variable genes 
sds <- sds[order(-sds)]

n <- length(sds)
keep <- head(sds, round(0.10*n)) #prova con e senza
keep_genes <- names(keep)
desd <- de[rownames(de) %in% keep_genes,]

desd1 <- desd

samples <- names(desd1)

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

method <- 'pearson'

get_cor <- function(pair, data, method="pearson") {
  # we are sure that we will always find a match because we have build samples coming from our data itself
  expr_one <- data[, names(data) == pair[1]]
  expr_two <- data[, names(data) == pair[2]]
  corrs <- cor.test(expr_one, expr_two, method=method)
  res <- c(corrs$estimate, corrs$p.value)
  attributes(res) <- NULL
  return(res)
}

right_pairs <- as.data.frame(t(apply(pairsdf, 1, get_cor, desd1, method)))

pairsdf_wrong <- pairsdf

for (i in seq(1, nrow(pairsdf_wrong))) {
  model <- substr(pairsdf_wrong[i, 'xeno'],0,7)
  random_model <- model
  while (random_model == model) {
    random <- sample(pairsdf$human[-i],1)
    random_model <- substr(random,0,7)
  }
  pairsdf_wrong[i,'human'] <- random
}

wrong_pairs <- as.data.frame(t(apply(pairsdf_wrong, 1, get_cor, desd1, method)))

colnames(right_pairs) <- c(method, 'pvalue')
colnames(wrong_pairs) <- c(method, 'pvalue')
wrong_pairs$pairing <- 'wrong'
right_pairs$pairing <- 'right'

m <- rbind(right_pairs, wrong_pairs)

ggplot(m, aes_string(x="pairing", y=method, fill="pairing"))+geom_violin()+theme_bw()+scale_fill_manual(values=c("darkgoldenrod","darkgreen"))+theme(text=element_text(size=15))+geom_signif(comparisons=list(c('right','wrong')))

ggplot(m, aes_string(x="pairing", y=method, fill="pairing"))+geom_boxplot(outlier.shape = NA)+geom_jitter()+theme_bw()+scale_fill_manual(values=c("darkgoldenrod","darkgreen"))+theme(text=element_text(size=15))+geom_signif(comparisons=list(c('right','wrong')))
