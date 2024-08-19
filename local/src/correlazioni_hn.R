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

plot <- mean

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

rownames(plot) <- plot$gene

for (i in rownames(plot)) {
  if (i %in% genidaeliminare) {
    plot[i, "eliminare"] <- "sì"
  } else {
    plot[i, "eliminare"] <- "no"
  }
}

ggplot(plot, aes(mean_lmx, mean_prh))+geom_point(aes(color=eliminare))+scale_color_manual(values=c("#000000", "#ff0000"))

#mean <- mean %>% filter(diffsegn > 0)

d <- cbind(nuovo_vsd2, prh)
d$mean_lmx <- NULL
d$mean_prh <- NULL

d$geni <- rownames(nuovo_vsd2)

d <- d %>% filter(!geni %in% genidaeliminare)

d$geni <- NULL

### filtering expression data: we want high sd genes but not clear outliers / not expressed genes
### filter not expressed genes
# rownames(d) <- d$genes
# d$genes <- NULL
# means <- apply(d, 1, mean)
# med <- median(means)
# 
# de <- d[means > med,]
# 
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

prx <- desd1[grepl('PRX', names(desd1), , fixed=TRUE)]
prh <- desd1[grepl('PRH', names(desd1), , fixed=TRUE)]


# TODO check so we do not have replicates or they are sorted 'correctly' here?
#> length(simo$LMX_lineage)
#[1] 55
#> length(simo$LMO_lineage)
#[1] 55
#> length(unique(substr(simo$LMO_lineage,0,7)))
#[1] 55
#> length(unique(substr(simo$LMX_lineage,0,7)))
#[1] 55
# Since we do not have replicates at the smodel level in Simo's right pairs we
# can safely substr right now and know that we do not have repetitions that will
# mess things up.
colnames(prx) <- substr(colnames(prx),0,7)
colnames(prh) <- substr(colnames(prh),0,7)
cprx <- prx[, names(prh)]
cprh <- prh[, names(cprx)]
### check also for genes
all_genes <- intersect(rownames(prx), rownames(prh))
cprx <- cprx[all_genes,]
cprh <- cprh[all_genes,]


if (all(colnames(cprx)!=colnames(cprh)) & all(rownames(cprx)!=rownames(cprh))) {
  stop('Brutto llama!')
}

res <- cor(cprx, cprh)
dmatr <- res
nondiagl <- res[lower.tri(dmatr, diag = FALSE)]
nondiagu <- res[upper.tri(dmatr, diag = FALSE)]
nondiag <- c(nondiagl, nondiagu)

diag <- as.data.frame(diag(res))
diag$case <- rownames(diag)
diag <- diag[order(diag$`diag(res)`),]
diag$case <- NULL

write.table(genidaeliminare, "/home/mferri/genidaeliminarehn.tsv", quote = FALSE, sep = "/t", row.names = FALSE, col.names = FALSE)

site <- "/mnt/cold1/snaketree/prj/hn/local/share/data/SampleTablenoCRLF.txt"
site <- read.table(site, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
site$model <- substr(site$SampleID, 1, 7)
site <- site[,c("model", "Site")]

diag$model <- rownames(diag)

merged <- merge(diag, site, by="model")
merged <- merged[!duplicated(merged$model),]
merged$`diag(res)` <- NULL
merged <- merged[order(merged$Site),]
rownames(merged) <- merged$model
merged$model <- NULL

res <- as.data.frame(res)

res <- res[,rownames(merged)]
res <- res[rownames(merged),]

pheatmap(res, cluster_rows = FALSE, cluster_cols = FALSE, annotation_row = merged)

d <- res

minv <- -0.5
maxv <- 0.5
#d[d < -4] <- -4
#d[d > 4] <- 4

neutral_value <- 0
#bk1 <- c(seq(minv-0.1,neutral_value-0.1,by=0.2),neutral_value-0.0999)
bk1 <- seq(minv-0.001, neutral_value-0.0009, length.out=224)
#bk2 <- c(neutral_value+0.001, seq(neutral_value+0.1,maxv+0.1,by=0.2))
bk2 <- seq(neutral_value+0.0001, maxv+0.001, length.out=224)
bk <- c(bk1, bk2)
my_palette <- c(colorRampPalette(colors = c("darkblue",
                                            "white"))(n = length(bk1)-1),
                "#e1e1e1", "#e1e1e1",
                c(colorRampPalette(colors = c("white", "darkred"))(n
                                                                   = length(bk2)-1)))
#pheatmap(matrix, breaks = seq(-rg, rg, length.out = 100))
pheatmap(d, cluster_rows = F, cluster_cols=F,
         breaks = bk, color=my_palette, na_col = "#FFFFFF")

g_h <- read.table("/mnt/cold1/snaketree/prj/hn/dataset/V1/RNAseq/umani/Cetuximab_Response_cutoff0.05-NonResponder.vs.Responder.deseq2.tsv",
                  quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
g_x <- read.table("/mnt/cold1/snaketree/prj/hn/dataset/V1/RNAseq/xeno/Cetuximab_Response_cutoff0.05-NonResponder.vs.Responder.deseq2.tsv",
                  quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE) 

h <- rownames(g_h)
x <- rownames(g_x)

intersection <- intersect(h, x)


diff <- as.data.frame(setdiff(h, intersection))

g_h <- read.table("/mnt/cold1/snaketree/prj/hn/dataset/V1/RNAseq/umani/Cetuximab_Response_cutoff0.05-NonResponder.vs.Responder.deseq2.tsv",
                  quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
g_h$genes <- rownames(g_h)

g_h <- g_h %>% filter(genes %in% diff$`setdiff(h, intersection)`)


score <- "/mnt/cold1/snaketree/prj/hn/dataset/V1/RNAseq/tmm_lymphoma_scores.tsv.gz"
score <- read.table(score, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
score$type <- substr(rownames(score), 8,10)
score <- score %>% filter(type == "PRH")
score$type <- NULL
score$model <- substr(rownames(score), 1,7)

diag$model <- rownames(diag)

merged <- merge(diag, score, by="model")
colnames(merged)[2]<- "cor"

ggplot(merged, aes(x=cor, y=Endothelial)) + geom_point()+ geom_text(aes(label = model), vjust = -1, hjust = 0, size=3)
ggplot(merged, aes(x=cor, y=CAF)) + geom_point()+ geom_text(aes(label = model), vjust = -1, hjust = 0, size=3)
ggplot(merged, aes(x=cor, y=Leucocyte)) + geom_point()+ geom_text(aes(label = model), vjust = -1, hjust = 0, size=3)


