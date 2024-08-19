vsd <- "/mnt/cold1/snaketree/prj/DE_RNASeq/dataset/Biodiversa_up5_starOK_selected/vsd.tsv.gz"
vsd <- read.table(vsd, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
rownames(vsd) <- gsub("H_", "", rownames(vsd))

vsd <- as.data.frame(t(vsd))
vsd$type <- substr(rownames(vsd), 8, 10)
lmh <- vsd %>% filter(type == "LMH")
lmh$type <- NULL

lmh_model <- substr(rownames(lmh), 1, 7)

## gli lmh sono tutti senza repliche quindi si può girare il vsd ed è pronto

lmh <- as.data.frame(t(lmh))

## elimino subitooo gli LMX che non hanno un LMH
lmx <- vsd %>% filter(type == "LMX")
lmx$type <- NULL
lmx$model <- substr(rownames(lmx), 1 , 7)
lmx <- lmx %>% filter(model %in% lmh_model)


meta_f <- "/mnt/cold1/snaketree/prj/RNASeq_biod_metadata/dataset/july2020_starOK/selected_metadata_annot_final_nolinfo_nooutlier"
meta_df<-read.table(meta_f, quote = "", sep = "\t", header = TRUE)
meta_df$type <- gsub(".1", "", meta_df$type)
meta_df$type <- gsub(".2", "", meta_df$type)

meta <- meta_df %>% filter(type == "LMX_BASALE")

lmx_basali <- meta$sample_id_R

## ripetizioni per gli LMX

lmx$genealogy <- rownames(lmx)
lmx <- lmx %>% filter(genealogy %in% lmx_basali)
lmx$genealogy <- NULL

lmx$model <- substr(rownames(lmx), 1, 7)

modelli <- as.data.frame(table(lmx$model))
modelli <- modelli %>% filter(Freq == 1)
modellinorepliche <- as.character(modelli$Var1)

lmxnorepliche <- lmx %>% filter(model %in% modellinorepliche)
rownames(lmxnorepliche) <- paste0(lmxnorepliche$model, "LMX")
lmxnorepliche$model <- NULL
lmxnorepliche <- as.data.frame(t(lmxnorepliche))

lmx <- lmx %>% filter(!model %in% modellinorepliche)

unique_pairs <- unique(lmx$model)

vsd <- lmx

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

# calculate_correlation <- function(df) {
#   result <- cor.test(df[[1]], df[[2]])
#   estimate <- result$estimate
#   pvalue <- result$p.value
#   return(c(estimate = estimate, pvalue = pvalue))
# }
# 
# correlations <- sapply(t_df_list, calculate_correlation)

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
nuovo_vsd <- cbind(nuovo_vsd, lmxnorepliche)

colnames(lmh) <- substr(colnames(lmh), 1, 10)
### "CRC1451" "CRC0771" "CRC1605" non c'è LMX fuori dal cazzo il LMH
lmh$CRC1451LMH <- NULL
lmh$CRC0771LMH <- NULL
lmh$CRC1605LMH <- NULL

d <- cbind(nuovo_vsd, lmh)

### filtering expression data: we want high sd genes but not clear outliers / not expressed genes
### filter not expressed genes
rownames(d) <- d$genes
d$genes <- NULL
means <- apply(d, 1, mean)
med <- median(means)

de <- d[means > med,]

sds <- apply(de, 1, sd)
### there are some very high sd?
### Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
### 0.1086  0.3540  0.4702  0.5313  0.6411  3.7110

#noisy_genes_thr <- quantile(sds, 0.9) # tenere o no? non � una differenza cos� enorme come prima...toglilo

#de <- de[sds < noisy_genes_thr,]
#sds <- sds[sds < noisy_genes_thr]

### now we keep thet top 10% variable genes 
sds <- sds[order(-sds)]

n <- length(sds)
keep <- head(sds, round(0.10*n)) #prova con e senza
keep_genes <- names(keep)
desd <- de[rownames(de) %in% keep_genes,]

desd1 <- desd

lmx <- desd1[grepl('LMX', names(desd1),, fixed=TRUE)]
lmh <- desd1[grepl('LMH', names(desd1),, fixed=TRUE)]

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
colnames(lmx) <- substr(colnames(lmx),0,7)
colnames(lmh) <- substr(colnames(lmh),0,7)
clmx <- lmx[, names(lmh)]
clmh <- lmh[, names(clmx)]
### check also for genes
all_genes <- intersect(rownames(lmx), rownames(lmh))
clmx <- clmx[all_genes,]
clmh <- clmh[all_genes,]

if (all(colnames(clmx)!=colnames(clmh)) & all(rownames(clmx)!=rownames(clmh))) {
  stop('Brutto llama!')
}

res <- cor(clmx, clmh)

diag <- as.data.frame(diag(res))
diag$case <- rownames(diag)
diag <- diag[order(diag$`diag(res)`),]
diag$case <- NULL

#site <- "/mnt/cold1/snaketree/prj/hn/local/share/data/SampleTablenoCRLF.txt"
#site <- read.table(site, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
#site$model <- substr(site$SampleID, 1, 7)
#site <- site[,c("model", "Site")]

#diag$model <- rownames(diag)
# 
# merged <- merge(diag, site, by="model")
# merged <- merged[!duplicated(merged$model),]
# merged$`diag(res)` <- NULL
# merged <- merged[order(merged$Site),]
# rownames(merged) <- merged$model
# merged$model <- NULL
# 
res <- as.data.frame(res)

res <- res[,rownames(diag)]
res <- res[rownames(diag),]

pheatmap(res, cluster_rows = FALSE, cluster_cols = FALSE)

d <- res

minv <- -1
maxv <- 1
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


score <- "/mnt/cold1/snaketree/prj/hn/dataset/V1/RNAseq/xeno/tmm_lymphoma_scores.tsv.gz"
score <- read.table(score, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
