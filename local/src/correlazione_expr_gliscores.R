## correlazione espressione gli e score fra

vsd <- "/mnt/cold1/snaketree/prj/hn/dataset/V1/RNAseq/xeno/vsd.tsv.gz"
vsd <- read.table(vsd, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
vsd$gene <- rownames(vsd)
vsd <- vsd %>% filter(gene == "GLI2")
vsd$gene <- NULL

vsd <- as.data.frame(t(vsd))
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

fra <- "/mnt/cold1/snaketree/prj/hn/local/share/data/gli_scores_fra.xlsx"
fra <- read.xlsx(fra)
colnames(fra) <- c("lmodel", "gli_score")
fra$model <- substr(fra$lmodel, 1,7)
fra$lmodel <- NULL
fra$gli_score[is.na(fra$gli_score)] <- 0

summaryfra <- fra %>%
  group_by(model) %>%
  summarize(gli_score = mean(gli_score))

summaryfra$model <- paste0(summaryfra$model, "PRX")

vsd2 <- nuovo_vsd
vsd2 <- as.data.frame(t(vsd2))
vsd2$model <- rownames(vsd2)

merged <- merge(vsd2, summaryfra, by="model")

cor.test(merged$GLI2, merged$gli_score)

ggplot(merged, aes(GLI2, gli_score))+geom_point()+geom_smooth(method = "lm")

jj <- "/mnt/cold1/snaketree/prj/hn/local/share/data/gli_scores_fra.xlsx"
jj <- read.xlsx(jj, sheet = 2)
jj <- jj[,c(1, 6)]
colnames(jj) <- c("lmodel", "gli_score")
jj$model <- substr(jj$lmodel, 1,7)
jj$lmodel <- NULL
jj$gli_score[is.na(jj$gli_score)] <- 0

summaryjj <- jj %>%
  group_by(model) %>%
  summarize(gli_score = mean(gli_score))

summaryjj$model <- paste0(summaryjj$model, "PRX")

mergedjj <- merge(vsd2, summaryjj, by="model")

cor.test(mergedjj$GLI2, mergedjj$gli_score)


samples <- read.table("/mnt/cold1/snaketree/prj/hn/local/share/data/SampleTablenoCRLF.txt", quote = "",
                      sep = "\t", header = TRUE, stringsAsFactors = FALSE)
table(samples$Cetuximab_Response)
