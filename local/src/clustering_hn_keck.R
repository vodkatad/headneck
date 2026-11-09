library(networkD3)
library(htmlwidgets)

vsd <- "/mnt/cold1/snaketree/prj/hn/dataset/V1/RNAseq_likeMK/HUMAN/fpkm.tsv.gz"
vsd <- read.table(vsd, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
vsd <- log2(vsd+1)

vsd2 <- as.data.frame(t(vsd))
centered_data <- scale(vsd2, scale = FALSE)
vsd3 <- as.data.frame(t(centered_data))

centroids <- "/mnt/cold1/snaketree/prj/hn/local/share/data/centroids_hn.xlsx"
centroids <- read.xlsx(centroids, colNames = TRUE)
colnames(centroids) <- centroids[1,]
centroids <- centroids[-1,]
centroids$IMS <- as.numeric(centroids$IMS)
centroids$BA <- as.numeric(centroids$BA)
centroids$CL <- as.numeric(centroids$CL)
ims <- centroids
ims$BA <- NULL
ims$CL <- NULL
ba <- centroids
ba$IMS <- NULL
ba$CL <- NULL
cl <- centroids
cl$IMS <- NULL
cl$BA <- NULL

vsd3 <- as.data.frame(vsd3)
vsd3$GENE <- rownames(vsd3)

imsvsd <- merge(ims, vsd3, by="GENE")

# Initialize a dataframe to store results
ims_correlations <- data.frame(
  column = character(),
  estimate = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

rownames(imsvsd) <- imsvsd$GENE
imsvsd$GENE <- NULL

for (i in colnames(imsvsd)[2:48]) {
  # Perform correlation test with column 21
  j <- imsvsd[,i]
  cor_test <- cor.test(j, imsvsd$IMS, method = "pearson")
  
  # Extract estimate and p-value
  estimate <- cor_test$estimate
  p_value <- cor_test$p.value
  
  # Add results to the dataframe
  ims_correlations <- rbind(ims_correlations, data.frame(
    column = i,
    estimate = estimate,
    p_value = p_value,
    stringsAsFactors = FALSE
  ))
}

bavsd <- merge(ba, vsd3, by="GENE")

# Initialize a dataframe to store results
ba_correlations <- data.frame(
  column = character(),
  estimate = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

rownames(bavsd) <- bavsd$GENE
bavsd$GENE <- NULL

# Loop through columns 1 to 20
for (i in colnames(bavsd)[2:48]) {
  # Perform correlation test with column 21
  j <- bavsd[,i]
  cor_test <- cor.test(j, bavsd$BA, method = "pearson")
  
  # Extract estimate and p-value
  estimate <- cor_test$estimate
  p_value <- cor_test$p.value
  
  # Add results to the dataframe
  ba_correlations <- rbind(ba_correlations, data.frame(
    column = i,
    estimate = estimate,
    p_value = p_value,
    stringsAsFactors = FALSE
  ))
}

clvsd <- merge(cl, vsd3, by="GENE")

# Initialize a dataframe to store results
cl_correlations <- data.frame(
  column = character(),
  estimate = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

rownames(clvsd) <- clvsd$GENE
clvsd$GENE <- NULL

# Loop through columns 1 to 20
for (i in colnames(clvsd)[2:48]) {
  # Perform correlation test with column 21
  j <- clvsd[,i]
  cor_test <- cor.test(j, clvsd$CL, method = "pearson")
  
  # Extract estimate and p-value
  estimate <- cor_test$estimate
  p_value <- cor_test$p.value
  
  # Add results to the dataframe
  cl_correlations <- rbind(cl_correlations, data.frame(
    column = i,
    estimate = estimate,
    p_value = p_value,
    stringsAsFactors = FALSE
  ))
}

names(ba_correlations)[names(ba_correlations) =="estimate"] <- "cor_BA"
names(cl_correlations)[names(cl_correlations) =="estimate"] <- "cor_CL"
names(ims_correlations)[names(ims_correlations) =="estimate"] <- "cor_IMS"

ba_correlations[ ba_correlations$p_value > 0.05, 'cor_BA'] <- -1
cl_correlations[ cl_correlations$p_value > 0.05, 'cor_CL'] <- -1
ims_correlations[ ims_correlations$p_value > 0.05, 'cor_IMS'] <- -1

res <- merge(ba_correlations, cl_correlations, by="column")
res$p_value.x <- NULL
res$p_value.y <- NULL
res <- merge(res, ims_correlations, by="column")
res$p_value <- NULL


dp <- res
rownames(dp) <- dp$column
dp$column <- NULL
pheatmap(dp)

table(colnames(dp)[apply(dp, 1, which.max)])

max_col <- apply(dp, 1, function(x) colnames(dp)[which.max(x)])

assign_df <- data.frame(
  model = rownames(dp),
  assigned_cor = max_col,
  stringsAsFactors = FALSE
)


## prendere TCGA e vedere se i numeri sono confrontabili con il paper

# tcga <- "/home/mferri/TCGA.HNSC.sampleMap_HiSeqV2.gz"
# tcga <- read.table(tcga, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
# rownames(tcga) <- tcga$sample
# tcga$sample <- NULL
# 
# tcga <- as.data.frame(t(tcga))
# centered_data <- scale(tcga, scale = FALSE)
# tcga <- as.data.frame(t(centered_data))
# tcga <- scale(tcga, scale = TRUE)
# tcga <- as.data.frame(tcga)
# tcga$GENE <- rownames(tcga)
# 
# 
# imsvsd <- merge(ims, tcga, by="GENE")
# 
# # Initialize a dataframe to store results
# ims_correlations <- data.frame(
#   column = character(),
#   estimate = numeric(),
#   p_value = numeric(),
#   stringsAsFactors = FALSE
# )
# 
# rownames(imsvsd) <- imsvsd$GENE
# imsvsd$GENE <- NULL
# 
# 
# for (i in colnames(imsvsd)[2:567]) {
#   # Perform correlation test with column 21
#   j <- imsvsd[,i]
#   cor_test <- cor.test(j, imsvsd$IMS, method = "pearson")
#   
#   # Extract estimate and p-value
#   estimate <- cor_test$estimate
#   p_value <- cor_test$p.value
#   
#   # Add results to the dataframe
#   ims_correlations <- rbind(ims_correlations, data.frame(
#     column = i,
#     estimate = estimate,
#     p_value = p_value,
#     stringsAsFactors = FALSE
#   ))
# }
# 
# bavsd <- merge(ba, tcga, by="GENE")
# 
# # Initialize a dataframe to store results
# ba_correlations <- data.frame(
#   column = character(),
#   estimate = numeric(),
#   p_value = numeric(),
#   stringsAsFactors = FALSE
# )
# 
# rownames(bavsd) <- bavsd$GENE
# bavsd$GENE <- NULL
# 
# 
# for (i in colnames(bavsd)[2:567]) {
#   # Perform correlation test with column 21
#   j <- bavsd[,i]
#   cor_test <- cor.test(j, bavsd$BA, method = "pearson")
#   
#   # Extract estimate and p-value
#   estimate <- cor_test$estimate
#   p_value <- cor_test$p.value
#   
#   # Add results to the dataframe
#   ba_correlations <- rbind(ba_correlations, data.frame(
#     column = i,
#     estimate = estimate,
#     p_value = p_value,
#     stringsAsFactors = FALSE
#   ))
# }
# 
# clvsd <- merge(cl, tcga, by="GENE")
# 
# # Initialize a dataframe to store results
# cl_correlations <- data.frame(
#   column = character(),
#   estimate = numeric(),
#   p_value = numeric(),
#   stringsAsFactors = FALSE
# )
# 
# rownames(clvsd) <- clvsd$GENE
# clvsd$GENE <- NULL
# 
# 
# for (i in colnames(clvsd)[2:567]) {
#   # Perform correlation test with column 21
#   j <- clvsd[,i]
#   cor_test <- cor.test(j, clvsd$CL, method = "pearson")
#   
#   # Extract estimate and p-value
#   estimate <- cor_test$estimate
#   p_value <- cor_test$p.value
#   
#   # Add results to the dataframe
#   cl_correlations <- rbind(cl_correlations, data.frame(
#     column = i,
#     estimate = estimate,
#     p_value = p_value,
#     stringsAsFactors = FALSE
#   ))
# }
# 
# names(ba_correlations)[names(ba_correlations) =="estimate"] <- "cor_BA"
# names(cl_correlations)[names(cl_correlations) =="estimate"] <- "cor_CL"
# names(ims_correlations)[names(ims_correlations) =="estimate"] <- "cor_IMS"
# 
# 
# ba_correlations[ ba_correlations$p_value > 0.05, 'cor_BA'] <- -1
# cl_correlations[ cl_correlations$p_value > 0.05, 'cor_CL'] <- -1
# ims_correlations[ ims_correlations$p_value > 0.05, 'cor_IMS'] <- -1
# 
# res <- merge(ba_correlations, cl_correlations, by="column")
# res$p_value.x <- NULL
# res$p_value.y <- NULL
# res <- merge(res, ims_correlations, by="column")
# res$p_value <- NULL
# 
# 
# dp <- res
# dp$column <- NULL
# pheatmap(dp)
# 
# table(colnames(dp)[apply(dp, 1, which.max)])


## prendere i PDX e vedere se sono confrontabili con gli umani

vsd_xeno <- "/mnt/cold1/snaketree/prj/hn/dataset/V1/RNAseq_likeMK/PDX/fpkm.tsv.gz"
vsd_xeno <- read.table(vsd_xeno, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
vsd_xeno <- log2(vsd_xeno+1)

vsd_xeno2 <- as.data.frame(t(vsd_xeno))
vsd_xeno2$model <- substr(rownames(vsd_xeno2), 1, 7)
vsd_xeno2 <- vsd_xeno2 %>%
  group_by(model) %>%
  summarise(across(everything(), ~ mean(.x, na.rm = TRUE)))
vsd_xeno2 <- as.data.frame(vsd_xeno2)
rownames(vsd_xeno2) <- vsd_xeno2$model
vsd_xeno2$model <- NULL
centered_data_xeno <- scale(vsd_xeno2, scale = FALSE)
vsd_xeno3 <- as.data.frame(t(centered_data_xeno))

centroids <- "/mnt/cold1/snaketree/prj/hn/local/share/data/centroids_hn.xlsx"
centroids <- read.xlsx(centroids, colNames = TRUE)
colnames(centroids) <- centroids[1,]
centroids <- centroids[-1,]
centroids$IMS <- as.numeric(centroids$IMS)
centroids$BA <- as.numeric(centroids$BA)
centroids$CL <- as.numeric(centroids$CL)
ims <- centroids
ims$BA <- NULL
ims$CL <- NULL
ba <- centroids
ba$IMS <- NULL
ba$CL <- NULL
cl <- centroids
cl$IMS <- NULL
cl$BA <- NULL

vsd_xeno3 <- as.data.frame(vsd_xeno3)
vsd_xeno3$GENE <- rownames(vsd_xeno3)

imsvsd_xeno <- merge(ims, vsd_xeno3, by="GENE")

# Initialize a dataframe to store results
ims_xeno_correlations <- data.frame(
  column = character(),
  estimate = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

rownames(imsvsd_xeno) <- imsvsd_xeno$GENE
imsvsd_xeno$GENE <- NULL

for (i in colnames(imsvsd_xeno)[2:53]) {
  j <- imsvsd_xeno[,i]
  cor_test <- cor.test(j, imsvsd_xeno$IMS, method = "pearson")
  estimate <- cor_test$estimate
  p_value <- cor_test$p.value
  ims_xeno_correlations <- rbind(ims_xeno_correlations, data.frame(
    column = i,
    estimate = estimate,
    p_value = p_value,
    stringsAsFactors = FALSE
  ))
}

bavsd_xeno <- merge(ba, vsd_xeno3, by="GENE")

# Initialize a dataframe to store results
ba_xeno_correlations <- data.frame(
  column = character(),
  estimate = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

rownames(bavsd_xeno) <- bavsd_xeno$GENE
bavsd_xeno$GENE <- NULL

for (i in colnames(bavsd_xeno)[2:53]) {
  j <- bavsd_xeno[,i]
  cor_test <- cor.test(j, bavsd_xeno$BA, method = "pearson")
  estimate <- cor_test$estimate
  p_value <- cor_test$p.value
  ba_xeno_correlations <- rbind(ba_xeno_correlations, data.frame(
    column = i,
    estimate = estimate,
    p_value = p_value,
    stringsAsFactors = FALSE
  ))
}

clvsd_xeno <- merge(cl, vsd_xeno3, by="GENE")

# Initialize a dataframe to store results
cl_xeno_correlations <- data.frame(
  column = character(),
  estimate = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

rownames(clvsd_xeno) <- clvsd_xeno$GENE
clvsd_xeno$GENE <- NULL

for (i in colnames(clvsd_xeno)[2:53]) {
  j <- clvsd_xeno[,i]
  cor_test <- cor.test(j, clvsd_xeno$CL, method = "pearson")
  estimate <- cor_test$estimate
  p_value <- cor_test$p.value
  cl_xeno_correlations <- rbind(cl_xeno_correlations, data.frame(
    column = i,
    estimate = estimate,
    p_value = p_value,
    stringsAsFactors = FALSE
  ))
}

names(ba_xeno_correlations)[names(ba_xeno_correlations) =="estimate"] <- "cor_BA"
names(cl_xeno_correlations)[names(cl_xeno_correlations) =="estimate"] <- "cor_CL"
names(ims_xeno_correlations)[names(ims_xeno_correlations) =="estimate"] <- "cor_IMS"


ba_xeno_correlations[ ba_xeno_correlations$p_value > 0.05, 'cor_BA'] <- -1
cl_xeno_correlations[ cl_xeno_correlations$p_value > 0.05, 'cor_CL'] <- -1
ims_xeno_correlations[ ims_xeno_correlations$p_value > 0.05, 'cor_IMS'] <- -1

res_xeno <- merge(ba_xeno_correlations, cl_xeno_correlations, by="column")
res_xeno$p_value.x <- NULL
res_xeno$p_value.y <- NULL
res_xeno <- merge(res_xeno, ims_xeno_correlations, by="column")
res_xeno$p_value <- NULL

dp_xeno <- res_xeno
rownames(dp_xeno) <- dp_xeno$column
dp_xeno$column <- NULL
pheatmap(dp_xeno)

table(colnames(dp_xeno)[apply(dp_xeno, 1, which.max)])

max_col_xeno <- apply(dp_xeno, 1, function(x) colnames(dp_xeno)[which.max(x)])

assign_df_xeno <- data.frame(
  model = rownames(dp_xeno),
  assigned_cor = max_col_xeno,
  stringsAsFactors = FALSE
)

write.table(assign_df_xeno, file="/mnt/cold1/snaketree/prj/hn/local/share/data/classificazione_xeno.tsv",
            quote = FALSE, sep="\t", col.names = TRUE, row.names = FALSE)

assign_df$model <- substr(assign_df$model, 1, 7)
names(assign_df)[names(assign_df)=="assigned_cor"] <- "Class_Human"
assign_df_xeno$model <- substr(assign_df_xeno$model, 1, 7)
names(assign_df_xeno)[names(assign_df_xeno)=="assigned_cor"] <- "Class_Xeno"

controllo <- merge(assign_df, assign_df_xeno, by="model")

controllo$same_diff <- ifelse(
  controllo$Class_Human == controllo$Class_Xeno,
  "same",
  "different"
)


ggplot(controllo, aes(x = Class_Human, y = Class_Xeno, color = same_diff)) +
  geom_jitter(width = 0.2, height = 0.2, size = 3, alpha = 0.7) +
  scale_color_manual(values = c("same" = "steelblue", "different" = "tomato")) +
  theme_minimal() +
  labs(x = "Class Human", y = "Class Xeno", color = "Concordanza")

controllo <- as.data.frame(controllo)
write.table(controllo, file="/mnt/cold1/snaketree/prj/hn/dataset/classificazione_human_xeno.tsv",
            quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)

tab <- as.data.frame(table(controllo$Class_Human, controllo$Class_Xeno))
colnames(tab) <- c("Class_Human", "Class_Xeno", "n")

ggplot(tab, aes(x = Class_Human, y = Class_Xeno, fill = n)) +
  geom_tile(color = "white") +
  geom_text(aes(label = n), color = "black") +
  scale_fill_gradient(low = "white", high = "steelblue") +
  theme_minimal() +
  labs(x = "Class Human", y = "Class Xeno", fill = "Frequenza")

controllo_long <- controllo %>%
  pivot_longer(cols = c(Class_Human, Class_Xeno),
               names_to = "class_type", values_to = "class")

ggplot(controllo_long,
       aes(x = class_type, stratum = class,
           alluvium = model, y = 1, fill = same_diff)) +
  geom_flow(stat = "alluvium", lode.guidance = "forward", color = "grey40") +
  geom_stratum(width = 0.2, fill = "grey80", color = "black") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
  scale_fill_manual(values = c("same" = "steelblue", "different" = "tomato")) +
  theme_minimal() +
  labs(title = "Transizione Classi Human → Xeno",
       x = "", y = "Numero modelli")
ggplot(controllo, aes(x = Class_Human, y = model, fill = same_diff)) +
  geom_tile(color = "white") +
  geom_text(aes(label = Class_Xeno), color = "black", size = 3) +
  scale_fill_manual(values = c("same" = "steelblue", "different" = "tomato")) +
  theme_minimal() +
  labs(title = "Classi Human vs Xeno per modello",
       x = "Class Human", y = "Model")

df <- controllo
df$same_diff <- NULL
rownames(df) <- df$model
df$model <- NULL
df[,1] <- as.factor(df[,1])
df[,2] <- as.factor(df[,2])

links <-
  df %>% 
  mutate(row = row_number()) %>%  # add a row id
  gather('col', 'source', -row) %>%  # gather all columns
  mutate(col = match(col, names(df))) %>%  # convert col names to col nums
  mutate(source = paste0(source, '_', col)) %>%  # add col num to node names
  group_by(row) %>%
  arrange(col) %>%
  mutate(target = lead(source)) %>%  # get target from following node in row
  ungroup() %>% 
  filter(!is.na(target)) %>%  # remove links from last column in original data
  select(source, target) %>% 
  group_by(source, target) %>% 
  summarise(value = n())

# create nodes data frame from unque nodes found in links data frame
nodes <- data.frame(id = unique(c(links$source, links$target)),
                    stringsAsFactors = TRUE)

# remove column id from names
nodes$name <- as.factor(sub('_[0-9]*$', '', nodes$id))

# set links data to the 0-based index of the nodes in the nodes data frame
links$source <- match(links$source, nodes$id) - 1
links$target <- match(links$target, nodes$id) - 1

sn <- sankeyNetwork(Links = links, Nodes = nodes, Source = 'source',
                    Target = 'target', Value = 'value', NodeID = 'name',
                    fontSize = 12, nodeWidth = 30, iterations = 0)

saveWidget(sn, "/home/mferri/sankey_hn.html")


## jaccard
jaccard_overall <- sum(df$Class_Human == df$Class_Xeno) / nrow(df)

subtypes <- c("cor_BA", "cor_CL", "cor_IMS")

jaccard_by_subtype <- sapply(subtypes, function(sub) {
  idx <- which(df$Class_Human == sub | df$Class_Xeno == sub)
  intersect <- sum(df$Class_Human[idx] == sub & df$Class_Xeno[idx] == sub)
  union <- length(unique(c(
    which(df$Class_Human[idx] == sub),
    which(df$Class_Xeno[idx] == sub)
  )))
  return(intersect / union)
})

jaccard_overall
jaccard_by_subtype

## arricchimenti con la risposta
risposta <- "/mnt/cold1/snaketree/prj/hn/local/share/data/def_cohort/tot_resp"
risposta <- read.table(risposta, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
risposta <- risposta[risposta$Resp != "", ]
names(risposta)[names(risposta)=="Case"] <- "model"

merged <- merge(assign_df_xeno, risposta, by="model")
merged$Resp <- ifelse(merged$Resp == "Responder", "Resp", "Non_resp")
tab_xeno <- table(merged$Class_Xeno, merged$Resp)

fisher_xeno <- fisher.test(tab_xeno)
fisher_xeno

df2 <- merge(assign_df_xeno, risposta, by="model")
df2$Class_group <- ifelse(df2$Class_Xeno == "cor_BA", "cor_BA", "else")

cor_ba_tab <- table(df2$Class_group, df2$Resp)
fisher_ba <- fisher.test(cor_ba_tab)
fisher_ba

df3 <- merge(assign_df_xeno, risposta, by="model")
df3$Class_group <- ifelse(df3$Class_Xeno == "cor_CL", "cor_CL", "else")

cor_cl_tab <- table(df3$Class_group, df3$Resp)
fisher_cl <- fisher.test(cor_cl_tab)
fisher_cl

df4 <- merge(assign_df_xeno, risposta, by="model")
df4$Class_group <- ifelse(df4$Class_Xeno == "cor_IMS", "cor_IMS", "else")

cor_ims_tab <- table(df4$Class_group, df4$Resp)
fisher_ims <- fisher.test(cor_ims_tab)
fisher_ims

## calcolo percentuali
table(controllo$Class_Human)

ims_same <- sum(controllo$Class_Human == "cor_IMS" & controllo$Class_Xeno == "cor_IMS")
ims_same
ims_tot <- sum(controllo$Class_Human == "cor_IMS")
ims_tot
ims_perc <- round(ims_same / ims_tot * 100, 1)
ims_perc

ims_reas <- ims_tot - ims_same
ims_reas

ba_reas <- sum(controllo$Class_Human == "cor_IMS" & controllo$Class_Xeno == "cor_BA")
ba_reas
cl_reas <- sum(controllo$Class_Human == "cor_IMS" & controllo$Class_Xeno == "cor_CL")
cl_reas

cl_same <- sum(controllo$Class_Human == "cor_CL" & controllo$Class_Xeno == "cor_CL")
cl_same
cl_tot <- sum(controllo$Class_Human == "cor_CL")
cl_tot
cl_perc <- round(cl_same / cl_tot * 100, 1)
cl_perc

ba_same <- sum(controllo$Class_Human == "cor_BA" & controllo$Class_Xeno == "cor_BA")
ba_same
ba_tot <- sum(controllo$Class_Human == "cor_BA")
ba_tot
ba_perc <- round(ba_same / ba_tot * 100, 1)
ba_perc

## barplot con risposta 
m <- merge(risposta, assign_df_xeno, by="model")
table(m$Class_Xeno)
table(m$Resp)
m$Class_Xeno <- gsub("cor_", "", m$Class_Xeno)

pdf <- as.data.frame(table(m$Class_Xeno, m$Resp))

ggplot(data=pdf, aes(x=Var1, y=Freq, fill=Var2))+geom_col(position='dodge')

write.table(m, file="/mnt/cold1/snaketree/prj/hn/dataset/classificazione_xeno_risposta.tsv", quote = FALSE, sep = "\t",
            col.names = TRUE, row.names = FALSE)

pdf_plot <- m %>%
  dplyr::count(Class_Xeno, Resp, name = "Freq") %>%
  group_by(Class_Xeno) %>%
  mutate(Perc = Freq / sum(Freq) * 100) %>%
  ungroup() %>%
  dplyr::rename(Var1 = Class_Xeno, Var2 = Resp)

ggplot(pdf_plot, aes(x = Var1, y = Perc, fill = Var2)) +
  geom_col(width = 0.7) +
  geom_text(aes(label = ifelse(Freq > 0, paste0(round(Perc, 1), "%\n(n=", Freq, ")"), "")),
            position = position_stack(vjust = 0.5), size = 3.5) +
  scale_y_continuous(labels = function(x) paste0(x, "%"),
                     expand = expansion(mult = c(0, 0.02))) +
  scale_fill_manual(values = c("Responder" = "#2C7BB6", "Non Responder" = "#D7191C")) +
  labs(x = "Class_Xeno", y = "Percentuale di modelli", fill = NULL) +
  theme_classic(base_size = 14)

pdf_plot <- pdf %>%
  group_by(Var2) %>%
  mutate(Perc = Freq / sum(Freq) * 100) %>%
  ungroup()

pdf_plot$Var2 <- gsub("Non Responder", "Non-responders", pdf_plot$Var2)
pdf_plot$Var2 <- gsub("Responder", "Responders", pdf_plot$Var2)

ggplot(pdf_plot, aes(x = Var1, y = Perc, fill = Var2)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = paste0(round(Perc, 1), "%\n(n=", Freq, ")")),
            position = position_dodge(width = 0.8), vjust = -0.3, size = 3.5) +
  scale_y_continuous(labels = function(x) paste0(x, "%"),
                     expand = expansion(mult = c(0, 0.15))) +
  scale_fill_manual(values = c("Responders" = "#1d7937", "Non-responders" = "#772b84"),
                    breaks = c("Responders", "Non-responders")) +
  labs(x = NULL, y = "% of group total", fill = NULL) +
  theme_classic(base_size = 14)
