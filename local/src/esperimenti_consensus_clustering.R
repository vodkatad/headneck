vsd <- "/mnt/cold1/snaketree/prj/hn/dataset/V1/RNAseq/xeno/vsd.tsv.gz"
vsd <- read.table(vsd, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

library(ALL) 
data(ALL) 
d=exprs(ALL) 
mads=apply(d,1,mad) 
d=d[rev(order(mads))[1:5000],]
d = sweep(d,1, apply(d,1,median,na.rm=T))

library(ConsensusClusterPlus)
title="/home/mferri/provahn"
results = ConsensusClusterPlus(d,
                               maxK=6,
                               reps=50,
                               pItem=0.8,
                               pFeature=1,
                               title=title,
                               clusterAlg="hc",
                               distance="pearson",
                               seed=1262118388.71279,
                               plot = "pdf")
results[[2]][["consensusMatrix"]][1:5,1:5]
results[[2]][["consensusTree"]]

d <- as.matrix(vsd)

results = ConsensusClusterPlus(d, # normalized expression matrix
                               maxK=6, # maximum k to make clusters
                               reps=50, # resampling
                               pItem=0.8, # item resampling
                               pFeature=1, # gene resampling
                               title= title, # output title 
                               clusterAlg="pam", # clustering algorithm
                               distance="spearman", # distance measure
                               seed=1262118388.71279, # random number
                               plot="pdf") # output file format


##  nuovo round

vsd <- "/mnt/cold1/snaketree/prj/hn/dataset/V1/RNAseq/umani/vsd.tsv.gz"
vsd <- read.table(vsd, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
vsd <- log2(vsd)

centered_data <- scale(vsd, scale = FALSE)

title <- "/home/mferri/prova"
results <- ConsensusClusterPlus(
  d = centered_data,
  maxK = 6,
  reps = 100,
  pItem = 0.8,
  pFeature = 1,
  title = title,
  clusterAlg = "hc",
  distance = "euclidean",
  seed = 123,
  plot = "pdf"
)

# Determine the optimal number of clusters (for example, k=3)
optimal_k <- 5

# Extract consensus clustering results for the optimal number of clusters
clusters <- results[[optimal_k]]$consensusClass

# Combine data with clusters for aggregation
data_with_clusters <- data.frame(centered_data, cluster = clusters)

# Compute centroids for each cluster
centroids <- aggregate(. ~ cluster, data = data_with_clusters, FUN = mean)

# Remove the cluster column to get the centroids matrix
centroids_matrix <- as.matrix(centroids[,-1])

# Function to assign samples to clusters based on nearest centroid
assign_to_cluster <- function(sample, centroids_matrix) {
  distances <- apply(centroids_matrix, 1, function(centroid) {
    sum((sample - centroid)^2)
  })
  return(which.min(distances))
}

# Classify each sample
predicted_clusters <- apply(centered_data, 1, function(sample) {
  assign_to_cluster(sample, centroids_matrix)
})

# Add predicted clusters to the data frame
classified_data <- data.frame(centered_data, cluster = predicted_clusters)

# View the classified data with clusters
print(classified_data)


## cosa stupida per ripartire

vsd <- "/mnt/cold1/snaketree/prj/hn/dataset/V1/RNAseq/umani/vsd.tsv.gz"
vsd <- read.table(vsd, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
vsd <- log2(vsd)

centered_data <- scale(vsd, scale = FALSE)

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

vsd <- as.data.frame(centered_data)
vsd$GENE <- rownames(vsd)

imsvsd <- merge(ims, vsd, by="GENE")

# Initialize a dataframe to store results
ims_correlations <- data.frame(
  column = character(),
  estimate = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

rownames(imsvsd) <- imsvsd$GENE
imsvsd$GENE <- NULL

# Loop through columns 1 to 20
for (i in colnames(imsvsd)[2:37]) {
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

bavsd <- merge(ba, vsd, by="GENE")

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
for (i in colnames(bavsd)[2:37]) {
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

clvsd <- merge(cl, vsd, by="GENE")

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
for (i in colnames(clvsd)[2:37]) {
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

res <- merge(ba_correlations, cl_correlations, by="column")
res$p_value.x <- NULL
res$p_value.y <- NULL
res <- merge(res, ims_correlations, by="column")
res$p_value <- NULL
