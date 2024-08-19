w <- "/home/mferri/hnc_new_data_conditi.tsv"
w <- read.table(w, quote = "", header = TRUE, sep = "\t", stringsAsFactors = FALSE)
# w <- w[,c(1,9)]
# w <- w %>% filter(!Average.3.6WKS == "")
casi <- w$Genealogy
# w$Average.3.6WKS <- gsub("%", "", w$Average.3.6WKS)
# w$Average.3.6WKS <- as.numeric(w$Average.3.6WKS)
w <- w[,c(1,7, 10)]
w <- w %>% filter(seq == "YES")
w <- w[order(w$Average.3.6WKS, decreasing = TRUE),]
w$seq <- NULL

write.xlsx(w, file = "/home/mferri/oncoprint_hn.xlsx")

xeno <- "/mnt/cold1/snaketree/prj/hn/dataset/V1/WES/xeno_gene_matrix.tsv"
xeno <- read.table(xeno, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
geni <- c("TP53", "NOTCH1", "KMT2D", "CDKN2A", "FAT1", "HRAS", "KRAS", "NRAS")
xeno <- xeno %>% filter(Hugo_Symbol %in% geni)
rownames(xeno) <- xeno$Hugo_Symbol
xeno$Hugo_Symbol <- NULL
xeno <- as.data.frame(t(xeno))

xeno <- as.data.frame(apply(xeno, 2, function(x) ifelse(x != 0, 1, 0)))

xeno$cases <- substr(rownames(xeno), 1, 7)
xeno <- xeno %>% filter(cases %in% casi)
rownames(xeno) <- xeno$cases
xeno$cases <- NULL
xeno$HRAS <- 0
xeno$NRAS <- 0

mat_list <- list(HR=as.matrix(t(xeno)))

col = c("wt"="#969595", "Mut"= "firebrick")
alter_fun = list(
  wt = function(x, y, w, h) {
    grid.rect(x, y, w*0.5, h, 
              gp = gpar(fill = col["wt"], col = NA))
  },
  Mut = function(x, y, w, h) {
    grid.rect(x, y, w*0.5, h, 
              gp = gpar(fill = col["Mut"], col = NA))
  }
)

# to order genes depending on total muts
s <- t(mat_list[[1]])
su <- colSums(s, na.rm = TRUE)
su <- su[order(-su)]

mat_list <-mat_list[[1]]
mat_list[is.na(mat_list)] <- 3
mat_list <- ifelse(mat_list==1, "Mut", "wt")

ordine_casi <- w$Genealogy
ordered_mat_list <- mat_list[, ordine_casi]
print(length(ordine_casi))

#oncoPrint(mat_list2, alter_fun = alter_fun, col = col)
op <- oncoPrint(ordered_mat_list, alter_fun = alter_fun, col = col, row_order = names(su),
                remove_empty_columns = FALSE, remove_empty_rows = FALSE, column_names_gp=gpar(fontsize=8),
                top_annotation = HeatmapAnnotation(cbar = anno_oncoprint_barplot(), annotation_name_gp=gpar(fontsize=8)), show_pct=FALSE)

get_recist <- function(x) {
  res <- vector(mode="character", length=length(x))
  #res <- ifelse(x < -50, 'OR', ifelse(x > 35, 'PD', 'SD'))
  res <- ifelse(x < -50, "springgreen3", ifelse(x > 35, "firebrick2", "darkgoldenrod3"))
  #res[is.na(x)] <- 'black'
  return(res)
}
op2 <- oncoPrint(ordered_mat_list, alter_fun = alter_fun, col = col, row_order = names(su), column_order=w$Case,
                 remove_empty_columns = FALSE, remove_empty_rows = FALSE, column_names_gp=gpar(fontsize=8),
                 top_annotation = HeatmapAnnotation(cbar = anno_oncoprint_barplot(),
                                                    Cetuximab = anno_barplot(w$Average.3.6WKS,  gp = gpar(fill = get_recist(w$Average.3.6WKS))),
                                                    annotation_name_gp=gpar(fontsize=8), height = unit(4, "cm")), show_pct=FALSE)

pdf("/mnt/cold1/snaketree/prj/hn/dataset/V1/WES/oncoprint_waterfall.pdf")
op2
dev.off()

## fisher

colnames(w) <- c('smodel', 'perc')

data_merged <- merge(xeno, w, by.x="row.names", by.y='smodel')
# if (snakemake@wildcards[['week']] == "w3") {
#   stopifnot(nrow(data_merged)==wanted)
# }

#save.image('pi.Rdata')

fisher <- function(gene, mydata) {
  mydata$class <- ifelse(mydata$perc > 35, "Non_Responder", "Responder")
  ct <- table(mydata$class, mydata[,gene] != 0)
  if (ncol(ct) == 2) {
    ft <- fisher.test(ct)
    return(ft$p.value)
  } else {
    return(NA)
  }
}

mydata <- data_merged
mydata$class <- ifelse(mydata$perc > 35, "Non_Responder", "Responder")

fish <- as.data.frame(sapply(colnames(xeno), fisher, data_merged))
colnames(fish) <- c('pvalue')
fish <- fish[!is.na(fish$pvalue),, drop=FALSE]
fish$padj <- p.adjust(fish$pvalue, method="BH")
fish <- fish[order(fish$pvalue),]

waterfall <- read_excel("/mnt/cold1/snaketree/prj/hn/local/share/data/hnc_response_waterfall.xlsx")
waterfall <- as.data.frame(waterfall[c(1:80), c("Genealogy", "Average 3-6WKS (log)", "Average 3-6WKS response class")])
colnames(waterfall) <- c("Genealogy", "logAverage3_6", "response_class")

p <- ggplot(waterfall, aes(x = reorder(Genealogy, -logAverage3_6), y = logAverage3_6, fill=response_class)) + geom_bar(stat = "identity")+
  theme(axis.text.x = element_text(angle = 90, size = 5))+scale_fill_manual(values = c("firebrick2", "springgreen4", "darkgoldenrod3"))

pdf("/mnt/cold1/snaketree/prj/hn/dataset/V1/WES/waterfall_all.pdf")
p
dev.off()
