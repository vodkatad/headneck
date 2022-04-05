library(ComplexHeatmap)
library(ggplot2)

osnakemake <- snakemake

load(snakemake@input[['Rimage']])
op_f <- osnakemake@output[['op']]
op_data_f <- osnakemake@output[['op_data']]
pie_f <- osnakemake@output[['pie']]
log_f <- osnakemake@log[['log']]
n_muts_f <- osnakemake@output[['n_muts']]

# TODO keep only top mutated genes, top 100 then union
n_o <- apply(pdobing, 1, sum)
n_x <- apply(pdxbing, 1, sum)
n_o <- n_o[order(n_o)]
n_x <- n_x[order(n_x)]
chosen_genes <- union(names(tail(n_o, n=60)), names(tail(n_x, n=60)))
n_muts_o <- as.data.frame(n_o)
n_muts_x <- as.data.frame(n_x)
n_muts <- merge(n_muts_o, n_muts_x, by="row.names")
n_muts$top100 <- ifelse(n_muts$Row.names %in% chosen_genes, 'yes', 'no')
write.table(n_muts, n_muts_f, sep="\t", quote=FALSE)
pdobing <- pdobing[rownames(pdobing) %in% chosen_genes,]
pdxbing <- pdxbing[rownames(pdxbing) %in% chosen_genes,]

sink(log_f)
print('Never mutated h')
table(apply(pdobing, 2, function(x) {!any(x)}))
print('Never mutated x')
table(apply(pdxbing, 2, function(x) {!any(x)}))
sink()

mergemut <- merge(pdobing, pdxbing, by="row.names")#, all.x=TRUE, all.y=TRUE, fill=FALSE) # no need all the same genes:
#> setdiff(rownames(pdo_df), rownames(xeno_df))
#character(0)
# TODO TRUE?
rownames(mergemut) <- mergemut$Row.names
mergemut$Row.names <- NULL

d2 <- t(mergemut)
d2 <- ifelse(d2, 1, 0)
#d <- read.table(treat_f, sep="\t", header=FALSE, stringsAsFactors = TRUE, row.names=1)


#rimuovere i non mutati
#> sum(rowSums(d2) == 0)


d2pdo <- d2[grepl('.x', rownames(d2)),]
d2pdx <- d2[grepl('.y', rownames(d2)),]

td2x <- t(d2pdx)
td2o <- t(d2pdo)
colnames(td2x) <- substr(colnames(td2x), 0, 7)
colnames(td2o) <- substr(colnames(td2o), 0, 7)

both <- td2o & td2x
both2 <- t(apply(both, 1, as.numeric))
rownames(both2) <- rownames(both)
mat_list <- list(Both=both, Human= td2o-both, Xeno=td2x-both  ) 

#col = c("Both" = "blue", "PDO" = "red", "PDX"= "#0b7015")
col = c("Both" = "darkgoldenrod3", "Human" = "darkblue", "Xeno"= "firebrick1")
alter_fun = list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w*0.5, h, 
              gp = gpar(fill = "#CCCCCC", col = NA))
  },
  Both = function(x, y, w, h) {
    grid.rect(x, y, w*0.5, h, 
              gp = gpar(fill = col["Both"], col = NA))
  },
  Human = function(x, y, w, h) {
    grid.rect(x, y, w*0.5, h, 
              gp = gpar(fill = col["Human"], col = NA))
  },
  Xeno = function(x, y, w, h) {
    grid.rect(x, y, w*0.5, h, 
              gp = gpar(fill = col["Xeno"], col = NA))
  }
)

s <- mat_list[[1]]+mat_list[[2]]+mat_list[[3]]
su <- colSums(s)
su <- su[order(-su)]

op <- oncoPrint(mat_list, alter_fun = alter_fun, col = col, column_order = names(su),
                remove_empty_columns = TRUE, remove_empty_rows = TRUE, pct_gp=gpar(fontsize=3), column_names_gp=gpar(fontsize=3),
                top_annotation = HeatmapAnnotation(cbar = anno_oncoprint_barplot(), annotation_name_gp=gpar(fontsize=3)))

pdf(op_f, width=4.7, height=4.7, family="sans")
print(op)
graphics.off()

save.image(op_data_f)

pd <- as.data.frame(sapply(mat_list, sum))
colnames(pd) <- "alterations"
pd$class <- rownames(pd)
ggplot(pd, aes(x="", y=alterations, fill=class))+ geom_bar(width = 1, stat = "identity")+ coord_polar("y", start=0)+scale_fill_manual(values=c("darkgoldenrod3", "darkblue", "firebrick1"))+theme_bw()
ggsave(pie_f)