h <- "/mnt/cold1/snaketree/prj/hn/local/share/data/Merged_raw_counts_human.txt"
h <- read.table(h, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
colnames(h) <- gsub("X.", "", colnames(h))
colnames(h) <- gsub('|.$', '', colnames(h))
h$GENE <- gsub("\\.' ", "", h$GENE)
h$GENE <- gsub('^"|"$', '', h$GENE)

meta <- "/mnt/cold1/snaketree/prj/hn/local/share/data/SampleTablenoCRLF.txt"
meta <- read.table(meta, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
meta$Count_Filename <- NULL
meta$ssample <- substr(meta$SampleID, 1, 7)

umani <- as.data.frame(colnames(h)[2:37])
umani$ssample <- substr(umani$`colnames(h)[2:37]`, 1, 7)

meta <- merge(meta, umani, by="ssample")
meta <- meta[!duplicated(meta$`colnames(h)[2:37]`),]
meta$SampleID <- meta$`colnames(h)[2:37]`
meta$`colnames(h)[2:37]` <- NULL
meta$ssample <- NULL

x <- "/mnt/cold1/snaketree/prj/hn/local/share/data/Merged_raw_counts_noCRLF.txt"
x <- read.table(x, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
colnames(x) <- gsub("X\\.", "", colnames(x))
colnames(x) <- gsub('|.$', '', colnames(x))
x$GENE <- gsub("\\.' ", "", x$GENE)
x$GENE <- gsub('^"|"$', '', x$GENE)

metax <- "/mnt/cold1/snaketree/prj/hn/local/share/data/SampleTablenoCRLF.txt"
metax <- read.table(metax, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
metax$Count_Filename <- NULL
metax$ssample <- metax$SampleID

xeno <- as.data.frame(colnames(x)[2:61])
xeno$ssample <- xeno$`colnames(x)[2:61]`

metax <- merge(metax, xeno, by="ssample")
metax <- metax[!duplicated(metax$`colnames(x)[2:61]`),]
metax$SampleID <- metax$`colnames(x)[2:61]`
metax$`colnames(x)[2:61]` <- NULL
metax$ssample <- NULL

meta_all <- rbind(metax, meta)
meta_all <- meta_all[,c(3,1,2)]
meta_all$type <- substr(meta_all$SampleID, 8, 10)

write.table(meta_all, file="/mnt/cold1/snaketree/prj/hn/local/share/data/SampleTable_hx.txt", quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)

counts <- merge(h, x, by="GENE")

write.table(counts, file="/mnt/cold1/snaketree/prj/hn/local/share/data/Merged_all_counts_hx.txt", quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
