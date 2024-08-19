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

write.table(meta, "/mnt/cold1/snaketree/prj/hn/local/share/data/SampleTableHuman.txt", quote = FALSE, sep = "\t", col.names = TRUE, row.names = FALSE)
