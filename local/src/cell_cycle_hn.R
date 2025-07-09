## head and neck cellular cycle

signatures <- c("HALLMARK_E2F_TARGETS",
                "HALLMARK_G2M_CHECKPOINT",
                "HALLMARK_MITOTIC_SPINDLE",
                "HALLMARK_MYC_TARGETS_V1",
                "HALLMARK_MYC_TARGETS_V2",
                "KIM_MYC_AMPLIFICATION_TARGETS_DN")

h <- "/mnt/cold1/snaketree/prj/hn/dataset/V1/RNAseq/xeno/GSEA_H_results_Cetuximab_Response_cutoff0.05-NonResponder.vs.Responder.tsv"
h <- read.table(h, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
h <- h %>% filter(ID %in% signatures)

c2 <- "/mnt/cold1/snaketree/prj/hn/dataset/V1/RNAseq/xeno/GSEA_C2_results_Cetuximab_Response_cutoff0.05-NonResponder.vs.Responder.tsv"
c2 <- read.table(c2, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
c2 <- c2 %>% filter(ID %in% signatures)

e2f <- c("AK2","ANP32E","ASF1A","ASF1B","ATAD2","AURKA","AURKB","BARD1","BIRC5","BRCA1","BRCA2","BRMS1L","BUB1B","CBX5",
         "CCNB2","CCNE1","CCP110","CDC20","CDC25A","CDC25B","CDCA3","CDCA8","CDK1","CDK4","CDKN1A","CDKN1B","CDKN2A","CDKN2C",
         "CDKN3","CENPE","CENPM","CHEK1","CHEK2","CIT","CKS1B","CKS2","CNOT9","CSE1L","CTCF","CTPS1","DCK","DCLRE1B","DCTPP1","DDX39A",
         "DEK","DEPDC1","DIAPH3","DLGAP5","DNMT1","DONSON","DSCC1","DUT","E2F8","EED","EIF2S1","ESPL1","EXOSC8","EZH2","GINS1","GINS3",
         "GINS4","GSPT1","H2AX","H2AZ1","HELLS","HMGA1","HMGB2","HMGB3","HMMR","HNRNPD","HUS1","ILF3","ING3","IPO7","JPT1","KIF18B","KIF22",
         "KIF2C","KIF4A","KPNA2","LBR","LIG1","LMNB1","LUC7L3","LYAR","MAD2L1","MCM2","MCM3","MCM4","MCM5","MCM6","MCM7","MELK","MKI67","MLH1",
         "MMS22L","MRE11","MSH2","MTHFD2","MXD3","MYBL2","MYC","NAA38","NAP1L1","NASP","NBN","NCAPD2","NME1","NOLC1","NOP56","NUDT21","NUP107","NUP153",
         "NUP205","ORC2","ORC6","PA2G4","PAICS","PAN2","PCNA","PDS5B","PHF5A","PLK1","PLK4","PMS2","PNN","POLA2","POLD1","POLD2","POLD3","POLE","POLE4","POP7",
         "PPM1D","PPP1R8","PRDX4","PRIM2","PRKDC","PRPS1","PSIP1","PSMC3IP","PTTG1","RACGAP1","RAD1","RAD21","RAD50","RAD51AP1","RAD51C","RAN","RANBP1","RBBP7",
         "RFC1","RFC2","RFC3","RNASEH2A","RPA1","RPA2","RPA3","RRM2","SHMT1","SLBP","SMC1A","SMC3","SMC4","SMC6","SNRPB","SPAG5","SPC24","SPC25","SRSF1","SRSF2",
         "SSRP1","STAG1","STMN1","SUV39H1","SYNCRIP","TACC3","TBRG4","TCF19","TFRC","TIMELESS","TIPIN","TK1","TMPO","TOP2A","TP53",
         "TRA2B","TRIP13","TUBB","TUBG1","UBE2S","UBE2T","UBR7","UNG","USP1","WDR90","WEE1","XPO1","XRCC6","ZW10")

g2m <- c("ABL1","AMD1","ARID4A","ATF5","ATRX","AURKA","AURKB","BARD1","BCL3","BIRC5","BRCA2","BUB1","BUB3","CASP8AP2","CBX1",
         "CCNA2","CCNB2","CCND1","CCNF","CCNT1","CDC20","CDC25A","CDC25B","CDC27","CDC45","CDC6","CDC7","CDK1","CDK4","CDKN1B",
         "CDKN2C","CDKN3","CENPA","CENPE","CENPF","CHAF1A","CHEK1","CHMP1A","CKS1B","CKS2","CTCF","CUL1","CUL3","CUL4A","CUL5",
         "DBF4","DDX39A","DKC1","DMD","DR1","DTYMK","E2F1","E2F2","E2F3","E2F4","EFNA5","EGF","ESPL1","EWSR1","EXO1","EZH2",
         "FANCC","FBXO5","FOXN3","G3BP1","GINS2","GSPT1","H2AX","H2AZ1","H2AZ2","H2BC12","HIF1A","HIRA","HMGA1","HMGB3","HMGN2",
         "HMMR","HNRNPD","HNRNPU","HOXC10","HSPA8","HUS1","ILF3","INCENP","JPT1","KATNA1","KIF11","KIF15","KIF20B","KIF22","KIF23",
         "KIF2C","KIF4A","KIF5B","KMT5A","KNL1","KPNA2","KPNB1","LBR","LIG3","LMNB1","MAD2L1","MAP3K20","MAPK14","MARCKS","MCM2",
         "MCM3","MCM5","MCM6","MEIS1","MEIS2","MKI67","MNAT1","MT2A","MTF2","MYBL2","MYC","NASP","NCL","NDC80","NEK2","NOLC1","NOTCH2",
         "NSD2","NUMA1","NUP50","NUP98","NUSAP1","ODC1","ODF2","ORC5","ORC6","PAFAH1B1","PBK","PDS5B","PLK1","PLK4","PML","POLA2",
         "POLE","POLQ","PRC1","PRIM2","PRMT5","PRPF4B","PTTG1","PTTG3P","PURA","RACGAP1","RAD21","RAD23B","RAD54L","RASAL2","RBL1",
         "RBM14","RPA2","RPS6KA5","SAP30","SFPQ","SLC12A2","SLC38A1","SLC7A1","SLC7A5","SMAD3","SMARCC1","SMC1A","SMC2","SMC4","SNRPD1",
         "SQLE","SRSF1","SRSF10","SRSF2","SS18","STAG1","STIL","STMN1","SUV39H1","SYNCRIP","TACC3","TENT4A","TFDP1","TGFB1","TLE3",
         "TMPO","TNPO2","TOP1","TOP2A","TPX2","TRA2B","TRAIP","TROAP","TTK","UBE2C","UBE2S","UCK2","UPF1","WRN","XPO1","YTHDC1")

genes <- unique(c(e2f, g2m))

vsd <- "/mnt/cold1/snaketree/prj/hn/dataset/V1/RNAseq/xeno/vsd.tsv.gz"
vsd <- read.table(vsd, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
vsd$geni <- rownames(vsd)
vsd <- vsd %>% filter(geni %in% genes)
vsd$geni <- NULL
vsd <- as.data.frame(t(vsd))

samples <- "/mnt/cold1/snaketree/prj/hn/local/share/data/SampleTablenoCRLF.txt"
samples <- read.table(samples, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
samples$Count_Filename <- NULL
samples$Site <- NULL
rownames(samples) <- samples$SampleID
samples <- samples[order(samples$Cetuximab_Response),]
samples$SampleID <- NULL

vsd <- vsd[rownames(samples), ]
vsdt <- as.data.frame(t(vsd))

pheatmap(vsd, cluster_cols = FALSE, cluster_rows = FALSE, annotation_row = samples)

pheatmap(vsdt, cluster_cols = FALSE, cluster_rows = FALSE, annotation_col = samples)

lfc <- "/mnt/cold1/snaketree/prj/hn/dataset/V1/RNAseq/xeno/Cetuximab_Response_cutoff0.05-NonResponder.vs.Responder.deseq2.tsv"
lfc <- read.table(lfc, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)
lfc$geni <- rownames(lfc)
lfc <- lfc %>% filter(geni %in% genes)

lfc <- lfc[,c("log2FoldChange", "geni")]
lfc$geni <- NULL

pheatmap(lfc, cluster_rows = FALSE, cluster_cols = FALSE)

d <- lfc

minv <- -1.5
maxv <- 1.5
#d[d < -4] <- -4
#d[d > 4] <- 4

neutral_value <- 0
#bk1 <- c(seq(minv-0.1,neutral_value-0.1,by=0.2),neutral_value-0.0999)
bk1 <- seq(minv-0.001, neutral_value-0.0009, length.out=224)
#bk2 <- c(neutral_value+0.001, seq(neutral_value+0.1,maxv+0.1,by=0.2))
bk2 <- seq(neutral_value+0.0001, maxv+0.001, length.out=224)
bk <- c(bk1, bk2)
my_palette <- c(colorRampPalette(colors = c("darkblue",
                                            "lightblue"))(n = length(bk1)-1),
                "#e1e1e1", "#e1e1e1",
                c(colorRampPalette(colors = c("tomato1", "darkred"))(n
                                                                     = length(bk2)-1)))
#pheatmap(matrix, breaks = seq(-rg, rg, length.out = 100))
pheatmap(d, cluster_rows = F, cluster_cols=F,
         breaks = bk, color=my_palette)

lfc <- lfc %>% filter(abs(log2FoldChange) > 1)
