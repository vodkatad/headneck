maryk <- read.table('/mnt/cold1/snaketree/prj/hn/local/share/data/HNSCC_All_PDX_Data.tsv', sep="\t", header=TRUE, stringsAsFactors = FALSE)
eugy <- read.table('/home/mferri/hnc_new_data_conditi.tsv', sep="\t", header=TRUE, stringsAsFactors = FALSE)

noteugy <- setdiff(maryk$Case, eugy$Genealogy)
mah <- maryk[maryk$Case %in% noteugy,]
mmmi <- maryk[maryk$Case %in% noteugy & maryk$Average.3.6WKS.response.class != "", ]

notmary <- setdiff(eugy$Genealogy, maryk$Case)
mah2 <- eugy[eugy$Genealogy %in% notmary,]

length(intersect(maryk$Case, eugy$Genealogy))

m <- merge(maryk, eugy, by.x="Case", by.y="Genealogy")
assert(nrow(m)==75)

table(m$Average.3.6WKS.response.class.x, m$Average.3.6WKS.response.class.y)

all(abs(m$Aver.VolVar3WKS.....y - m$Aver.VolVar3WKS.....x) < 0.0000001)
de <- abs(m$Aver.VolVar6WKS.....y - m$Aver.VolVar6WKS.....x)
all(de[!is.na(de)] < 0.0000001)
all(abs(m$Average.3.6WKS.response.class.y - m$Average.3.6WKS.response.class.y) < 0.0000001)

### 3 - 6
table(eugy[is.na(eugy$Aver.VolVar6WKS....),'Average.3.6WKS.response.class'])

mah3 <- eugy[is.na(eugy$Aver.VolVar6WKS....) & eugy$Average.3.6WKS.response.class=="SD",]

table(maryk[is.na(maryk$Aver.VolVar6WKS....),'Average.3.6WKS.response.class'])

