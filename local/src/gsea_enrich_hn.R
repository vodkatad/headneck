gseaScores <- function(geneList, geneSet, exponent=1, fortify=FALSE) {
  ###################################################################
  ##    geneList                                                   ##
  ##                                                               ##
  ## 1. Rank order the N genes in D to form L = { g_1, ... , g_N}  ##
  ##    according to the correlation, r(g_j)=r_j,                  ##
  ##    of their expression profiles with C.                       ##
  ##                                                               ##
  ###################################################################
  
  ###################################################################
  ##    exponent                                                   ##
  ##                                                               ##
  ## An exponent p to control the weight of the step.              ##
  ##   When p = 0, Enrichment Score ( ES(S) ) reduces to           ##
  ##   the standard Kolmogorov-Smirnov statistic.                  ##
  ##   When p = 1, we are weighting the genes in S                 ##
  ##   by their correlation with C normalized                      ##
  ##   by the sum of the correlations over all of the genes in S.  ##
  ##                                                               ##
  ###################################################################
  
  ## genes defined in geneSet should appear in geneList.
  ## this is a must, see https://github.com/GuangchuangYu/DOSE/issues/23
  geneSet <- intersect(geneSet, names(geneList))
  
  N <- length(geneList)
  Nh <- length(geneSet)
  
  Phit <- Pmiss <- numeric(N)
  hits <- names(geneList) %in% geneSet ## logical
  
  Phit[hits] <- abs(geneList[hits])^exponent
  NR <- sum(Phit)
  Phit <- cumsum(Phit/NR)
  
  Pmiss[!hits] <-  1/(N-Nh)
  Pmiss <- cumsum(Pmiss)
  
  runningES <- Phit - Pmiss
  
  ## ES is the maximum deviation from zero of Phit-Pmiss
  max.ES <- max(runningES)
  min.ES <- min(runningES)
  if( abs(max.ES) > abs(min.ES) ) {
    ES <- max.ES
  } else {
    ES <- min.ES
  }
  
  df <- data.frame(x=seq_along(runningES),
                   runningScore=runningES,
                   position=as.integer(hits)
  )
  
  if(fortify==TRUE) {
    return(df)
  }
  
  df$gene = names(geneList)
  res <- list(ES=ES, runningES = df)
  return(res)
}
gsInfo <- function(object, geneSetID) {
  geneList <- object@geneList
  
  if (is.numeric(geneSetID))
    geneSetID <- object@result[geneSetID, "ID"]
  
  geneSet <- object@geneSets[[geneSetID]]
  exponent <- object@params[["exponent"]]
  df <- gseaScores(geneList, geneSet, exponent, fortify=TRUE)
  df$ymin <- 0
  df$ymax <- 0
  pos <- df$position == 1
  h <- diff(range(df$runningScore))/20
  df$ymin[pos] <- -h
  df$ymax[pos] <- h
  df$geneList <- geneList
  
  df$Description <- object@result[geneSetID, "Description"]
  return(df)
}

mygseaplot <- function (x, geneSetID, title = "", color = "green", base_size = 11, 
                        rel_heights = c(1.5, 0.5, 1), subplots = 1:3, pvalue_table = FALSE, 
                        ES_geom = "line") 
{
  ES_geom <- match.arg(ES_geom, c("line", "dot"))
  geneList <- position <- NULL
  if (length(geneSetID) == 1) {
    gsdata <- gsInfo(x, geneSetID)
  }
  else {
    gsdata <- do.call(rbind, lapply(geneSetID, gsInfo, object = x))
  }
  p <- ggplot(gsdata, aes_(x = ~x)) + xlab(NULL) + theme_classic(base_size) + 
    theme(panel.grid.major = element_line(colour = "grey92"), 
          panel.grid.minor = element_line(colour = "grey92"), 
          panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank()) + 
    scale_x_continuous(expand = c(0, 0))
  if (ES_geom == "line") {
    es_layer <- geom_line(aes_(y = ~runningScore, color = ~Description), 
                          size = 1)
  }
  else {
    es_layer <- geom_point(aes_(y = ~runningScore, color = ~Description), 
                           size = 1, data = subset(gsdata, position == 1))
  }
  p.res <- p + es_layer + theme(legend.position = c(0.8, 0.8), 
                                legend.title = element_blank(), legend.background = element_rect(fill = "transparent"))
  p.res <- p.res + ylab("Running Enrichment Score") + theme(axis.text.x = element_blank(), 
                                                            axis.ticks.x = element_blank(), axis.line.x = element_blank(), 
                                                            plot.margin = margin(t = 0.2, r = 0.2, b = 0, l = 0.2, 
                                                                                 unit = "cm"))
  i <- 0
  for (term in unique(gsdata$Description)) {
    idx <- which(gsdata$ymin != 0 & gsdata$Description == 
                   term)
    gsdata[idx, "ymin"] <- i
    gsdata[idx, "ymax"] <- i + 1
    i <- i + 1
  }
  p2 <- ggplot(gsdata, aes_(x = ~x)) + geom_linerange(aes_(ymin = ~ymin, 
                                                           ymax = ~ymax, color = ~Description)) + xlab(NULL) + ylab(NULL) + 
    theme_classic(base_size) + theme(legend.position = "none", 
                                     plot.margin = margin(t = -0.1, b = 0, unit = "cm"), axis.ticks = element_blank(), 
                                     axis.text = element_blank(), axis.line.x = element_blank()) + 
    scale_x_continuous(expand = c(0, 0)) + scale_y_continuous(expand = c(0, 
                                                                         0))
  if (length(geneSetID) == 1) {
    v <- seq(1, sum(gsdata$position), length.out = 9)
    inv <- findInterval(rev(cumsum(gsdata$position)), v)
    if (min(inv) == 0) 
      inv <- inv + 1
    col = c(rev(brewer.pal(5, "Blues")), brewer.pal(5, "Reds"))
    ymin <- min(p2$data$ymin)
    yy <- max(p2$data$ymax - p2$data$ymin) * 0.3
    xmin <- which(!duplicated(inv))
    xmax <- xmin + as.numeric(table(inv)[unique(inv)])
    d <- data.frame(ymin = ymin, ymax = yy, xmin = xmin, 
                    xmax = xmax, col = col[unique(inv)])
    p2 <- p2 + geom_rect(aes_(xmin = ~xmin, xmax = ~xmax, 
                              ymin = ~ymin, ymax = ~ymax, fill = ~I(col)), data = d, 
                         alpha = 0.9, inherit.aes = FALSE)
  }
  df2 <- p$data
  df2$y <- p$data$geneList[df2$x]
  p.pos <- p + geom_segment(data = df2, aes_(x = ~x, xend = ~x, 
                                             y = ~y, yend = 0), color = "grey")
  p.pos <- p.pos + ylab("Ranked list metric") + xlab("Rank in Ordered Dataset") + 
    theme(plot.margin = margin(t = -0.1, r = 0.2, b = 0.2, 
                               l = 0.2, unit = "cm"))
  if (!is.null(title) && !is.na(title) && title != "") 
    p.res <- p.res + ggtitle(title)
  if (length(color) == length(geneSetID)) {
    p.res <- p.res + scale_color_manual(values = color)
    if (length(color) == 1) {
      p.res <- p.res + theme(legend.position = "none")
      p2 <- p2 + scale_color_manual(values = "black")
    }
    else {
      p2 <- p2 + scale_color_manual(values = color)
    }
  }
  if (pvalue_table) {
    pd <- x[geneSetID, c("Description", "pvalue", "p.adjust")]
    pd <- pd[order(pd[, 1], decreasing = FALSE), ]
    rownames(pd) <- pd$Description
    pd <- pd[, -1]
    pd <- round(pd, 4)
    tp <- tableGrob2(pd, p.res)
    p.res <- p.res + theme(legend.position = "none") + annotation_custom(tp, 
                                                                         xmin = quantile(p.res$data$x, 0.5), xmax = quantile(p.res$data$x, 
                                                                                                                             0.95), ymin = quantile(p.res$data$runningScore, 
                                                                                                                                                    0.75), ymax = quantile(p.res$data$runningScore, 
                                                                                                                                                                           0.9))
  }
  p.res <- p.res + theme(legend.position="none", axis.text.y=element_blank(), axis.title.y=element_blank())
  p.pos <- p.pos + theme(axis.text.x=element_blank(), axis.text.y=element_blank())
  p2 <- p2 + theme(axis.text.x=element_blank(), axis.text.y=element_blank())
  plotlist <- list(p.res, p2, p.pos)[subplots]
  n <- length(plotlist)
  plotlist[[n]] <- plotlist[[n]] + theme(axis.line.x = element_line(), 
                                         axis.ticks.x = element_line(), axis.text.x = element_blank())
  if (length(subplots) == 1) 
    return(plotlist[[1]] + theme(plot.margin = margin(t = 0.2, 
                                                      r = 0.2, b = 0.2, l = 0.2, unit = "cm")))
  if (length(rel_heights) > length(subplots)) 
    rel_heights <- rel_heights[subplots]
  
  plot_grid(plotlist = plotlist, ncol = 1, align = "v", rel_heights = rel_heights)
}


load("/mnt/cold1/snaketree/prj/hn/dataset/V1/RNAseq/gsea_results_C5.Rdata")
load(rdata)

tot <- read.table("/mnt/cold1/snaketree/prj/hn/dataset/V1/RNAseq/GSEA_C5_results_Cetuximab_Response_cutoff0.05-NonResponder.vs.Responder.tsv", quote = "",sep = "\t", header = TRUE)
#tot <- read.table(data_f, quote = "", sep = "\t", header = TRUE, stringsAsFactors = FALSE)

blu <- c("REACTOME_REPRESSION_OF_WNT_TARGET_GENES", "PID_HEDGEHOG_2PATHWAY")
blu <- c("HALLMARK_HEDGEHOG_SIGNALING")
blu <- c("GOBP_KERATINIZATION")

blucheck <- tot %>% filter(ID %in% blu)
blucheck <- blucheck %>% filter(pvalue < 0.05)
res <- blucheck[c(1,4,5,6,8)]
wi <- c()

for (w in blu) {
  i <- which(em@result$Description==w)
  wi <- c(wi, i)
  #print(gseaplot(em, geneSetID = i, title = em$Description[i], by = "runningScore"))
}

## cambiare palette sul blu
p <- gseaplot2(em, geneSetID=wi, color = c("#74A9CF", "#3690C0", "#0570B0"), subplots = 1:2)

pdf(p, "/mnt/cold1/snaketree/prj/hn/dataset/V1/RNAseq/plot_gsea_Keratinization.pdf")
p
dev.off()

ggsave(p, file=blu_s, width=200, height=89, units="mm", device="jpeg", dpi=300)
p <- mygseaplot(em, geneSetID=wi, color = c("#74A9CF", "#3690C0", "#0570B0"), subplots = 1:2)
ggsave(p, file=blu_m, width=200, height=89, units="mm", device="jpeg", dpi=300)

bordeaux <- c("HALLMARK_FATTY_ACID_METABOLISM", "REACTOME_FATTY_ACID_METABOLISM")
bordeauxcheck <- tot %>% filter(ID %in% bordeaux)
bordeauxcheck <- bordeauxcheck %>% filter(pvalue < 0.05)
res <- rbind(res, bordeauxcheck[c(1,4,5,6,8)])

wi <- c()

for (w in bordeaux) {
  i <- which(em@result$Description==w)
  wi <- c(wi, i)
  #print(gseaplot(em, geneSetID = i, title = em$Description[i], by = "runningScore"))
}

## cambiare palette sul bordeaux
p <- gseaplot2(em, geneSetID=wi, color = c("#A50026", "#D73027"), subplots = 1:2)
ggsave(p, file=gialli_s, width=200, height=89, units="mm", device="jpeg", dpi=300)
p <- mygseaplot(em, geneSetID=wi, color = c("#A50026", "#D73027"), subplots = 1:2)
ggsave(p, file=gialli_m, width=200, height=89, units="mm", device="jpeg", dpi=300)


grigio <- c("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION")
grigiocheck <- tot %>% filter(ID %in% grigio)
grigiocheck <- grigiocheck %>% filter(pvalue < 0.05)
res <- rbind(res, grigiocheck[c(1,4,5,6,8)])

wi <- c()

for (w in grigio) {
  i <- which(em@result$Description==w)
  wi <- c(wi, i)
  #print(gseaplot(em, geneSetID = i, title = em$Description[i], by = "runningScore"))
}

## cambiare palette sul grigio
p <- gseaplot2(em, geneSetID=wi, color = c("black"), subplots = 1:2)
ggsave(p, file=grigio_s, width=200, height=89, units="mm", device="pdf", dpi=300)
p <- mygseaplot(em, geneSetID=wi, color = c("black"), subplots = 1:2)
ggsave(p, file=grigio_m, width=200, height=89, units="mm", device="pdf", dpi=300)


verdi1 <- c("HALLMARK_KRAS_SIGNALING_UP", "LIN_APC_TARGETS", "KIM_MYC_AMPLIFICATION_TARGETS_DN",
            "SANSOM_APC_TARGETS_DN", "LEF1_UP.V1_UP", "LEF1_UP.V1_DN", "BCAT.100_UP.V1_UP")
verdi1check <- tot %>% filter(ID %in% verdi1)
verdi1check <- verdi1check %>% filter(pvalue < 0.05)
res <- rbind(res, verdi1check[c(1,4,5,6,8)])

verdipos <- verdi1check %>% filter(enrichmentScore > 0)
verdipos <- verdipos$ID

wi <- c()

for (w in verdipos) {
  i <- which(em@result$Description==w)
  wi <- c(wi, i)
  #print(gseaplot(em, geneSetID = i, title = em$Description[i], by = "runningScore"))
}

## cambiare palette sul verde
p <- gseaplot2(em, geneSetID=wi, color = c("#41AB5D", "#238B45", "#006D2C"), subplots = 1:2)
ggsave(p, file=verdi_s_p, width=200, height=89, units="mm", device="jpeg", dpi=300)
p <- mygseaplot(em, geneSetID=wi, color = c("#41AB5D", "#238B45", "#006D2C"), subplots = 1:2)
ggsave(p, file=verdi_m_p, width=200, height=89, units="mm", device="jpeg", dpi=300)

#"#74C476", "#41AB5D", "#238B45", "#006D2C"

verdineg <- verdi1check %>% filter(enrichmentScore < 0)
verdineg <- verdineg$ID

wi <- c()

for (w in verdineg) {
  i <- which(em@result$Description==w)
  wi <- c(wi, i)
  #print(gseaplot(em, geneSetID = i, title = em$Description[i], by = "runningScore"))
}

## cambiare palette sul verde
p <- gseaplot2(em, geneSetID=wi, color = c("#74C476", "#41AB5D", "#238B45", "#006D2C"), subplots = 1:2)
ggsave(p, file=verdi_s_n, width=200, height=89, units="mm", device="jpeg", dpi=300)
p <- mygseaplot(em, geneSetID=wi, color = c("#74C476", "#41AB5D", "#238B45", "#006D2C"), subplots = 1:2)
ggsave(p, file=verdi_m_n, width=200, height=89, units="mm", device="jpeg", dpi=300)

write.table(res, file=tsv_f, sep = "\t", quote = FALSE, col.names = TRUE, row.names = FALSE)