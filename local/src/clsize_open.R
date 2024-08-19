library(tidyverse)


expr_f <- snakemake@input[["expr"]]
samples_data_f <- snakemake@input[["samples"]]
cls_f <- snakemake@output[["clsf"]]

what <- snakemake@wildcards[['what']]
nom <- snakemake@wildcards[['nom']]
den <- snakemake@wildcards[['den']]

###
#expr_f <- '/scratch/trcanmed/DE_RNASeq/dataset/kras_deg/tmm.tsv.gz'
#samples_data_f <- '/scratch/trcanmed/DE_RNASeq/dataset/kras_deg/samples_data'
#what <- 'type'
#nom <- 'Ind_4Q'
#den <- 'Dep_1Q'
#nom <- "PD"
#den <- "PR"


samples_data <- read.table(samples_data_f, header=TRUE, sep="\t", stringsAsFactors = FALSE)
expr_data <- read.table(expr_f, header=TRUE, sep="\t", stringsAsFactors = FALSE)

classes <- unique(samples_data[,what])
sink(cls_f)
cat(paste(nrow(samples_data), length(classes), 1))
cat("\n")
sink()

rownames(samples_data) <- samples_data$SampleID

first_class <- samples_data[rownames(samples_data)==colnames(expr_data)[1], what]
ordered_classes <- c(first_class, classes[classes != first_class]) # TODO test if more than two classes

sink(cls_f, append=TRUE)
cat(paste(c('#', ordered_classes), collapse=" "))
cat("\n")
sink()

ordered_samples <- samples_data[match(colnames(expr_data), rownames(samples_data)),]

if (! all(colnames(expr_data)==rownames(ordered_samples))) {
  stop('Something off in the match')
}

ordered_samples$labels <- ifelse(ordered_samples[, what] == first_class, 0, 1)

sink(cls_f, append=TRUE)
cat(paste(ordered_samples$labels, collapse=" "))
cat("\n")
sink()

