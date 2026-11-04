#!/usr/bin/Rscript
library(S4Vectors)
## manual inputs of ck for lazyness, Fra S has prism files (are there libraries?)
ck <- list()
ck[[1]] <- matrix(c(34,13, 5,10), nrow=2, byrow=FALSE) # CK1
ck[[2]] <- matrix(c(36,16,3,7), nrow=2, byrow=FALSE) # CK5
ck[[3]] <- matrix(c(23,8,16,15), nrow=2, byrow=FALSE) # CK10
ck[[4]] <- matrix(c(20,15,19,8), nrow=2, byrow=FALSE) # CK8
ck[[5]] <- matrix(c(5,5,34,18), nrow=2, byrow=FALSE) # CK18
names(ck) <- c('CK1', 'CK5', 'CK10', 'CK8', 'CK18')

for (i in seq(1, length(ck))) {
  rownames(ck[[i]]) <- c('Resp', 'Non-Resp')
  colnames(ck[[i]]) <- c('High', 'Low')
}

standard_or_ci_correctedzero <- function(a, b, c, d, conf.level = 0.95) {
  # 1. Calculate standard sample odds ratio
  if (any(c(a==0, b==0, c==0, d==0))) {
    a <- a + 0.5
    b <- b + 0.5
    c <- c + 0.5
    d <- d + 0.5
  }
  OR <- (a * d) / (b * c)
  # 2. Calculate standard error of the log OR
  se_log_OR <- sqrt(1/a + 1/b + 1/c + 1/d)
  # 3. Get the critical Z value
  alpha <- 1 - conf.level
  z <- qnorm(1 - alpha/2)
  # 4. Calculate confidence bounds on the log scale, then exponentiate
  lower <- exp(log(OR) - z * se_log_OR)
  upper <- exp(log(OR) + z * se_log_OR)
  # Return results
  return(list(
    odds_ratio = OR,
    cl_conf.int = c(lower, upper),
    cl_conf.level = conf.level
  ))
}

our_fisher_standard_or_and_ic <- function(mat, conf.level=0.95) {
  a <- mat[1,1]
  b <- mat[1,2]
  c <- mat[2,1]
  d <- mat[2,2]
  f <- fisher.test(mat, conf.level = conf.level)
  or <- standard_or_ci_correctedzero(a=a, b=b, c=c, d=d, conf.level=0.95)
  return(c(f, or))
}
res <- data.frame(matrix(NA, nrow=length(ck), ncol=8))
colnames(res) <- c('fisher_pvalue', 'OR_cl', 'OR_low_cl', 'OR_high_cl', 'OR_MLE', 'OR_low_MLE', 'OR_high_MLE', 'padj' )
for (i in seq(1, length(ck))) {
  fi <- our_fisher_standard_or_and_ic(ck[[i]])
  res[i, 'fisher_pvalue'] <- fi$p.value
  res[i, 'OR_cl'] <- fi$odds_ratio
  res[i, 'OR_low_cl'] <- fi$cl_conf.int[1]
  res[i, 'OR_high_cl'] <- fi$cl_conf.int[2]
  res[i, 'OR_MLE'] <- fi$estimate
  res[i, 'OR_low_MLE'] <- fi$conf.int[1]
  res[i, 'OR_high_MLE'] <- fi$conf.int[2]
}
res$padj <- p.adjust(res$fisher_pvalue, method='BH') 
res <- endoapply(res, function(x) {round(x, digits=3)})
rownames(res) <- names(ck)

write.table(res, row.names=TRUE, quote=FALSE, sep="\t")
##
q(save="no")

# previous pvalues
library(DescTools)
MHChisqTest(ck$CK8)
MHChisqTest(ck$CK18)
our_fisher_standard_or_and_ic(ck[[1]]) # in the old text we had 'classical' odds ratio not MLE ones. I will start plugging in them.
###

MHChisqTest(matrix(c(5,17,18,4), nrow=2, byrow=F)) # check chemio
chisq.test(matrix(c(5,17,18,4), nrow=2, byrow=F))
chisq.test(matrix(c(5,10,34,13), nrow=2, byrow=F))

## ck nuove revision
MHChisqTest(matrix(c(5,10,34,13), nrow=2, byrow=F)) 
MHChisqTest(matrix(c(3,7,36,16), nrow=2, byrow=F))
MHChisqTest(matrix(c(19,8,20,15), nrow=2, byrow=F))
MHChisqTest(matrix(c(16,15,23,8), nrow=2, byrow=F))
MHChisqTest(matrix(c(34,18,5,5), nrow=2, byrow=F))

fisher.test(matrix(c(5,10,34,13), nrow=2, byrow=F)) 
fisher.test(matrix(c(3,7,36,16), nrow=2, byrow=F))
fisher.test(matrix(c(19,8,20,15), nrow=2, byrow=F))
fisher.test(matrix(c(16,15,23,8), nrow=2, byrow=F))
fisher.test(matrix(c(34,18,5,5), nrow=2, byrow=F))


MHChisqTest(matrix(c(5,10,34,13), nrow=2, byrow=F)) 
chisq.test(matrix(c(5,10,34,13), nrow=2, byrow=F))


## FAT1
fisher.test(matrix(c(8, 0, 24, 13), nrow=2))
MHChisqTest(matrix(c(8, 0, 24, 13), nrow=2))
chisq.test(matrix(c(8, 0, 24, 13), nrow=2))
chisq.test(matrix(c(8, 0, 24, 13), nrow=2), correct=F )
MHChisqTest(matrix(c(8, 0, 24, 13), nrow=2), srow=rep(1,4), scol=rep(1,4))
# ma non ho F/M nella figura quindi noia

# EGFR AMPL OR 14.44 p 0.015
#Correzione di Haldane-Anscombe (pseudo-conteggi):   ((a+0.5)*(d+0.5)) / ((b+0.5)*(c+0.5))
mat <- matrix(c(11,0,21,13), nrow=2)
MHChisqTest(mat) # is this 0.015?
chisq.test(mat)
fisher.test(mat)

#or torna, pvalue mah.
((mat[1,1]+0.5)*(mat[2,2]+0.5)) / ((mat[1,2]+0.5)*(mat[2,1]+0.5))


f <- matrix(c(2, 0, 7, 3), nrow=2)
m <- matrix(c(9, 0, 14, 10), nrow=2)

x <- array(c(f, m), dim=c(2,2,2))
mantelhaen.test(x) # 0.015 boh :()

# N / N-1 etc
chisq.test(matrix(c(80, 100, 50, 50), nrow=2), correct=F )
MHChisqTest(matrix(c(80, 100, 50, 50), nrow=2))

### IC for or - checked vs EGFR by MK
standard_or_ci <- function(a, b, c, d, conf.level = 0.95) {
  # 1. Calculate standard sample odds ratio
  OR <- (a * d) / (b * c)
  # 2. Calculate standard error of the log OR
  se_log_OR <- sqrt(1/a + 1/b + 1/c + 1/d)
  # 3. Get the critical Z value
  alpha <- 1 - conf.level
  z <- qnorm(1 - alpha/2)
  # 4. Calculate confidence bounds on the log scale, then exponentiate
  lower <- exp(log(OR) - z * se_log_OR)
  upper <- exp(log(OR) + z * se_log_OR)
  # Return results
  return(list(
    odds_ratio = OR,
    conf.int = c(lower, upper),
    conf.level = conf.level
  ))
}

our_fisher_standard_or_and_ic <- function(mat, conf.level=0.95) {
  a <- mat[1,1]
  b <- mat[1,2]
  c <- mat[2,1]
  d <- mat[2,2]
  f <- fisher.test(mat, conf.level = conf.level)
  or <- standard_or_ci_correctedzero(a=a, b=b, c=c, d=d, conf.level=0.95)
  return(c(f, or))
}

# ci atteniamo amettere OR del fisher.test classico anche noi

## TODO test one old from chemio or other
# FSIP2: 
## 10 11 1.121 P = 1.00 
our_fisher_standard_or_and_ic(matrix(c(10, 51, 11, 50), nrow=2, byrow=F))
our_fisher_standard_or_and_ic(matrix(c(51, 10, 50, 11), nrow=2, byrow=F))

#8 12 1.616 P = 0.464 # qui mk ha messo odds ratio del fisher dritto con MLE non quello stimato normale.
our_fisher_standard_or_and_ic(matrix(c(53, 8, 49, 12), nrow=2, byrow=F))

# old chemo
our_fisher_standard_or_and_ic(matrix(c(13, 28, 28, 13), nrow=2, byrow=F))

# 

# Example usage with cell counts: a=10, b=20, c=15, d=40
standard_or_ci(a = (mat[1,1]+0.5), b = (mat[1,2]+0.5), c = (mat[2,1]+0.5), d = (mat[2,2]+0.5))
standard_or_ci_correctedzero(a = mat[1,1], b = mat[1,2], c = mat[2,1], d = mat[2,2])

#((a+0.5)*(d+0.5)) / ((b+0.5)*(c+0.5))

#((mat[1,1]+0.5)*(mat[2,2]+0.5)) / ((mat[1,2]+0.5)*(mat[2,1]+0.5))
# 0.79-265.67;

#                          Non Responder  Responder
#BA                    2                         19
#else                  8                         11

#Non Responder Responder
#IMS                  0                        2
#else                 10                       28

#         Non Responder Responder
#cor_CL             8         9
#else               2        21

matBA <- matrix(c(19,11,2,8), ncol=2)
matIMS <- matrix(c(2,28,0,10), ncol=2)
matCL <- matrix(c(9,21,8,2), ncol=2)
our_fisher_standard_or_and_ic(matBA)
our_fisher_standard_or_and_ic(matIMS)
our_fisher_standard_or_and_ic(matCL)

MHChisqTest(matBA)
MHChisqTest(matIMS)
MHChisqTest(matCL)

MHChisqTest(matBA, srow = c(1,2), scol = c(1,2))
            

m <- rbind(matBA[1,], matCL[1,], matIMS[1,])
MHChisqTest(m)