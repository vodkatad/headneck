library(DescTools)
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

# Example usage with cell counts: a=10, b=20, c=15, d=40
standard_or_ci(a = (mat[1,1]+0.5), b = (mat[1,2]+0.5), c = (mat[2,1]+0.5), d = (mat[2,2]+0.5))
#((a+0.5)*(d+0.5)) / ((b+0.5)*(c+0.5))
#((mat[1,1]+0.5)*(mat[2,2]+0.5)) / ((mat[1,2]+0.5)*(mat[2,1]+0.5))
# 0.79-265.67;