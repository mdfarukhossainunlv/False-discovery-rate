if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("multtest")
BiocManager::install("qvalue")
BiocManager::install("siggenes")
require(multtest)
data(golub)
class(golub)
dim(golub)
head(golub)
dim(golub.gnames)
golub.gnames[1:4, ]
golub.cl

###  Computing Welch t-test
teststat = mt.teststat(golub, golub.cl)

require(ggplot2)
plt = ggplot(data.frame(teststat), aes(sample = teststat)) + stat_qq() + theme_bw()
plt

### Adjusting p values
##First we will compute raw nominal two-sided p-values
rawp = 2 * (1 - pnorm(abs(teststat)))
## Adjust
procedures = c( "BH", "BY")
adjusted = mt.rawp2adjp(rawp, procedures)
adjusted$adjp[1:10, ]

# Display them in original data order
adjusted$adj[order(adjusted$index)[1:10], ]

## Calculation of p values
resT = mt.maxT(golub, golub.cl, B = 10000)
names(resT)

ord = order(resT$index)
rawp = resT$rawp[ord]
maxT = resT$adjp[ord]
teststat = resT$teststat

#### The p-values are sorted in decreasing order of the absolute values of the test statistics.
head(resT$teststat)


### plot (histogram)
library(ggplot2)
perms = as.numeric(mt.sample.teststat(golub[1, ], golub.cl, B = 10000))
plt = ggplot(data.frame(perms), aes(x = perms)) + stat_bin(position = "identity")
plt

### q-value using the p.adjust function with the FDR option ( FDR < .05)

library(multtest)
library(siggenes)
data(golub)

# Perform a SAM analysis.
sam.out<-sam(golub,golub.cl,B=100,rand=123)

# Estimate the prior probability that a gene is not significant.
pi0<-pi0.est(sam.out@p.value)$p0

# Compute the q-values of the genes.
q.value<-qvalue.cal(sam.out@p.value,pi0)
head(q.value, n=10)
# how many 1 value
sum(q.value<0.05)

#### control FDR
fdr = MTP(matrix(golub[1, ], nrow = 1), Y = golub.cl, typeone = "fdr")
summary(fdr)
