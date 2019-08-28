library(testthat)
library(gwsem)
library(MASS)

set.seed(1)

dir <- '/scratch/gw-sem/'

pheno <- read.table(file.path(dir, "hapmap1.fam"), stringsAsFactors = FALSE)
colnames(pheno)[1] <- c('id')
colnames(pheno)[6] <- c('disease')

numIndicators <- 5
loadings <- rbeta(numIndicators, 4, 3)
resid <- rbeta(numIndicators, 4, 3)^2
indicators <- (pheno$disease-1) %*% t(loadings) +
  mvrnorm(nrow(pheno), mu=rep(0, numIndicators), Sigma=diag(numIndicators))
colnames(indicators) <- paste0("i", 1:numIndicators)
pheno <- cbind(pheno, indicators)

oneFacGWAS(pheno, file.path(dir,"hapmap1.bed"),
           paste0("i", 1:numIndicators), nSNP=10, snpFileType="pgen")
