library(testthat)
library(gwsem)
library(MASS)

set.seed(1)

dir <- system.file("extdata", package = "gwsem")

pheno <- read.table(file.path(dir, "example.psam"),
                    stringsAsFactors = FALSE,header=TRUE, comment.char="")
colnames(pheno)[1] <- "FID"

numIndicators <- 7
loadings <- rbeta(numIndicators, 4, 3)
resid <- rbeta(numIndicators, 4, 3)^2
indicators <- pheno$phenotype %*% t(loadings) +
  mvrnorm(nrow(pheno), mu=rep(0, numIndicators), Sigma=diag(numIndicators))
colnames(indicators) <- paste0("i", 1:numIndicators)
pheno <- cbind(pheno, indicators)
origPheno <- pheno
pheno$i1 <- cut(pheno$i1, c(-Inf, quantile(pheno$i1, 1:2/3), Inf), ordered_result = TRUE)
pheno$i2 <- cut(pheno$i2, c(-Inf, quantile(pheno$i2, .5), Inf), ordered_result = TRUE)

oneFacGWAS(pheno, file.path(dir,"example.pgen"), paste0("i", 1:numIndicators), SNP=c(3:4,6))

pgen <- read.table("out.log", stringsAsFactors = FALSE, header=TRUE,
                   sep="\t", check.names=FALSE, quote="", comment.char="")

expect_equal(pgen$ID, paste0("RSID_", c(4,5,7)))
expect_equal(pgen$snpReg, c(.2,.06, .02), tolerance=.02)
expect_equivalent(pgen$i1_Thr_1, rep(quantile(origPheno$i1, 1/3),3), tolerance=.05)
expect_equivalent(pgen$i1_Thr_2, rep(quantile(origPheno$i1, 2/3),3), tolerance=.12)
expect_equivalent(pgen$i2_Thr_1, rep(0,3), tolerance=0.05)

if(0) {
  m1 <- buildOneFac(pheno, file.path(dir,"example.pgen"), paste0("i", 1:numIndicators), SNP=3)
  m1 <- mxRun(m1)
}
