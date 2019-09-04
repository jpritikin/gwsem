library(testthat)
library(gwsem)
library(MASS)

suppressWarnings(RNGversion("3.5"))
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

numCovariate <- 2
for (cx in 1:numCovariate) {
  pheno[[paste0("covar", cx)]] <- rnorm(nrow(pheno))
}

expect_error(oneFacGWAS(pheno, file.path(dir,"example.pgen"), paste0("i", 1:numIndicators),
                        covariates = paste0("covar",1:numCovariate)),
             "record 200 requested but only 199 in file")

pgen <- read.table("out.log", stringsAsFactors = FALSE, header=TRUE,
                   sep="\t", check.names=FALSE, quote="", comment.char="")

mask <- (pgen$catch1 == "" & pgen$statusCode=="OK" & !is.na(pgen$snpRegSE))
pgen <- pgen[mask,]

mask <- (abs(pgen$snpReg) < 2.6*mad(pgen$snpReg) & (abs(pgen$snpReg / pgen$snpRegSE) > .5))
pgen <- pgen[mask,]

cvNames <- paste(rep(paste0("covar",1:numCovariate), each = numIndicators),
      paste0("i", 1:numIndicators), sep = "_")
expect_equivalent(colMeans(pgen[,cvNames]), rep(0, length(cvNames)), tolerance=.1)
