library(testthat)
library(gwsem)
library(MASS)

suppressWarnings(RNGversion("3.5"))
set.seed(1)

dir <- system.file("extdata", package = "gwsem")
tdir <- tempdir()

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

GWAS(buildOneFac(pheno, paste0("i", 1:numIndicators),
                 exogenousCovariates = paste0("covar",1:numCovariate)),
     file.path(dir,"example.pgen"),
     file.path(tdir,"out.log"))

pgen <- read.table(file.path(tdir,"out.log"), stringsAsFactors = FALSE, header=TRUE,
                   sep="\t", check.names=FALSE, quote="", comment.char="")

mask <- (pgen$catch1 == "" & pgen$statusCode=="OK" & !is.na(pgen$snp2FSE))
pgen <- pgen[mask,]

mask <- (abs(pgen$snp2F) < 2.6*mad(pgen$snp2F) & (abs(pgen$snp2F / pgen$snp2FSE) > .5))
pgen <- pgen[mask,]

cvNames <- paste(rep(paste0("covar",1:numCovariate), each = numIndicators),
      paste0("i", 1:numIndicators), sep = "_")
expect_equivalent(colMeans(pgen[,cvNames]), rep(0, length(cvNames)), tolerance=.1)

pgen <- loadResults(file.path(tdir,"out.log"), "snp2F")
expect_error(plot(pgen, y=1),
             "plot does not accept a y= argument")

# -----

GWAS(buildOneFac(pheno, paste0("i", 1:numIndicators),
                 covariates = paste0("covar",1:numCovariate)),
     file.path(dir,"example.pgen"),
     file.path(tdir,"out.log"))

pgen2 <- loadResults(file.path(tdir,"out.log"), "snp2F")

expect_equal(cor(pgen$Z, pgen2$Z, use = "pairwise"), 1, tolerance=.15)
