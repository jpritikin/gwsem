library(testthat)
library(gwsem)
library(MASS)

set.seed(1)

dir <- system.file("extdata", package = "gwsem")

pheno <- read.table(file.path(dir, "example.psam"),
                    stringsAsFactors = FALSE,header=TRUE, comment.char="")
colnames(pheno)[1] <- "FID"

numIndicators <- 5
loadings <- rbeta(numIndicators, 4, 3)
resid <- rbeta(numIndicators, 4, 3)^2
indicators <- pheno$phenotype %*% t(loadings) +
  mvrnorm(nrow(pheno), mu=rep(0, numIndicators), Sigma=diag(numIndicators))
colnames(indicators) <- paste0("i", 1:numIndicators)
pheno <- cbind(pheno, indicators)

oneFacGWAS(pheno, file.path(dir,"example"),
           paste0("i", 1:numIndicators), snpFileType="pgen")

result <- read.table("out.log", stringsAsFactors = FALSE, header=TRUE,
                     sep="\t", check.names=FALSE, quote="", comment.char="")
result <- subset(result, statusCode=='OK')
expect_equal(nrow(result), 185)
