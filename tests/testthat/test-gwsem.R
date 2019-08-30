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

expect_error(oneFacGWAS(pheno, file.path(dir,"example.pgen"),
           paste0("i", 1:numIndicators)),
           "record 200 requested but only 199 in file")

pgen <- read.table("out.log", stringsAsFactors = FALSE, header=TRUE,
                     sep="\t", check.names=FALSE, quote="", comment.char="")
pgen <- subset(pgen, statusCode=='OK')
#expect_equal(nrow(pgen), 195)

expect_error(oneFacGWAS(pheno, file.path(dir,"example.bed"),
                        paste0("i", 1:numIndicators)),
             "record 200 requested but only 199 in file")
bed <- read.table("out.log", stringsAsFactors = FALSE, header=TRUE,
                   sep="\t", check.names=FALSE, quote="", comment.char="")
bed <- subset(bed, statusCode=='OK')
#expect_equal(nrow(bed), 195)

expect_error(oneFacGWAS(pheno, file.path(dir,"example.bgen"),
                        paste0("i", 1:numIndicators)),
             "has no more varients")
bgen <- read.table("out.log", stringsAsFactors = FALSE, header=TRUE,
                  sep="\t", check.names=FALSE, quote="", comment.char="")
bgen <- subset(bgen, statusCode=='OK')
#expect_equal(nrow(bgen), 197)

mask <- abs(pgen$snpReg) < 2.6*mad(pgen$snp_res)
#cor(pgen$snpReg[mask], bed$snpReg[mask])

#hist(pgen$snpReg[mask] / bed$snpReg[mask])


# loadCounter TODO
