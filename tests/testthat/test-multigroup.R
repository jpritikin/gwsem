library(testthat)
library(gwsem)
library(MASS)

suppressWarnings(RNGversion("3.5"))
set.seed(1)
#mxOption(key="default optimizer", value="SLSQP")

dir <- system.file("extdata", package = "gwsem")
tdir <- tempdir()

pheno <- read.table(file.path(dir, "example.psam"),
                    as.is = TRUE,header=TRUE, comment.char="")
colnames(pheno)[1] <- "FID"

numIndicators <- 7

createIndicators <- function(pheno) {
  loadings <- rbeta(numIndicators, 4, 3)
  resid <- rbeta(numIndicators, 4, 3)^2
  indicators <- pheno$phenotype %*% t(loadings) +
    mvrnorm(nrow(pheno), mu=rep(0, numIndicators), Sigma=diag(numIndicators))
  colnames(indicators) <- paste0("i", 1:numIndicators)
  indicators <- as.data.frame(indicators)
  indicators$i1 <- cut(indicators$i1, c(-Inf, quantile(indicators$i1, 1:2/3), Inf),
                  ordered_result = TRUE)
  indicators$i2 <- cut(indicators$i2, c(-Inf, quantile(indicators$i2, .5), Inf),
                  ordered_result = TRUE)
  cbind(pheno, indicators)
}

# -----------------

oi <- mxModel(
  "mg",
  mxModel(name="g1", buildOneFac(createIndicators(pheno),
                                 paste0("i", 1:numIndicators))),
  mxModel(name="g2", buildOneFac(createIndicators(pheno),
                                 paste0("i", 1:numIndicators))),
  mxFitFunctionMultigroup(paste0('g',1:2)))

# clear labels
oi <- omxSetParameters(oi, names(coef(oi)),
                       newlabels = rep(NA_character_, length(coef(oi))))

# run both groups together

z1 <- GWAS(
  oi,
  c(g1=file.path(dir,"example.pgen"), g2=file.path(dir,"example.pgen")),
  file.path(tdir, "out.log"), SNP=1:10)

r1 <- loadResults(file.path(tdir, "out.log"), "g1.A[9,1]", signAdj='g1.A[2,9]')
r2 <- loadResults(file.path(tdir, "out.log"), "g2.A[9,1]", signAdj='g2.A[2,9]')

# run each group separately

z1 <- GWAS(oi$g1, file.path(dir,"example.pgen"),
           file.path(tdir, "out1.log"), SNP=1:10)
s1 <- loadResults(file.path(tdir, "out1.log"), "g1.A[9,1]", signAdj='g1.A[2,9]')
expect_equal(s1$Z - r1$Z, rep(0, nrow(s1)), 1e-5)

z1 <- GWAS(oi$g2, file.path(dir,"example.pgen"),
           file.path(tdir, "out2.log"), SNP=1:10)
s2 <- loadResults(file.path(tdir, "out2.log"), "g2.A[9,1]", signAdj='g2.A[2,9]')
expect_equal(s2$Z - r2$Z, rep(0, nrow(s1)), 1e-5)
