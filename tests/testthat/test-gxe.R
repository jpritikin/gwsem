library(testthat)
library(gwsem)
library(MASS)

suppressWarnings(RNGversion("3.5"))
set.seed(1)
#mxOption(key="default optimizer", value="SLSQP")

dir <- system.file("extdata", package = "gwsem")
tdir <- tempdir()

pheno <- read.table(file.path(dir, "example.psam"),
                    stringsAsFactors = FALSE,header=TRUE, comment.char="")
colnames(pheno)[1] <- "FID"
colnames(pheno)[3] <- "sex"

numIndicators <- 7
loadings <- rbeta(numIndicators, 4, 3)
resid <- rbeta(numIndicators, 4, 3)^2
indicators <- pheno$phenotype %*% t(loadings) +
  mvrnorm(nrow(pheno), mu=rep(0, numIndicators), Sigma=diag(numIndicators))
colnames(indicators) <- paste0("i", 1:numIndicators)
pheno <- cbind(pheno, indicators)
#pheno$sex <- pheno$sex - min(pheno$sex) # will blow up

# -----

test_that("gxe data roundtrip", {
  skip_if_not_installed("OpenMx", "2.16.0.1")

  oi <- buildItem(pheno, paste0("i", 3), gxe=c("sex","i4"))
  expect_true(!is.null(oi$snp_sexAlg))
  expect_true(!is.null(oi$snp_i4Alg))
  expect_equal(oi$data$algebra, c("snp_sexAlg", "snp_i4Alg"))

  fit <- GWAS(oi,
              file.path(dir,"example.pgen"),
              file.path(tdir, "out.log"), SNP=c(3))
  ob <- fit$data$observed
  expect_true(all(ob$snp * ob$sex == ob$snp_sex))
  expect_true(all(ob$snp * ob$i4 == ob$snp_i4))

  pheno$snp <- runif(nrow(pheno), 0, 2)
  oi <- buildItem(pheno, paste0("i", 3), gxe=c("sex","i4"))
  fit2 <- mxRun(oi)
  ob <- fit2$data$observed
  expect_true(all(ob$snp * ob$sex == ob$snp_sex))
  expect_true(all(ob$snp * ob$i4 == ob$snp_i4))
})
