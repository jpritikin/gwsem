library(testthat)
library(gwsem)
library(MASS)

skip_if(Sys.info()[["machine"]] != 'x86_64')

suppressWarnings(RNGversion("3.5"))
set.seed(1)
#mxOption(key="default optimizer", value="SLSQP")

dir <- system.file("extdata", package = "gwsem")
tdir <- tempdir()

pheno <- read.table(file.path(dir, "example.psam"),
                    as.is = TRUE,header=TRUE, comment.char="")
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

test_that("moderator level", {
  oi <- buildItem(pheno, paste0("i", 3), gxe=c("i4"))
  fit <- GWAS(oi,
              file.path(dir,"example.pgen"),
              file.path(tdir, "out.log"))
  m1 <- loadResults(file.path(tdir, "out.log"), 'snp_i4_to_i3')
  m1.5 <- signifGxE(m1, 'snp_i4_to_i3', level = .5)
  m1.7 <- signifGxE(m1, 'snp_i4_to_i3', level = .7)
  expect_equal(fivenum(m1.5$Z - m1.7$Z),
               c(-0.55, -0.34, -0.27, -0.19, 0.04), .01)
  expect_equivalent(c(table(isSuspicious(m1, c('snp_to_i3', 'snp_i4_to_i3')))),
               c(197, 2))
})
