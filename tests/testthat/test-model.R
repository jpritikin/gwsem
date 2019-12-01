library(testthat)
library(gwsem)
library(MASS)

# win32 doesn't implement C++ exceptions correctly
skip_if(.Platform$OS.type == "windows" && .Machine$sizeof.pointer == 4)

suppressWarnings(RNGversion("3.5"))
set.seed(1)
#mxOption(key="default optimizer", value="SLSQP")

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

# -----------------

expect_error(buildOneItem(pheno, paste0("i", 1:5)),
             "buildOneItem provided with 5 dependent")

oi <- buildOneItem(pheno, paste0("i", 1))
expect_error(GWAS(oi, "example"),
             "rename snpData")

expect_error(GWAS(oi, file.path(dir,"example.xyz")),
             "Unrecognized file extension")

expect_error(GWAS(oi, something_else=1),
             "Rejected are any values")

GWAS(oi,
     file.path(dir,"example.pgen"),
     file.path(tdir, "out.log"), SNP=c(3:4,6))

pgen <- loadResults(file.path(tdir, "out.log"), "snp2i1")
rx <- which(min(pgen$P) == pgen$P)
expect_equal(rx, 2)
expect_equal(pgen$P[rx], 0.128, tolerance=1e-2)

expect_error(GWAS(oi,
                  file.path(dir,"example.pgen"),
                  file.path(tdir, "out.log"), SNP=c(250)),
             "out of data")

oi <- expect_warning(buildOneItem(pheno, paste0("i", 1),fitfun = "ML",
                   minMAF=.1),
                   "minMAF is ignored when fitfun")

GWAS(oi,
     file.path(dir,"example.pgen"),
     file.path(tdir, "out.log"), SNP=c(3:4,6))

ml <- loadResults(file.path(tdir, "out.log"), "snp2i1")
expect_equal(cor(ml$P, pgen$P), 1, tolerance=1e-3)

expect_error(buildOneItem(pheno, paste0("i", 1), fitfun = "bob"),
             "should be one of")

# -----------------

if (.Platform$OS.type == "windows" && packageVersion("data.table") <= '1.12.2') {
  skip("data.table has a bug on windows")
}

GWAS(buildOneFac(pheno, paste0("i", 1:numIndicators)),
     file.path(dir,"example.pgen"),
     file.path(tdir, "out.log"), SNP=c(3:4,6))

pgen <- loadResults(file.path(tdir, "out.log"), "snp2F")
pgen2 <- loadResults(rep(file.path(tdir, "out.log"), 2), 'snp2F')

expect_equal(pgen$SNP, paste0("RSID_", c(4,5,7)))
expect_equal(pgen2$SNP, rep(paste0("RSID_", c(4,5,7)), 2))
expect_equal(pgen$snp2F, c(-.188, -.066, .006), tolerance=.001)

pgen <- read.table(file.path(tdir, "out.log"), stringsAsFactors = FALSE, header=TRUE,
                   sep="\t", check.names=FALSE, quote="", comment.char="")

l2 <- pgen[,paste0('lambda_i',1:7)]
expect_equivalent(colMeans(l2) / loadings, rep(1, numIndicators), tolerance=.2)

if(0) {
  m1 <- buildOneFac(pheno, file.path(dir,"example.pgen"), paste0("i", 1:numIndicators), SNP=3)
  m1 <- mxRun(m1)
}

# -----------------

GWAS(buildOneFacRes(pheno, paste0("i", 1:numIndicators)),
     file.path(dir,"example.pgen"),
     file.path(tdir, "out.log"), SNP=c(3:4,6))

pgen <- read.table(file.path(tdir, "out.log"), stringsAsFactors = FALSE, header=TRUE,
                   sep="\t", check.names=FALSE, quote="", comment.char="")
expect_equivalent(pgen$i1_Thr_1, rep(quantile(origPheno$i1, 1/3),3), tolerance=.18)
expect_equivalent(pgen$i2_Thr_1, rep(0,3), tolerance=0.15)

l2 <- pgen[,paste0('lambda_i',1:7)]
expect_equivalent(colMeans(l2) / loadings, rep(1, numIndicators), tolerance=.2)

pgen <- loadResults(file.path(tdir, "out.log"), "snp2i2")
expect_equal(pgen$P, rep(1,3), tolerance=1e-5)

# -----------------

GWAS(buildTwoFac(pheno, F1itemNames = paste0("i",1:4),
                 F2itemNames = paste0("i",4:7)),
     file.path(dir,"example.pgen"),
     file.path(tdir, "out.log"), SNP=c(3,4,6))

pgen <- read.table(file.path(tdir, "out.log"), stringsAsFactors = FALSE, header=TRUE,
                   sep="\t", check.names=FALSE, quote="", comment.char="")
l2 <- pgen[,paste0('lambda_i',1:7)]
expect_equivalent((colMeans(l2) / loadings)[-4], rep(1, numIndicators-1),
                  tolerance=.15)
expect_equal(pgen[['TwoFac.S[9,10]']], rep(1.17,3), .02)

# ------------

for (f in c("buildOneFac", "buildOneFacRes", "buildOneItem",
               "buildTwoFac")) {
  expect_error(do.call(f, args=list(pheno, another.arg="xyz")),
               "Rejected are any values")
}

