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
                    as.is = TRUE,header=TRUE, comment.char="")
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

oi <- buildItem(pheno, paste0("i", 3), exogenous = TRUE)
expect_equivalent(oi$M$labels[1,'snp'], "data.snp")

oi <- buildItem(pheno, paste0("i", 3))
expect_true(is.null(oi$M))
expect_true(oi$A$free['i3', 'snp'])
expect_error(GWAS(oi, "example"),
             "rename snpData")

expect_error(GWAS(oi, file.path(dir,"example.xyz")),
             "Unrecognized file extension")

expect_error(GWAS(oi, something_else=1),
             "Rejected are any values")

got6 <- GWAS(oi,
     file.path(dir,"example.pgen"),
     file.path(tdir, "out.log"), SNP=c(3:5))
expect_true(is.null(got6$data$observedStats$means))

pgen <- loadResults(file.path(tdir, "out.log"), "snp_to_i3")
pgen <- signif(pgen, "snp_to_i3")
rx <- which(min(pgen$P) == pgen$P)
expect_equal(rx, 3)
expect_equal(pgen$P[rx], .155, tolerance=1e-2)

test_that("out of data", {
  skip_on_os('mac') # exceptions are broken
  expect_error(GWAS(oi,
                    file.path(dir,"example.pgen"),
                    file.path(tdir, "out.log"), SNP=c(250)),
               "out of data")
})

oi <- expect_warning(buildItem(pheno, paste0("i", 3), fitfun = "ML",
                   minMAF=.1),
                   "minMAF is ignored when fitfun")
expect_true(!is.null(oi$M))

GWAS(oi,
     file.path(dir,"example.pgen"),
     file.path(tdir, "out.log"), SNP=c(3:5))

ml <- loadResults(file.path(tdir, "out.log"), "snp_to_i3")
ml <- signif(ml, "snp_to_i3")
expect_equal(cor(ml$P, pgen$P), 1, tolerance=1e-3)

expect_error(buildItem(pheno, paste0("i", 1), fitfun = "bob"),
             "should be one of")

# -----------------

oi <- buildItem(pheno, paste0("i", 2:3))
expect_true(oi$S$free['i2','i3'])

GWAS(oi,
     file.path(dir,"example.pgen"),
     file.path(tdir, "out.log"), SNP=c(11:20))

pgen <- loadResults(file.path(tdir, "out.log"), "snp_to_i2")
pgen <- signif(pgen, "snp_to_i2")
rx <- which(min(pgen$P) == pgen$P)
expect_equal(rx, 10)
expect_equal(pgen$P[rx], .0251, tolerance=1e-2)

pgen <- loadResults(file.path(tdir, "out.log"), "snp_to_i3")
pgen <- signif(pgen, "snp_to_i3")
rx <- which(min(pgen$P) == pgen$P)
expect_equal(rx, 4)
expect_equal(pgen$P[rx], .086, tolerance=1e-2)

# -----------------

if (.Platform$OS.type == "windows" && packageVersion("data.table") <= '1.12.2') {
  skip("data.table has a bug on windows")
}

z1 = GWAS(buildOneFac(pheno, paste0("i", 1:numIndicators)),
     file.path(dir,"example.pgen"),
     file.path(tdir, "out.log"), SNP=c(3:5))

pgen <- loadResults(file.path(tdir, "out.log"), c("snp_to_F", 'lambda_i1'))
pgen <- signif(pgen, "snp_to_F", signAdj = "lambda_i1")
pgen2 <- loadResults(rep(file.path(tdir, "out.log"), 2), c('snp_to_F', 'lambda_i1'))
pgen2 <- signif(pgen2, "snp_to_F", signAdj = "lambda_i1")

expect_equal(pgen$SNP, paste0("RSID_", 4:6))
expect_equal(pgen2$SNP, rep(paste0("RSID_", 4:6), 2))
expect_equal(pgen$snp_to_F, c(.188, .066, -.081), tolerance=.001)

pgen <- read.table(file.path(tdir, "out.log"), as.is = TRUE, header=TRUE,
                   sep="\t", check.names=FALSE, quote="", comment.char="")

l2 <- pgen[,paste0('lambda_i',1:7)]
expect_equivalent(colMeans(l2) / loadings, rep(1, numIndicators), tolerance=.2)

if(0) {
  m1 <- buildOneFac(pheno, file.path(dir,"example.pgen"), paste0("i", 1:numIndicators), SNP=3)
  m1 <- mxRun(m1)
}

# -----------------

m1 <- buildOneFacRes(pheno, paste0("i", 1:numIndicators), factor=TRUE)
expect_true(m1$A$free['F','snp'])

m2 <- buildOneFacRes(pheno, paste0("i", 1:numIndicators))
m2 <- GWAS(m2,
     file.path(dir,"example.pgen"),
     file.path(tdir, "out.log"), SNP=c(3:5))
expect_false(m2$A$free['F','snp'])

pgen <- read.table(file.path(tdir, "out.log"), as.is = TRUE, header=TRUE,
                   sep="\t", check.names=FALSE, quote="", comment.char="")
expect_equivalent(pgen$i1_Thr_1, rep(quantile(origPheno$i1, 1/3),3), tolerance=.18)
expect_equivalent(pgen$i2_Thr_1, rep(0,3), tolerance=0.15)

l2 <- pgen[,paste0('lambda_i',1:7)]
expect_equivalent(colMeans(l2) / loadings, rep(1, numIndicators), tolerance=.2)

pgen <- loadResults(file.path(tdir, "out.log"), "snp_to_i2")
pgen <- signif(pgen, "snp_to_i2")
expect_equal(pgen$P, c(0.481, 0.896, 0.412), tolerance=1e-2)

# -----------------

GWAS(buildTwoFac(pheno, F1itemNames = paste0("i",1:4),
                 F2itemNames = paste0("i",4:7)),
     file.path(dir,"example.pgen"),
     file.path(tdir, "out.log"), SNP=3:5)

pgen <- read.table(file.path(tdir, "out.log"), as.is = TRUE, header=TRUE,
                   sep="\t", check.names=FALSE, quote="", comment.char="")
l2 <- pgen[,c(paste0('F1_lambda_i',1:3), paste0('F2_lambda_i',5:7))]
expect_equivalent((colMeans(l2) / loadings[-4]), rep(1, numIndicators-1),
                  tolerance=.15)
expect_equal(pgen[['facCov']], rep(1.17,3), .02)

# ------------

for (f in c("buildOneFac", "buildOneFacRes", "buildItem",
               "buildTwoFac")) {
  expect_error(do.call(f, args=list(pheno, another.arg="xyz")),
               "Rejected are any values")
}

