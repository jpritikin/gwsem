library(testthat)
library(gwsem)
library(MASS)

set.seed(1)

dir <- system.file("extdata", package = "gwsem")
tdir <- tempdir()

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

# -----

m1 <- GWAS(buildOneFac(pheno, paste0("i", 1:numIndicators)),
           file.path(dir,"example.pgen"),
           file.path(tdir, "out.log"))
rawSNP <- m1$data$observed$snp
pgen <- read.table(file.path(tdir, "out.log"), stringsAsFactors = FALSE, header=TRUE,
                     sep="\t", check.names=FALSE, quote="", comment.char="")
expect_equal(nrow(pgen), 199)
expect_equal(m1$compute$steps$LD$debug$loadCounter, 1)
expect_equal(range(m1$data$observed$snp), c(0,2), .01)

m1 <- GWAS(buildOneFac(pheno, paste0("i", 1:numIndicators)),
     file.path(dir,"example.bed"),
     file.path(tdir, "out.log"))
# bed can't store dosages so we get ordinal values with some NAs
mask <- !is.na(m1$data$observed$snp)
expect_equal(m1$data$observed$snp[mask] - rawSNP[mask],
             rep(0, sum(mask)), .1)
bed <- read.table(file.path(tdir, "out.log"), stringsAsFactors = FALSE, header=TRUE,
                   sep="\t", check.names=FALSE, quote="", comment.char="")
expect_equal(nrow(bed), 199)
expect_equal(m1$compute$steps$LD$debug$loadCounter, 1)
lr <- loadResults(file.path(tdir, "out.log"), focus = "snp2F")
expect_equal(colnames(lr), c("MxComputeLoop1", "CHR", "BP", "SNP", "A1", "A2",
                             "statusCode", "catch1", "snp2F", "Z", "P"))
expect_equal(range(m1$data$observed$snp, na.rm = TRUE), c(0,2), .01)

m1 <- GWAS(buildOneFac(pheno, paste0("i", 1:numIndicators)),
     file.path(dir,"example.bgen"),
     file.path(tdir, "out.log"))
expect_equal(rawSNP, m1$data$observed$snp, tolerance=5e-5)
bgen <- read.table(file.path(tdir, "out.log"), stringsAsFactors = FALSE, header=TRUE,
                  sep="\t", check.names=FALSE, quote="", comment.char="")
expect_equal(nrow(bgen), 199)
expect_equal(m1$compute$steps$LD$debug$loadCounter, 1)
lr <- loadResults(file.path(tdir, "out.log"), "snp2F")
expect_equal(colnames(lr), c("MxComputeLoop1", "CHR", "BP", "SNP", "A1", "A2",
                             "statusCode", "catch1", "snp2F", "Z", "P"))
expect_equal(range(m1$data$observed$snp), c(0,2), .01)

# -----

m2 <- GWAS(buildOneFac(pheno, paste0("i", 1:numIndicators)),
           file.path(dir,"example.bgen"),
           file.path(tdir, "out.log"), startFrom=190)
last <- read.table(file.path(tdir, "out.log"), stringsAsFactors = FALSE, header=TRUE,
                   sep="\t", check.names=FALSE, quote="", comment.char="")
expect_equal(nrow(last), 10)
expect_equal(m2$compute$steps$LD$debug$loadCounter, 1)
expect_equal(last[['example.bgen:SNP']],
             bgen[['example.bgen:SNP']][190:199])
lr <- loadResults(file.path(tdir, "out.log"), "snp2F")
expect_equal(colnames(lr), c("MxComputeLoop1", "CHR", "BP", "SNP", "A1", "A2",
                             "statusCode", "catch1", "snp2F", "Z", "P"))

# -----
# These are the same SNP after permutation:

m3 <- GWAS(buildOneFac(pheno, paste0("i", 1:numIndicators)),
           file.path(dir,"example.pgen"),
           file.path(tdir, "out.log"), SNP=1)
expect_true(is.na(m3$data$observed$snp[1]))

m3 <- GWAS(buildOneFac(pheno, paste0("i", 1:numIndicators)),
           file.path(dir,"example.bgen"),
           file.path(tdir, "out.log"), SNP=2)
expect_true(is.na(m3$data$observed$snp[1]))

# ------------------

expect_equal(pgen$A1, bgen$A1)
expect_equal(pgen$A1, bed$A1)
expect_equal(pgen$A2, bgen$A2)
expect_equal(pgen$A2, bed$A2)

# should not match exactly because bed can't store dosages
expect_equal(sum(abs(pgen$MAF - bed$MAF)), 1.447, .01)
expect_equal(pgen$MAF, bed$MAF, .06)

expect_equal(match(pgen$SNP, bed$SNP), 1:199)  # same order

bgen <- bgen[match(pgen$SNP, sub("^SNP","RS",bgen[['SNP']])),]

expect_equal(bgen$MAF, pgen$MAF, 1e-4)

mask <- (bgen$catch1=="" & pgen$catch1 == "" &
           bgen$statusCode=="OK" & pgen$statusCode=="OK")
bgen <- bgen[mask,]
pgen <- pgen[mask,]
expect_equal(sum(mask), 195)

rmse <- function(x,y) sqrt(mean((x-y)^2))
expect_equal(rmse(bgen$snp2F, pgen$snp2F), 0, tolerance=.4)

#cat(deparse(pgen$ID))
# c("RSID_4", "RSID_5", "RSID_7", "RSID_11", "RSID_12", "RSID_15",  "RSID_27", "RSID_28", "RSID_35", "RSID_37", "RSID_39", "RSID_40",  "RSID_41", "RSID_46", "RSID_47", "RSID_49", "RSID_51", "RSID_60",  "RSID_63", "RSID_65", "RSID_66", "RSID_67", "RSID_68", "RSID_69",  "RSID_71", "RSID_85", "RSID_90", "RSID_97", "RSID_98", "RSID_100",  "RSID_101", "RSID_105", "RSID_111", "RSID_112", "RSID_115", "RSID_116",  "RSID_125", "RSID_127", "RSID_128", "RSID_135", "RSID_139", "RSID_140",  "RSID_141", "RSID_147", "RSID_151", "RSID_160", "RSID_163", "RSID_166",  "RSID_167", "RSID_168", "RSID_171", "RSID_182", "RSID_185", "RSID_190",  "RSID_192", "RSID_197", "RSID_200")

