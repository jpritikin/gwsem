library(testthat)
library(gwsem)
library(digest)

skip_if(Sys.info()[["machine"]] != 'x86_64')

suppressWarnings(RNGversion("3.5"))
set.seed(2)
totalRep <- 25

dir <- system.file("extdata", package = "gwsem")
tdir <- tempdir()

test_fail <- function(m1, snpData) {
  toSkip <- as.logical(rbinom(nrow(m1$data$observed)-5, 1, .2))
  m1$data$observed <- m1$data$observed[!toSkip,]
  expect_error(GWAS(m1,
             rowFilter = list(OneItem=toSkip),
             snpData,
             file.path(tdir, "out.log"), SNP=1),
             "does not match the number of rows of observed")
}

test_subset <- function(m1, snpData, fullSnp) {
  toSkip <- as.logical(rbinom(nrow(m1$data$observed), 1, .2))
  m1$data$observed <- m1$data$observed[!toSkip,]
  m1 <- GWAS(m1,
             rowFilter = list(OneItem=toSkip),
             snpData,
             file.path(tdir, "out.log"), SNP=1)
  expect_equal(m1$data$observed$snp, fullSnp[!toSkip])
}

# ---------

pheno <- read.table(file.path(dir, "example.psam"),
                    as.is = TRUE,header=TRUE, comment.char="")

oi <- buildItem(pheno, 'phenotype')
oi$data$naAction <- "pass"
oi <- GWAS(oi,
           file.path(dir,"example.pgen"),
           file.path(tdir, "out.log"), SNP=1)
expect_equal(sha1(oi$data$observed$snp),
             "aa5cc722df62aff622850ba00cc0dd844ff8e0d6")
snp1 <- oi$data$observed$snp

test_fail(oi, file.path(dir,"example.pgen"))
for (rep in 1:totalRep) test_subset(oi, file.path(dir,"example.pgen"), snp1)

# -------

oi <- buildItem(pheno, 'phenotype')
oi$data$naAction <- "pass"
oi <- GWAS(oi,
           file.path(dir,"example.bgen"),
           file.path(tdir, "out.log"), SNP=1)
expect_equal(sha1(oi$data$observed$snp),
             "0be5f6e2ff0c9fa3df2f3a0aed494185aec15e2a")
snp2 <- oi$data$observed$snp

for (rep in 1:totalRep) test_subset(oi, file.path(dir,"example.bgen"), snp2)

# -------

pheno <- read.table(file.path(dir, "test2.psam"),
                    as.is = TRUE,header=TRUE, comment.char="")
oi <- buildItem(pheno, 'PHENO1')
oi$data$naAction <- "pass"
oi <- GWAS(oi,
           file.path(dir, "test2.pgen"),
           file.path(tdir, "out.log"), SNP=1)
expect_equal(sha1(oi$data$observed$snp),
             "98012460ea29b97caf5c256242222f0934dae750")
snp3 <- oi$data$observed$snp

for (rep in 1:totalRep) test_subset(oi, file.path(dir, "test2.pgen"), snp3)
