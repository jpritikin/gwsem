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

m1 <- GWAS(buildOneFac(pheno, paste0("i", 1:numIndicators)),
           file.path(dir,"example.pgen"))

pgen <- read.table("out.log", stringsAsFactors = FALSE, header=TRUE,
                     sep="\t", check.names=FALSE, quote="", comment.char="")
expect_equal(nrow(pgen), 199)
expect_equal(m1$compute$steps$LD$debug$loadCounter, 1)

m1 <- GWAS(buildOneFac(pheno, paste0("i", 1:numIndicators)),
     file.path(dir,"example.bed"))
bed <- read.table("out.log", stringsAsFactors = FALSE, header=TRUE,
                   sep="\t", check.names=FALSE, quote="", comment.char="")
expect_equal(nrow(bed), 199)
expect_equal(m1$compute$steps$LD$debug$loadCounter, 1)

m1 <- GWAS(buildOneFac(pheno, paste0("i", 1:numIndicators)),
     file.path(dir,"example.bgen"))
bgen <- read.table("out.log", stringsAsFactors = FALSE, header=TRUE,
                  sep="\t", check.names=FALSE, quote="", comment.char="")
expect_equal(nrow(bgen), 199)
expect_equal(m1$compute$steps$LD$debug$loadCounter, 1)

# ------------------

expect_equal(match(pgen$ID, bed$SNP), 1:199)  # same order

bgen <- bgen[match(pgen$ID, sub("^SNP","RS",bgen[['example.bgen:SNP']])),]

mask <- (bgen$catch1=="" & pgen$catch1 == "" & bed$catch1=="" &
           bgen$statusCode=="OK" & pgen$statusCode=="OK" & bed$statusCode=="OK" &
           !is.na(bgen$snpRegSE) & !is.na(pgen$snpRegSE) & !is.na(bed$snpRegSE))
bgen <- bgen[mask,]
bed <- bed[mask,]
pgen <- pgen[mask,]

mask <- (abs(pgen$snpReg) < 2.6*mad(pgen$snpReg) &
           (abs(pgen$snpReg / pgen$snpRegSE) > .5) &
           (abs(bgen$snpReg / bgen$snpRegSE) > .5) &
           (abs(bed$snpReg / bed$snpRegSE) > .5))
bgen <- bgen[mask,]
bed <- bed[mask,]
pgen <- pgen[mask,]

rmse <- function(x,y) sqrt(mean((x-y)^2))
expect_equal(rmse(pgen$snpReg, bed$snpReg), 0, tolerance=.2)
expect_equal(rmse(bgen$snpReg, bed$snpReg), 0, tolerance=.2)
expect_equal(rmse(bgen$snpReg, pgen$snpReg), 0, tolerance=.2)

#cat(deparse(pgen$ID))
# c("RSID_4", "RSID_5", "RSID_7", "RSID_11", "RSID_12", "RSID_15",  "RSID_27", "RSID_28", "RSID_35", "RSID_37", "RSID_39", "RSID_40",  "RSID_41", "RSID_46", "RSID_47", "RSID_49", "RSID_51", "RSID_60",  "RSID_63", "RSID_65", "RSID_66", "RSID_67", "RSID_68", "RSID_69",  "RSID_71", "RSID_85", "RSID_90", "RSID_97", "RSID_98", "RSID_100",  "RSID_101", "RSID_105", "RSID_111", "RSID_112", "RSID_115", "RSID_116",  "RSID_125", "RSID_127", "RSID_128", "RSID_135", "RSID_139", "RSID_140",  "RSID_141", "RSID_147", "RSID_151", "RSID_160", "RSID_163", "RSID_166",  "RSID_167", "RSID_168", "RSID_171", "RSID_182", "RSID_185", "RSID_190",  "RSID_192", "RSID_197", "RSID_200")