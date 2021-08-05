library(testthat)
library(gwsem)

skip_if(Sys.info()[["machine"]] != 'x86_64')

suppressWarnings(RNGversion("3.5"))
set.seed(1)

dir <- system.file("extdata", package = "gwsem")
tdir <- tempdir()

pheno1 <- read.table(file.path(dir, "group1.fam"),
                     as.is = TRUE, comment.char="")

g1 <- mxRename(buildItem(pheno1, depVar = 'V6'), "sample1")
g1 <- omxSetParameters(g1, labels = names(coef(g1)),
                newlabels = paste0('S1',names(coef(g1))))

f1 <- GWAS(g1, file.path(dir, "group1.bed"),
           out=file.path(tdir, "out.log"), SNP=1)

pheno2 <- read.table(file.path(dir, "group2.fam"),
                     as.is = TRUE, comment.char="")

g2 <- mxRename(buildItem(pheno2, depVar = 'V6'), "sample2")
g2 <- omxSetParameters(g2, labels = names(coef(g2)),
                       newlabels = paste0('S2',names(coef(g2))))

f2 <- GWAS(g2, file.path(dir, "group2.bed"),
           out=file.path(tdir, "out.log"), SNP=1)

mg <- mxModel("mg", g1, g2, mxFitFunctionMultigroup(paste0("sample", 1:2)))

f3 <- GWAS(mg, c('sample1'=file.path(dir, "group1.bed"),
                 'sample2'=file.path(dir, "group2.bed")),
           out=file.path(tdir, "out.log"), SNP=1)

expect_equal(c(coef(f1), coef(f2)), coef(f3), 1e-4)
