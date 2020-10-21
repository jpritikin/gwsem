library(gwsem)
library(MASS)

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

oi <- buildItem(pheno, paste0("i", 3))

got6 <- GWAS(oi,
     file.path(dir,"example.pgen"),
     file.path(tdir, "out1.log"),
     SNP=rep(1:199, 50))

# -----------------

oi <- buildOneFac(pheno[,paste0("i", 3:5)], paste0("i", 3:5))

z1 = GWAS(oi,
     file.path(dir,"example.pgen"),
     file.path(tdir, "out2.log"),
     SNP=rep(1:199, 50))
