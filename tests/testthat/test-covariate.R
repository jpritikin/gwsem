library(testthat)
library(gwsem)
library(MASS)

skip_if(Sys.info()[["machine"]] != 'x86_64')

suppressWarnings(RNGversion("3.5"))
set.seed(1)

dir <- system.file("extdata", package = "gwsem")
tdir <- tempdir()

pheno <- read.table(file.path(dir, "example.psam"),
                    as.is = TRUE, header=TRUE, comment.char="")
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

numCovariate <- 2
for (cx in 1:numCovariate) {
  pheno[[paste0("covar", cx)]] <- rnorm(nrow(pheno))
}

# -----

m1 <- buildOneFac(pheno, paste0("i", 1:numIndicators),
                  covariates = paste0("covar",1:numCovariate), exogenous=TRUE)
expect_equivalent(m1$M$labels[1,'covar1'], 'data.covar1')
GWAS(m1, file.path(dir,"example.pgen"),
     file.path(tdir,"out.log"))

pgen <- read.delim(file.path(tdir,"out.log"), as.is = TRUE,
                   check.names=FALSE, quote="", comment.char="")

mask <- (pgen$catch1 == "" & pgen$statusCode=="OK" & !is.na(pgen[['Vsnp_to_F:snp_to_F']]))
pgen <- pgen[mask,]

mask <- (abs(pgen$snp_to_F) < 2.6*mad(pgen$snp_to_F) & (abs(pgen$snp_to_F / sqrt(pgen[['Vsnp_to_F:snp_to_F']])) > .5))
pgen <- pgen[mask,]

cvNames <- paste(rep(paste0("covar",1:numCovariate), each = numIndicators),
      paste0("i", 1:numIndicators), sep = "_to_")
expect_equivalent(colMeans(pgen[,cvNames]), rep(0, length(cvNames)), tolerance=.1)

pgen <- loadResults(file.path(tdir,"out.log"), c("snp_to_F", 'lambda_i1'))
pgen <- signif(pgen, "snp_to_F", signAdj='lambda_i1')
expect_equal(nrow(pgen), 197, 2)
expect_error(plot(pgen, y=1),
             "plot does not accept a y= argument")

bad <- sum(isSuspicious(pgen))
expect_equal(bad, 2, 1)
expect_true(any(grepl("observed variance less",
                      pgen[isSuspicious(pgen),'catch1'], fixed=TRUE)))
pgen <- pgen[!isSuspicious(pgen),]

# -----

m2 <- buildOneFac(pheno, paste0("i", 1:numIndicators),
                  covariates = paste0("covar",1:numCovariate),
                  exogenous = FALSE)
expect_equivalent(m2$M$labels[1,'covar1'], 'covar1_mean')
GWAS(m2,
     file.path(dir,"example.pgen"),
     file.path(tdir,"out2.log"))

pgen2 <- loadResults(file.path(tdir,"out2.log"), c("snp_to_F", 'lambda_i1'))
pgen2 <- signif(pgen2, "snp_to_F", signAdj='lambda_i1')
pgen2 <- pgen2[!isSuspicious(pgen2),]
expect_equal(nrow(pgen2), 196, 2)

both <- intersect(pgen$SNP, pgen2$SNP)
expect_equal(length(both), 196, 2)
pgen <- subset(pgen, SNP %in% both)
pgen2 <- subset(pgen2, SNP %in% both)

mask <- abs(pgen$Z) < 4 & abs(pgen2$Z) < 4

expect_equivalent(cor(pgen[mask,'Z'], pgen2[mask,'Z']), 1, tolerance=.05)

# ----- compare OneFacRes exo vs endo covariates

m1 <- buildOneFacRes(pheno, paste0("i", 1:numIndicators),
                  covariates = paste0("covar",1:numCovariate), exogenous=TRUE)
expect_equivalent(m1$M$labels[1,'covar1'], 'data.covar1')
GWAS(m1, file.path(dir,"example.pgen"),
     file.path(tdir,"outx.log"), SNP=1:100)

m2 <- buildOneFacRes(pheno, paste0("i", 1:numIndicators),
                     covariates = paste0("covar",1:numCovariate),
                     exogenous=FALSE)
expect_equivalent(m2$M$labels[1,'covar1'], 'covar1_mean')
GWAS(m2, file.path(dir,"example.pgen"),
     file.path(tdir,"out.log"), SNP=1:100)

for (ind in paste0("snp_to_i", 1:numIndicators)) {
  m1o <- loadResults(file.path(tdir,"outx.log"), ind)
  m1o <- signif(m1o, ind)
  m1o <- m1o[!isSuspicious(m1o),]
  m2o <- loadResults(file.path(tdir,"out.log"), ind)
  m2o <- signif(m2o, ind)
  m2o <- m2o[!isSuspicious(m2o),]
  both <- intersect(m1o$SNP, m2o$SNP)
  expect_equal(length(both), 200, 20)
  m1o <- subset(m1o, SNP %in% both)
  m2o <- subset(m2o, SNP %in% both)
  mask <- abs(m1o$Z) < 4 & abs(m2o$Z) < 4
  expect_equivalent(cor(m1o[mask,'Z'], m2o[mask,'Z']), 1, tolerance=.4)
}

# ----- compare TwoFac exo vs endo covariates

m1 <- buildTwoFac(pheno,
                  paste0("i", 1:(numIndicators-1)),
                  paste0("i", 2:numIndicators),
                  covariates = paste0("covar",1:numCovariate),
                  exogenous=TRUE)
expect_equivalent(m1$M$labels[1,'covar1'], 'data.covar1')
GWAS(m1, file.path(dir,"example.pgen"),
     file.path(tdir,"outx.log"), SNP=1:50)

m2 <- buildTwoFac(pheno,
                  paste0("i", 1:(numIndicators-1)),
                  paste0("i", 2:numIndicators),
                  covariates = paste0("covar",1:numCovariate),
                  exogenous=FALSE)
expect_equivalent(m2$M$labels[1,'covar1'], 'covar1_mean')
GWAS(m2, file.path(dir,"example.pgen"),
     file.path(tdir,"out.log"), SNP=1:50)

for (fx in 1:2) {
  ind <- paste0("snp_to_F", fx)
  sa <- paste0('F',fx,'_lambda_i2')
  m1o <- loadResults(file.path(tdir,"outx.log"), c(ind, sa))
  m1o <- signif(m1o, ind, signAdj=sa)
  m1o <- m1o[!isSuspicious(m1o),]
  m2o <- loadResults(file.path(tdir,"out.log"), c(ind, sa))
  m2o <- signif(m2o, ind, signAdj=sa)
  m2o <- m2o[!isSuspicious(m2o),]
  both <- intersect(m1o$SNP, m2o$SNP)
  expect_equal(length(both), 50, 5)
  m1o <- subset(m1o, SNP %in% both)
  m2o <- subset(m2o, SNP %in% both)
  mask <- abs(m1o$Z) < 4 & abs(m2o$Z) < 4
  # Wow, this tolerance is really terrible TODO
  expect_equivalent(cor(m1o[mask,'Z'], m2o[mask,'Z']), 1, tolerance=.9)
}

# ----- test continuous endogenous covariates

lastFit <- GWAS(buildTwoFac(pheno,
                            paste0("i", 1:(numIndicators-1)),
                            paste0("i", 2:numIndicators),
                 covariates = paste0("covar",1:numCovariate),
                 exogenous = FALSE), SNP=1,
     file.path(dir,"example.pgen"),
     file.path(tdir,"out.log"))

expect_equal(lastFit$A$labels['F1','covar1'], "covar1_to_F1")
expect_equal(lastFit$A$labels['i2','F1'], "F1_lambda_i2")

# ----- test ordinal endogenous covariates

lastFit <- GWAS(buildOneFac(pheno, paste0("i", 2:numIndicators),
                 covariates = 'i1', exogenous = FALSE), SNP=1,
     file.path(dir,"example.pgen"),
     file.path(tdir,"out.log"))

expect_true(!lastFit$M$free[,'i1'])
expect_equal(lastFit$S$values['i1','i1'], 1)
expect_equal(lastFit$A$labels['F','i1'], 'i1_to_F')
expect_equal(lastFit$A$values['F','i1'], .66, tolerance=.01)

# -----

m1 <- buildItem(pheno, "phenotype", covariates=paste0('i', 3:5))
expect_equal(length(coef(m1)), 9)

