require(gwsem)
require(MASS)

## Some basic parameters for the simulation
numSNP <- 2000        # number of SNPs in the data file
numPeople <- 6000     # number of individuals in the data
numPCs <- 6           # principle components


# Simulate the SNP data
# ./plink2 --dummy 6000 2000 scalar-pheno dosage-freq=.9
# The output files are renamed example.*


# Simulate Phenotypic Data for a One Factor Model
# These numbers can be changed to get different parameter estimates
loadings            <- as.matrix(c(.7, .7, .7))                                                            # The factor loadings
res                 <- c(1-loadings^2)                                                                     # The item residuals
sigma               <- loadings %*% t(loadings) + diag(res)                                                # The implied covariance matrix
phenoData           <- as.data.frame(mvrnorm(numPeople, mu=rep(0,3), Sigma = sigma, empirical=TRUE))       # Simulate the data
colnames(phenoData) <- c('tobacco','cannabis','alcohol')                                                   # add column names



# Make GWAS Hits

m1 <- buildItem(phenoData, depVar = "tobacco")          # We are going to construct a GW-SEM model (this is not a real model but we can use it to grab snps from the genetic file in their proper format )

##############################################################################################################################
##############################################################################################################################
##############################################################################################################################

# Select SNPs that will predict the latent variable

# Choose K SNPs that we will make associated with the latent variable

k <- 2                                                  # We are going to simulate 2 hits with the phenotype
hits <- list()                                          # Make an object to put the "hits" in
hits[['substanceUse']] <- sample.int(numSNP, k)         # We are going to sample k snps (randomly) from the number of SNPs

# Select SNPs that will predict the latent variable
snp1 <- GWAS(m1, 'example.pgen', out = 'result.log', SNP = hits[['substanceUse']][1])
snp2 <- GWAS(m1, 'example.pgen', out = 'result.log', SNP = hits[['substanceUse']][2])

# Add the effect of the SNPs to the latent variable
phenoData <- (phenoData + .1 * snp1$data$observed$snp + .08 * snp2$data$observed$snp)

##############################################################################################################################
##############################################################################################################################
##############################################################################################################################

# Select SNPs that will predict the tobacco residual

# Choose K SNPs that we will make associated with the tobacco residual

k_tob <- 2                                                  # We are going to simulate 2 hits with the phenotype
hits_tob <- list()                                          # Make an object to put the "hits" in
hits_tob[['tob']] <- sample.int(numSNP, k_tob)              # We are going to sample k snps (randomly) from the number of SNPs

# Select SNPs that will predict the latent variable
snp3 <- GWAS(m1, 'example.pgen', out = 'result.log', SNP = hits_tob[['tob']][1])
snp4 <- GWAS(m1, 'example.pgen', out = 'result.log', SNP = hits_tob[['tob']][2])

# Add the effect of the SNPs to the latent variable
phenoData <- (phenoData + .1 * snp3$data$observed$snp + .08 * snp4$data$observed$snp)


##############################################################################################################################
##############################################################################################################################
##############################################################################################################################

# Select SNPs that will predict the cannabis residual

# Choose K SNPs that we will make associated with the cannabis residual

k_can <- 2                                                  # We are going to simulate 2 hits with the phenotype
hits_can <- list()                                          # Make an object to put the "hits" in
hits_can[['can']] <- sample.int(numSNP, k_can)              # We are going to sample k snps (randomly) from the number of SNPs

# Select SNPs that will predict the latent variable
snp5 <- GWAS(m1, 'example.pgen', out = 'result.log', SNP = hits_can[['can']][1])
snp6 <- GWAS(m1, 'example.pgen', out = 'result.log', SNP = hits_can[['can']][2])

# Add the effect of the SNPs to the latent variable
phenoData <- (phenoData + .1 * snp5$data$observed$snp + .08 * snp6$data$observed$snp)


##############################################################################################################################
##############################################################################################################################
##############################################################################################################################

# Select SNPs that will predict the cannabis residual

# Choose K SNPs that we will make associated with the cannabis residual

k_alc <- 2                                                  # We are going to simulate 2 hits with the phenotype
hits_alc <- list()                                          # Make an object to put the "hits" in
hits_alc[['alc']] <- sample.int(numSNP, k_alc)              # We are going to sample k snps (randomly) from the number of SNPs

# Select SNPs that will predict the latent variable
snp7 <- GWAS(m1, 'example.pgen', out = 'result.log', SNP = hits_alc[['alc']][1])
snp8 <- GWAS(m1, 'example.pgen', out = 'result.log', SNP = hits_alc[['alc']][2])

# Add the effect of the SNPs to the latent variable
phenoData <- (phenoData + .1 * snp7$data$observed$snp + .08 * snp8$data$observed$snp)


##############################################################################################################################
##############################################################################################################################
##############################################################################################################################

# Generate principle components (general covariates)

phenoData$pc1 <- rnorm(numPeople)
phenoData$pc2 <- rnorm(numPeople)
phenoData$pc3 <- rnorm(numPeople)
phenoData$pc4 <- rnorm(numPeople)
phenoData$pc5 <- rnorm(numPeople)

# Take a look at the data and see if we have built a data frame properly
head(phenoData)


write.table(phenoData, "oneFacPhenoData.txt", col.names = T, quote = F, row.names = F)

