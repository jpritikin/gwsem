require(gwsem)
require(MASS)

## Some basic parameters for the simulation
numSNP <- 2000        # number of SNPs in the data file
numPeople <- 6000     # number of individuals in the data
numPCs <- 6           # principle components


# Simulate the SNP data
# plink2 --dummy 25000 10000 scalar-pheno dosage-freq=.9
# The output files are renamed example.*


# Simulate Phenotypic Data for a One Factor Model
# These numbers can be changed to get different parameter estimates

loadings            <- matrix(c(.7, .7, .7, 0,0,0, 0,0,0, .7, .7, .7), 6, 2)                                                            # The factor loadings
facCov              <- matrix(c(1,.3,.3,1),2,2)
res                 <- c(1-diag(loadings %*% facCov %*%  t(loadings)))                                                                     # The item residuals

sigma               <- loadings %*% facCov %*%  t(loadings) + diag(res)                                                # The implied covariance matrix
phenoData           <- as.data.frame(mvrnorm(numPeople, mu=rep(0,6), Sigma = sigma, empirical=TRUE))       # Simulate the data
colnames(phenoData) <- c('A1','A2','A3','B1',"B2","B3")                                                    # add column names



# Make GWAS Hits

m1 <- buildItem(phenoData, depVar = "A1")          # We are going to construct a GW-SEM model (this is not a real model but we can use it to grab snps from the genetic file in their proper format )

##############################################################################################################################
##############################################################################################################################
##############################################################################################################################

# Select SNPs that will predict the latent variable

# Choose K SNPs that we will make associated with the latent variable
set.seed(12345)
k <- 3                                                  # We are going to simulate 4 hits: 2 with latent A and 2 with latent B
hits <- list()                                          # Make an object to put the "hits" in
hits[['two']] <- sample.int(numSNP, k)         # We are going to sample k snps (randomly) from the number of SNPs

# Select SNPs that will predict the latent variable
snp1 <- GWAS(m1, 'example.pgen', out = 'result.log', SNP = hits[['two']][1])
snp2 <- GWAS(m1, 'example.pgen', out = 'result.log', SNP = hits[['two']][2])
snp3 <- GWAS(m1, 'example.pgen', out = 'result.log', SNP = hits[['two']][3])

# Add the effect of the SNPs to the latent variables
phenoData[,1:3] <- (phenoData[,1:3] + .2 * snp1$data$observed$snp + .16 * snp2$data$observed$snp)
phenoData[,4:6] <- (phenoData[,4:6] + .2 * snp1$data$observed$snp + .16 * snp3$data$observed$snp)

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


write.table(phenoData, "phenoTwoData.txt", col.names = T, quote = F, row.names = F)




