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
sigma               <- matrix(c(1, .3, .4, .3,1, .2, .4,.2,1 ), 3, 3)                                           # The implied covariance matrix
phenoData           <- as.data.frame(mvrnorm(numPeople, mu=rep(0,3), Sigma = sigma, empirical=TRUE))       # Simulate the data
colnames(phenoData) <- c('tobacco','cannabis','alcohol')                                                   # add column names


# Make GWAS Hits

m1 <- buildItem(phenoData, depVar = "tobacco")          # We are going to construct a GW-SEM model (this is not a real model but we can use it to grab snps from the genetic file in their proper format )

##############################################################################################################################
##############################################################################################################################
##############################################################################################################################

# Select SNPs that will predict the latent variable

# Choose K SNPs that we will make associated with the latent variable

k <- 6                                                  # We are going to simulate 2 hits with the phenotype
hits <- list()                                          # Make an object to put the "hits" in
hits[['items']] <- sample.int(numSNP, k)         # We are going to sample k snps (randomly) from the number of SNPs

# Select SNPs that will predict the latent variable
snp1 <- GWAS(m1, 'example.pgen', out = 'result.log', SNP = hits[['items']][1])
snp2 <- GWAS(m1, 'example.pgen', out = 'result.log', SNP = hits[['items']][2])
snp3 <- GWAS(m1, 'example.pgen', out = 'result.log', SNP = hits[['items']][3])
snp4 <- GWAS(m1, 'example.pgen', out = 'result.log', SNP = hits[['items']][4])
snp5 <- GWAS(m1, 'example.pgen', out = 'result.log', SNP = hits[['items']][5])
snp6 <- GWAS(m1, 'example.pgen', out = 'result.log', SNP = hits[['items']][6])



# Add the effect of the SNPs to the items
phenoData$tobacco   <- phenoData$tobacco   + .17 * snp1$data$observed$snp + .15 * snp2$data$observed$snp
phenoData$cannabis  <- phenoData$cannabis  + .15 * snp3$data$observed$snp + .18 * snp4$data$observed$snp
phenoData$alcohol   <- phenoData$alcohol   + .17 * snp5$data$observed$snp + .17 * snp6$data$observed$snp


# Chop continuous variables into ordinal data 
# with nThresholds+1 approximately equal categories, based on 1st variable
q1<-quantile(phenoData$tobacco,  probs = c((1:2)/(2+1)))
q2<-quantile(phenoData$cannabis,  probs = c((1:3)/(3+1)))
q3<-quantile(phenoData$alcohol,  probs = c((1:4)/(4+1)))

phenoData$tobacco  <- cut(as.vector(phenoData$tobacco ),c(-Inf,q1,Inf), labels = 0:2)
phenoData$cannabis <- cut(as.vector(phenoData$cannabis),c(-Inf,q2,Inf), labels = 0:3)
phenoData$alcohol  <- cut(as.vector(phenoData$alcohol ),c(-Inf,q3,Inf), labels = 0:4)

table(phenoData$tobacco )
table(phenoData$cannabis)
table(phenoData$alcohol )

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

summary(polr(phenoData$tobacco  ~ snp1$data$observed$snp + snp2$data$observed$snp ))
summary(polr(phenoData$cannabis ~ snp3$data$observed$snp + snp4$data$observed$snp ))
summary(polr(phenoData$alcohol  ~ snp5$data$observed$snp + snp6$data$observed$snp ))



write.table(phenoData, "itemPhenoData.txt", col.names = T, quote = F, row.names = F)


