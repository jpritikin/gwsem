require(gwsem)
require(MASS)

set.seed(12345)

## Some basic parameters for the simulation
numSNP <- 2000        # number of SNPs in the data file
numPeople <- 6000     # number of individuals in the data
numPCs <- 6           # principle components


# Simulate the SNP data
# ./plink2 --dummy 6000 2000 scalar-pheno dosage-freq=.9
# The output files are renamed example.*


# Simulate Phenotypic Data for a One Factor Model
# These numbers can be changed to get different parameter estimates

phe <- rnorm(numPeople )       # Simulate the phenotype data
mod <- rnorm(numPeople )       
pc1 <- rnorm(numPeople)
pc2 <- rnorm(numPeople)
pc3 <- rnorm(numPeople)
pc4 <- rnorm(numPeople)
pc5 <- rnorm(numPeople)

phenoData <- as.data.frame(cbind(phe, mod, pc1, pc2, pc3, pc4, pc5))


# Take a look at the data and see if we have built a data frame properly
head(phenoData)


# Make GWAS Hits

m1 <- buildItem(phenoData, depVar = "phe")          # We are going to construct a GW-SEM model (this is not a real model but we can use it to grab snps from the genetic file in their proper format )

##############################################################################################################################
##############################################################################################################################
##############################################################################################################################

# Select SNPs that will predict the latent variable
# Choose K SNPs that we will make associated with the latent variable

k <- 6                                                  # We are going to simulate 2 hits with the phenotype
hits <- list()                                          # Make an object to put the "hits" in
hits[['gxe']] <- sample.int(numSNP, k)                  # We are going to sample k snps (randomly) from the number of SNPs

# Select SNPs that will predict the latent variable
snp1 <- GWAS(m1, 'example.pgen', out = 'result.log', SNP = hits[['gxe']][1])
snp2 <- GWAS(m1, 'example.pgen', out = 'result.log', SNP = hits[['gxe']][2])
snp3 <- GWAS(m1, 'example.pgen', out = 'result.log', SNP = hits[['gxe']][3])
snp4 <- GWAS(m1, 'example.pgen', out = 'result.log', SNP = hits[['gxe']][4])
snp5 <- GWAS(m1, 'example.pgen', out = 'result.log', SNP = hits[['gxe']][5])
snp6 <- GWAS(m1, 'example.pgen', out = 'result.log', SNP = hits[['gxe']][6])

# Add the effect of the SNPs to the latent variable
phenoData$phe <- phenoData$phe + 
	             .11 * snp1$data$observed$snp + .12 * snp2$data$observed$snp  +                                   # Main Effect Only
				 .15 * snp3$data$observed$snp + .14 * snp3$data$observed$snp * phenoData$mod   +                  # Main Effect and Interaction
				 .14 * snp4$data$observed$snp + .12 * snp4$data$observed$snp * phenoData$mod   +                  # Main Effect and Interaction
                 .14 * snp5$data$observed$snp * phenoData$mod + .18 * snp6$data$observed$snp * phenoData$mod      # Interaction Only

				 
				 
summary(lm(phenoData$phe ~ snp1$data$observed$snp + phenoData$mod + snp1$data$observed$snp:phenoData$mod))
summary(lm(phenoData$phe ~ snp2$data$observed$snp + phenoData$mod + snp2$data$observed$snp:phenoData$mod))
summary(lm(phenoData$phe ~ snp3$data$observed$snp + phenoData$mod + snp3$data$observed$snp:phenoData$mod))
summary(lm(phenoData$phe ~ snp4$data$observed$snp + phenoData$mod + snp4$data$observed$snp:phenoData$mod))
summary(lm(phenoData$phe ~ snp5$data$observed$snp + phenoData$mod + snp5$data$observed$snp:phenoData$mod))
summary(lm(phenoData$phe ~ snp6$data$observed$snp + phenoData$mod + snp6$data$observed$snp:phenoData$mod))

write.table(phenoData, "gxeData.txt", col.names = T, quote = F, row.names = F)


