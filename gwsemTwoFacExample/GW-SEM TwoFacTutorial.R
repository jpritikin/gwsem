library(gwsem)

# Read the phenotype data into R and look at the data
TwoDat <- read.table("phenoTwoData.txt", header=TRUE)
head(TwoDat)

# This would be a great time to recode the data, transform it, and tell R if we have ordinal or binary indicators using mxFactor()
# The data for this example are simulated and continuous, and therefore, we will not be doing anything now, but if you would like 
# to chop the indicators up into binary or ordinal variables, this would be a great time to do it.


# The first thing that we would want to do is build a general model. 
# This model takes a series of arguments: 

                                                                                     # You must tell GW-SEM:
twoFac <- buildTwoFac(phenoData = TwoDat,                                         # what the data object is (which you read in above)
                     F1itemNames = c("A1", "A2", "A3"),                # what the items of the latent factor are
                     F2itemNames = c("B1", "B2", "B3"),                # what the items of the latent factor are
                     covariates=c('pc1','pc2','pc3','pc4','pc5'),                    # what covariates that you want to include in the analysis
                     fitfun = "WLS",
					 exogenous = T)                                                 # and the fit function that you would like to use (WLS is much faster than ML)

# You can take this object (addFac) which is technically an OpenMx model, and simply fit it using mxRun() 
twoFacFit <- mxRun(twoFac)
summary(twoFacFit)

# This is strongly advised, as it is a great time to learn if your model is being specified in a way that you didn't expect or that 
# the model is not giving you the answers that you are expecting

# Now, provided that the model looks reasonable, you can plug the model that you have built into the GWAS function
                                                                                     # To fit this function, you must tell GW-SEM :
GWAS(model = twoFac,                                                                 # what model object you would like to fit
	snpData = 'example.pgen',                                                        # that path to the snpData file.
	out="twoFac.log")                                                                # the file that you would like to save the full results into
	                                                                       # the index of the snps (how many) you would like to fit


# Note about snpData: The path to your snpData will likely include switching to a different directory (as you will likely do your analysis in a different
# folder than your SNP data). All you need to do is point to the data using relative paths. Further, it is able to take plink bed/bim/fam or pgen/psam/pvar 
# data or bgen data (oxford format)

# Note about SNP: This can be used to run a limited number of SNP (i.e. not the whole snp file). This is particularly useful if you would like to run chop up a 
# chr into several parts without cutting you actual genotype data into seperate files.

# This will take a while and frequently be done on a cluster.  If you chop up this script and remove the comments and insert you data, you can use it on your cluster


# The next step is to read the data in. While you are unlikely to want to read all of the output into R, it is useful to do this on a small test GWAS analysis
# to ensure that you are getting reasonable results.

TwoResult <- read.delim("twoFac.log")                           # This is didactically useful, but it contains much more information than most people want


# Then, we can read the results into R for a specific parameter. This function takes takes two arguments: the path to the data and the column in the results 
# file for the parameter that you want to examine.

# For the two factor model we want to look at both regression paths as well as the covariance between the latent factors for jointly significant SNPs

succinct1 <- loadResults(path = "twoFac.log", focus =  "snp_to_F1")
succinct2 <- loadResults(path = "twoFac.log", focus =  "snp_to_F2")


# Now you can plot the results for each factor.
plot(succinct1)
plot(succinct2)


# Let's sort the results by the largest Z value and then examine whether there are any overlaping SNPs

head(succinct1[order(succinct1$Z, decreasing = T),])
head(succinct2[order(succinct2$Z, decreasing = T),])




# Let's look at the correlation between the factors for the SNPs that are jointly associated with both latent factors.

facCov <- as.data.frame(loadResults(path = "twoFac.log", focus =  "facCov", .retainSE = T))


facCov[facCov$SNP == "snp141",]

# And now for snps associately with only one latent factor
facCov[facCov$SNP == "snp50",]
facCov[facCov$SNP == "snp719",]

# And not for non significant SNPs (say the first 6 snps analyzed)
head(facCov)

