library(gwsem)

# Read the phenotype data into R and look at the data
phenoData <- read.table("oneFacphenoData.txt", header=TRUE)
head(phenoData)

# This would be a great time to recode the data, transform it, and tell R if we have ordinal or binary indicators using mxFactor()
# The data for this example are simulated and continuous, and therefore, we will not be doing anything now, but if you would like 
# to chop the indicators up into binary or ordinal variables, this would be a great time to do it.


# The first thing that we would want to do is build a general model. 
# This model takes a series of arguments: 

                                                                                     # You must tell GW-SEM:
addFac <- buildOneFac(phenoData = phenoData,                                         # what the data object is (which you read in above)
                     itemNames = c("tobacco", "cannabis", "alcohol"),                # what the items of the latent factor are
                     covariates=c('pc1','pc2','pc3','pc4','pc5'),                    # what covariates that you want to include in the analysis
                     fitfun = "WLS")                                                 # and the fit function that you would like to use (WLS is much faster than ML)

# You can take this object (addFac) which is technically an OpenMx model, and simply fit it using mxRun() 
addFacFit <- mxRun(addFac)
summary(addFacFit)

# This is strongly advised, as it is a great time to learn if your model is being specified in a way that you didn't expect or that 
# the model is not giving you the answers that you are expecting

# Now, provided that the model looks reasonable, you can plug the model that you have built into the GWAS function
                                                                                     # To fit this function, you must tell GW-SEM :
GWAS(model = addFac,                                                                 # what model object you would like to fit
	snpData = 'example.pgen',                                                        # that path to the snpData file.
	out="latFac.log",                                                                # the file that you would like to save the full results into
	SNP=1:200)                                                                       # the index of the snps (how many) you would like to fit

# Note about snpData: The path to your snpData will likely include switching to a different directory (as you will likely do your analysis in a different
# folder than your SNP data). All you need to do is point to the data using relative paths. Further, it is able to take plink bed/bim/fam or pgen/psam/pvar 
# data or bgen data (oxford format)

# Note about SNP: This can be used to run a limited number of SNP (i.e. not the whole snp file). This is particularly useful if you would like to run chop up a 
# chr into several parts without cutting you actual genotype data into seperate files.

# This will take a while and frequently be done on a cluster.  If you chop up this script and remove the comments and insert you data, you can use it on your cluster


# The next step is to read the data in. While you are unlikely to want to read all of the output into R, it is useful to do this on a small test GWAS analysis
# to ensure that you are getting reasonable results.

FullResult <- read.delim("latFac.log")                            # This is didactically useful, but it contains much more information than most people want


# Then, we can read the results into R for a specific parameter. This function takes takes two arguments: the path to the data and the column in the results 
# file for the parameter that you want to examine.

succinct <- loadResults(path = "latFac.log", focus =  "snp_to_F")

# I will give you a way to read several files into a single R object 

# Lastly, you can plot the results
plot(succinct)


### The Residuals Model

# The process is very similar for the Latent factor model except that we are regressing the individual items on the SNP

# Again, the first thing that we would want to do is build a general model. 
# The function tells GW-SEM that you want to run the residuals model. This model takes a series of arguments: 
                                                                                       # You must tell GW-SEM:
addFacRes <- buildOneFacRes(phenoData,                                                 # what the data object is (which you read in above)
                     c("tobacco", "cannabis", "alcohol"),                              # what the items of the latent factor are
                     covariates=c('pc1','pc2','pc3','pc4','pc5'),                      # what covariates that you want to include in the analysis
                     fitfun = "WLS")                                                   # and the fit function that you would like to use (WLS is much faster than ML)


# Again, You can take this object (addFac) which is technically an OpenMx model, and simply fit it using mxRun() 
addFacResFit <- mxRun(addFacRes)
summary(addFacResFit)


# Then you will use the same function that we used for the latent variable model, with a different model object

GWAS(model = addFacRes, snpData = 'example.pgen', out="facRes.log")                 

# It is instructive to look at the results again, but this might be easier to do in linux rather than R
# This is how you load it into R
ResResult <- read.delim("facRes.log")                            


# Now that we have run the residual's model, we have multiple parameters that we want to load into R (Three in this case)
# The function to do this is the same as for one parameter, but you need to do it several times (once for each paramter). This gives several R objects.

succinctTob <- loadResults(path = "facRes.log", focus =  "snp_to_tobacco")
succinctCan <- loadResults(path = "facRes.log", focus =  "snp_to_cannabis")
succinctAlc <- loadResults(path = "facRes.log", focus =  "snp_to_alcohol")

# Now we can plot all the residual manhattan plots

plot(succinctTob)
plot(succinctCan)
plot(succinctAlc)


# Additional considerations: 

# Because we are working with numerous chromosomes, which we have chopped up into different chunks, we have several log files for each chromosome. For example, 
# we may have a set of log files such as: fac1a.log, fac1b.log, fac1c.log, fac1d.log, fac1e.log, fac1f.log, fac1g.log, fac1h.log.
# To simplify things, we can construct a list of file names for each chr with the following code.

c1  <- paste0(paste0("fac1",  letters)[1:8], ".log")
c2  <- paste0(paste0("fac2",  letters)[1:9], ".log")
c3  <- paste0(paste0("fac3",  letters)[1:7], ".log")
c4  <- paste0(paste0("fac4",  letters)[1:8], ".log")
c5  <- paste0(paste0("fac5",  letters)[1:7], ".log")
c6  <- paste0(paste0("fac6",  letters)[1:7], ".log")
c7  <- paste0(paste0("fac7",  letters)[1:6], ".log")
c8  <- paste0(paste0("fac8",  letters)[1:6], ".log")
c9  <- paste0(paste0("fac9",  letters)[1:5], ".log")
c10 <- paste0(paste0("fac10", letters)[1:5], ".log") 
c11 <- paste0(paste0("fac11", letters)[1:5], ".log")
c12 <- paste0(paste0("fac12", letters)[1:5], ".log")
c13 <- paste0(paste0("fac13", letters)[1:4], ".log")
c14 <- paste0(paste0("fac14", letters)[1:4], ".log")
c15 <- paste0(paste0("fac15", letters)[1:3], ".log")
c16 <- paste0(paste0("fac16", letters)[1:3], ".log")
c17 <- paste0(paste0("fac17", letters)[1:3], ".log")
c18 <- paste0(paste0("fac18", letters)[1:3], ".log")
c19 <- paste0(paste0("fac19", letters)[1:3], ".log")
c20 <- paste0(paste0("fac20", letters)[1:2], ".log")
c21 <- paste0(paste0("fac21", letters)[1:2], ".log")
c22 <- paste0(paste0("fac22", letters)[1:2], ".log")

# Then we can plug these files into the loadResults function so that all the data is loaded into a single object.

res <- loadResults(c(c1, c2, c3, c4, c5, c6, c7, c8, c9, c10,
	                 c11, c12, c13, c14, c15, c16, c17, c18, c19, c20,
					 c21, c22), "snp_to_F")



