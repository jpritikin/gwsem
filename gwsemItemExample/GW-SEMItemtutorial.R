library(gwsem)

# Read the phenotype data into R and look at the data
phenoData <- read.table("itemPhenoData.txt", header=TRUE)
head(phenoData)

table(phenoData$tobacco)
table(phenoData$cannabis)
table(phenoData$alcohol)


# This is when we recode the data, transform it, and tell R that we have ordinal or binary indicators using mxFactor()

phenoData$tobacco  <- mxFactor(phenoData$tobacco  , levels = 0:2)
phenoData$cannabis <- mxFactor(phenoData$cannabis , levels = 0:3)
phenoData$alcohol  <- mxFactor(phenoData$alcohol  , levels = 0:4)



# First let's look at how to do a standard GWAS with only one item.

# The first step is to to build a general model. 
# This model takes a series of arguments: 

                                                                                     # You must tell GW-SEM:
tob <- buildItem(phenoData = phenoData,                                         # what the data object is (which you read in above)
                     depVar = c("tobacco"),                # what the items of the latent factor are
                     covariates=c('pc1','pc2','pc3','pc4','pc5'),                    # what covariates that you want to include in the analysis
                     fitfun = "WLS")                                                 # and the fit function that you would like to use (WLS is much faster than ML)

# You can take this object (addFac) which is technically an OpenMx model, and simply fit it using mxRun() 
tobFit <- mxRun(tob)
summary(tobFit)

# This is strongly advised, as it is a great time to learn if your model is being specified in a way that you didn't expect or that 
# the model is not giving you the answers that you are expecting

# Now, provided that the model looks reasonable, you can plug the model that you have built into the GWAS function
                                                                                     # To fit this function, you must tell GW-SEM :
GWAS(model = tob,                                                                 # what model object you would like to fit
	snpData = 'example.pgen',                                                        # that path to the snpData file.
	out="tob.log")                                                                # the file that you would like to save the full results into





# The next step is to read the data in. While you are unlikely to want to read all of the output into R, it is useful to do this on a small test GWAS analysis
# to ensure that you are getting reasonable results.

TobFullResult <- read.delim("tob.log")                            # This is didactically useful, but it contains much more information than most people want


# Then, we can read the results into R for a specific parameter. This function takes takes two arguments: the path to the data and the column in the results 
# file for the parameter that you want to examine.

succinct <- loadResults(path = "tob.log", focus =  "snp_to_tobacco")

# I will give you a way to read several files into a single R object 

# Lastly, you can plot the results
plot(succinct)


### The Multivariate GWAS Model

# The process is very similar for the standard GWAS model except that we are regressing multiple items on each SNP

# Again, the first thing that we would want to do is build a general model. 
# The function tells GW-SEM that you want to run the residuals model. This model takes a series of arguments: 
                                                                                       # You must tell GW-SEM:
multi   <- buildItem(phenoData,                                                 # what the data object is (which you read in above)
              depVar = c("tobacco", "cannabis", "alcohol"),                     # what the items of the latent factor are
              covariates=c('pc1','pc2','pc3','pc4','pc5'),                      # what covariates that you want to include in the analysis
              fitfun = "WLS")                                                   # and the fit function that you would like to use (WLS is much faster than ML)




# Again, You can take this object (multi) which is technically an OpenMx model, and simply fit it using mxRun() 
multiFit <- mxRun(multi)
summary(multiFit)


# Then you will use the same function that we used for the latent variable model, with a different model object

GWAS(model = multi, snpData = 'example.pgen', out="multi.log")                 

# It is instructive to look at the results again, but this might be easier to do in linux rather than R
# This is how you load it into R
multiResult <- read.delim("multi.log")                            


# Now that we have run the residual's model, we have multiple parameters that we want to load into R (Three in this case)
# The function to do this is the same as for one parameter, but you need to do it several times (once for each paramter). This gives several R objects.

succinctTob <- loadResults(path = "multi.log", focus =  "snp_to_tobacco")
succinctCan <- loadResults(path = "multi.log", focus =  "snp_to_cannabis")
succinctAlc <- loadResults(path = "multi.log", focus =  "snp_to_alcohol")

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



