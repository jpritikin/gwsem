require(gwsem)

# Read the phenotype data into R and look at the data
gxeData <- read.table("gxeData.txt", header=TRUE)
head(gxeData)

# The first thing that we would want to do is build a general model. 
# This model takes a series of arguments: 


gxe <- buildItem(phenoData = gxeData, depVar  = "phe",
                  covariates = c("mod", "snp_mod" ,"pc1", "pc2", "pc3", "pc4","pc5"),
                  fitfun = "WLS", exogenous = T, gxe = "mod")


# You can take this object (gxe) which is technically an OpenMx model, and simply fit it using mxRun() 

# This is strongly advised, as it is a great time to learn if your model is being specified in a way that you didn't expect or that 
# the model is not giving you the answers that you are expecting

gxeFit <- mxRun(gxe)
summary(gxeFit)

# Now, provided that the model looks reasonable, you can plug the model that you have built into the GWAS function

GWAS(model = gxe,                                                                 # what model object you would like to fit
	snpData = 'example.pgen',                                                     # that path to the snpData file.
	out="gxe.log")                                                                # the file that you would like to save the full results into



gxeResult <- read.delim("gxe.log")                            


# Then, we can read the results into R for a specific parameter. This function takes takes two arguments: the path to the data and the column in the results 
# file for the parameter that you want to examine.


succinctCond <- loadResults(path = "gxe.log", focus =  "snp_to_phe")
succinctInt  <- loadResults(path = "gxe.log", focus =  "snp_mod_to_phe")


head(succinctCond[order(succinctCond$Z, decreasing = T),])
head(succinctInt[order(succinctInt$Z, decreasing = T),])

gxeResult[hits$gxe, c("snp_to_phe", "snp_to_pheSE", "mod_to_phe",  "mod_to_pheSE", "snp_mod_to_phe", "snp_mod_to_pheSE")]



# Now it is possible to plot the interaction parameters
plot(succinctCond)
plot(succinctInt)





