#
# This model that we will build on for GW-SEM
# 
require(gwsem)

# load the LGC phenotypic data
lgcData <- read.table("lgcData.txt", header = T)


## Build the model
manifests<-c("time1","time2","time3","time4","snp")
latents<-c("intercept","slope")
path <- list(mxPath(from="intercept",to=c("time1","time2","time3","time4"), free=c(FALSE,FALSE,FALSE,FALSE), value=c(1.0,1.0,1.0,1.0) , arrows=1, label=c("intercept__time1","intercept__time2","intercept__time3","intercept__time4") ),
             mxPath(from="slope",to=c("time2","time3","time4"), free=c(FALSE,FALSE,FALSE), value=c(1.0,2.0,3.0) , arrows=1, label=c("slope__time2","slope__time3","slope__time4") ),
             mxPath(from="one",to=c("intercept","slope"), free=c(TRUE,TRUE), value=c(1.0,1.0) , arrows=1, label=c("const__intercept","const__slope") ),
             mxPath(from="snp",to=c("intercept","slope"), free=c(TRUE,TRUE), value=c(0.5,0.1) , arrows=1, label=c("snp__intercept","snp__slope") ),
             mxPath(from="intercept",to=c("intercept","slope"), free=c(TRUE,TRUE), value=c(1.0,0.1) , arrows=2, label=c("sigma_i","COV_intercept_slope") ),
             mxPath(from="slope",to=c("slope"), free=c(TRUE), value=c(1.0) , arrows=2, label=c("sigma_s") ),
             mxPath(from="time1",to=c("time1"), free=c(TRUE), value=c(1.0) , arrows=2, label=c("res") ),
             mxPath(from="time2",to=c("time2"), free=c(TRUE), value=c(1.0) , arrows=2, label=c("res") ),
             mxPath(from="time3",to=c("time3"), free=c(TRUE), value=c(1.0) , arrows=2, label=c("res") ),
             mxPath(from="time4",to=c("time4"), free=c(TRUE), value=c(1.0) , arrows=2, label=c("res") ),
             mxPath(from="snp",to=c("snp"), free=c(TRUE), value=c(1.0) , arrows=2, label=c("VAR_snp") ),
             mxPath(from="one",to=c("time1","time2","time3","time4","snp"), free=F, value=0, arrows=1))
				 

lgcGWAS <- mxModel("lgc", type="RAM",
        manifestVars = manifests,
        latentVars = c(latents, paste0('pc', 1:5)),
        path,
        mxExpectationRAM(M="M"),
        mxFitFunctionWLS(allContinuousMethod="marginals"),
        mxData(observed=lgcData, type="raw", minVariance=0.1, warnNPDacov=FALSE))
		
lgcGWAS <- setupExogenousCovariates(lgcGWAS, paste0('pc', 1:5), paste0('time',1:4))


LGCtest <- mxRun(lgcGWAS)
summary(LGCtest)


GWAS(lgcGWAS, "example.pgen", "lgc.log")



LGC <- read.delim("lgc.log", header = T)


LGCresInt <- loadResults('lgc.log', 'snp__intercept')
LGCresSlo <- loadResults('lgc.log', 'snp__slope')

# In this example, we simulated snp1491 and snp1945 to be related to the latent intercept, and  snp1339 and snp1598 to be related to the latent slope, and snp1901 to be related to both the intercept and the slope.

head(LGCresInt[order(LGCresInt$Z, decreasing = T),])
head(LGCresSlo[order(abs(LGCresSlo$Z), decreasing = T),])