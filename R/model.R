makeFitFunction <- function(fitfun, maxThr)
{
  if(maxThr==0 && fitfun == "WLS")     mxFitFunctionWLS(allContinuousMethod= "marginals")
  else if(maxThr>0 && fitfun == "WLS") mxFitFunctionWLS()
  else if(fitfun == "FIML")            mxFitFunctionML()
  else stop(paste("Unknown fitfun", omxQuotes(fitfun)))
}

calcMinVar <- function(minMAF) 2*minMAF*(1-minMAF)

# export TODO
makeComputePlan <- function(modelName, snpData, SNP, out)
{
  pieces <- strsplit(snpData, ".", fixed=TRUE)[[1]]
  if (length(pieces) < 2) {
    stop(paste("Please rename snpData",omxQuotes(snpData),
               "to the form file.ext where ext reflects the format of the data"))
  }
  snpFileExt <- pieces[length(pieces)]
  stem <- paste(pieces[-length(pieces)], collapse=".")

  if (snpFileExt == 'pgen' || snpFileExt == 'bed') method <- 'pgen'
  else if (snpFileExt == 'bgen') method <- 'bgen'
  else stop(paste("Unrecognized file extension", omxQuotes(snpFileExt),
                  "inferred from snpData", omxQuotes(snpData)))

  onesnp <- list(
    ST=mxComputeSetOriginalStarts(),
    LD=mxComputeLoadData(modelName, column='snp',
                         path=snpData, method=method))

  if (snpFileExt == "pgen") {
    # TODO doc column=1:3, sep='\t'
    onesnp <- c(
      onesnp,
      LC=mxComputeLoadContext(path=paste(stem, "pvar", sep = "."), column=1:3, sep='\t'))
  } else if (snpFileExt == "bed") {
    onesnp <- c(
      onesnp,
      LC=mxComputeLoadContext(path=paste(stem, "bim", sep = "."), column=1:3, sep='\t', header=FALSE))
  }

  onesnp <- c(
    onesnp,
    TC=mxComputeTryCatch(mxComputeSequence(list(
      GD=mxComputeGradientDescent(),
      SE=mxComputeStandardError()))),
    CK=mxComputeCheckpoint(path=paste(out, "log", sep = "."), standardErrors = TRUE))

  mxComputeLoop(onesnp, i=SNP)
}

# TODO: Replace runModel another layer of functions to build the models

#' Conduct a single factor genome-wide association study
#'
#' @template args-phenoData
#' @template args-snpData
#' @template args-snp
#' @template args-fitfun
#' @importFrom stats rbinom
#' @family GWAS
#' @export
oneFacGWAS <- function(phenoData, snpData, itemNames, covariates = NULL, SNP = NULL, fitfun = c("WLS","FIML"), minMAF = .01, out = "out", runModel=TRUE)
{
  fitfun <- match.arg(fitfun)
  minVar <- calcMinVar(minMAF)

  fac <- matrix(1, length(itemNames))
  thr <- matrix(1, length(itemNames))
  
  for(i in 1:length(itemNames)){
    fac[i] <-  is.factor(phenoData[,itemNames][,i])
    thr[i] <-  nlevels(phenoData[,itemNames][,i])-1
  }
  thr[thr< 0] <- 0
  maxThr <- max(thr)

  phenoData$snp <- rbinom(dim(phenoData)[1], 2, .5) # create placeholder
  latents   <- c("F")
  lambda    <- mxPath(from=latents, to=itemNames,values=1, free = T, labels = paste("lambda", itemNames, sep = "_")  )
  snpMu     <- mxPath(from = "one", to = "snp" , labels = "snpMean")
  snpBeta   <- mxPath(from = "snp", to = "F", labels = "snpReg", values = 0, free = T)
  snpres    <- mxPath(from = "snp", arrows=2, values=1, free = T, labels = paste("snp", "res", sep = "_"))
  resid     <- mxPath(from = c(itemNames), arrows=2, values=1, free = c(fac==0), labels = paste(c(itemNames), "res", sep = "_"))
  facRes    <- mxPath(from=latents, arrows=2,free=F, values=1.0, labels = "facRes")
  itemMean  <- mxPath(from = 'one', to = itemNames, free= c(fac==0), values = 0, labels = paste0(itemNames, "Mean"))

  if(maxThr>0) thresh    <- mxThreshold(itemNames[c(fac==1)], nThresh=c(thr[fac==1]), free = T , labels = paste(rep(itemNames[c(fac==1)], each = maxThr), "_Thr_", 1:maxThr, sep = ""), values=mxNormalQuantiles(1))
  
  
  dat       <- mxData(observed=phenoData, type="raw", minVariance=minVar, warnNPDacov=FALSE)

  fun <- makeFitFunction(fitfun, maxThr)

  modelName <- "OneFac"
  oneFacPre <- mxModel(model=modelName, type="RAM",
                       manifestVars = c("snp", itemNames),
                       latentVars = c(latents, covariates),
                       lambda, snpMu, snpBeta, snpres, resid, facRes,
                       itemMean, dat, fun  )

  if (length(covariates)) {
    covMean   <- mxPath(from = "one", to = covariates, free=FALSE, labels = paste0('data.',covariates)) 
    cov2item  <- mxPath(from = covariates, to = itemNames, connect = "all.pairs", labels = paste(rep(covariates, each = length(itemNames)), itemNames, sep = "_2_"))
    oneFacPre <- mxModel(oneFacPre, covMean, cov2item)
  }

  if(maxThr>0) oneFacPre <- mxModel(oneFacPre, name = "OneFac", thresh  )

  plan <- makeComputePlan(modelName, snpData, SNP, out)

  oneFac <- mxModel(oneFacPre, name = "OneFac", plan  )

  if (runModel) {
    oneFacFit <- mxRun(oneFac)
    summary(oneFacFit)
  } else { oneFac }
}

#' Conduct a single factor genome-wide association study with a focus on residuals
#' 
#' @template args-phenoData
#' @template args-snpData
#' @template args-snp
#' @template args-fitfun
#' @family GWAS
#' @export
oneFacResGWAS <- function(phenoData, snpData, itemNames , factor = F, res = itemNames, covariates = NULL, SNP = NULL, fitfun = c("WLS","FIML"), minMAF = .01, out = "out", runModel=TRUE) {
  fitfun <- match.arg(fitfun)
  minVar <- calcMinVar(minMAF)

  fac <- matrix(1, length(itemNames))
  thr <- matrix(1, length(itemNames))

  for(i in 1:length(itemNames)){
    fac[i] <-  is.factor(phenoData[,itemNames][,i])
    thr[i] <-  nlevels(phenoData[,itemNames][,i])-1
  }
  thr[thr< 0] <- 0
  maxThr <- max(thr)


  phenoData$snp <- rbinom(dim(phenoData)[1], 2, .5) # create placeholder
  latents   <- c("F")
  lambda    <- mxPath(from=latents, to=itemNames,values=1, free = T, labels = paste("lambda", itemNames, sep = "_")  )
  snpMu     <- mxPath(from = "one", to = "snp" , labels = "snpMean")
  snpFac  <- mxPath(from = "snp", to = "F", labels = "snp2fac", values = 0, free = factor)
  snpItemRes   <- mxPath(from = "snp", to = res, labels = paste("snp", res, sep = "2"), values = 0, free = T)
  snpres    <- mxPath(from = "snp", arrows=2, values=1, free = T, labels = paste("snp", "res", sep = "_"))
  resid     <- mxPath(from = c(itemNames), arrows=2, values=1, free = c(fac==0), labels = paste(c(itemNames), "res", sep = "_"))
  facRes    <- mxPath(from=latents, arrows=2,free=F, values=1.0, labels = "facRes")
  itemMean  <- mxPath(from = 'one', to = itemNames, free= c(fac==0), values = 0, labels = paste0(itemNames, "Mean"))

  
  if(maxThr>0) thresh    <- mxThreshold(itemNames[c(fac==1)], nThresh=c(thr[fac==1]), free = T , labels = paste(rep(itemNames[c(fac==1)], each = maxThr), "_Thr_", 1:maxThr, sep = ""), values=mxNormalQuantiles(1))
  dat       <- mxData(observed=phenoData, type="raw", minVariance=minVar, warnNPDacov=FALSE)

  fun <- makeFitFunction(fitfun, maxThr)

  modelName <- "OneFacRes"
  oneFacPre <- mxModel(model=modelName, type="RAM",
                       manifestVars = c("snp", itemNames),
                       latentVars = c(latents, covariates),
                       lambda, snpMu, snpFac, snpItemRes, snpres, resid, facRes,
                       itemMean, dat, fun  )

  if (length(covariates)) {
    covMean   <- mxPath(from = "one", to = covariates, free=FALSE, labels = paste0('data.',covariates)) 
    cov2item  <- mxPath(from = covariates, to = itemNames, connect = "all.pairs", labels = paste(rep(covariates, each = length(itemNames)), itemNames, sep = "_2_"))
    oneFacPre <- mxModel(oneFacPre, covMean, cov2item)
  }

  if(maxThr>0) oneFacPre <- mxModel(oneFacPre, name = "OneFacRes", thresh  )

  plan <- makeComputePlan(modelName, snpData, SNP, out)

  oneFac <- mxModel(oneFacPre, name = "OneFacRes", plan)
  if (runModel) {
    oneFacFit <- mxRun(oneFac)
    summary(oneFacFit)
  } else { oneFac }
}

#' Conduct a two factor genome-wide association study
#' 
#' @template args-phenoData
#' @template args-snpData
#' @template args-snp
#' @template args-fitfun
#' @family GWAS
#' @export
twoFacGWAS <- function(phenoData, snpData, F1itemNames, F2itemNames, covariates = NULL, SNP = NULL, fitfun = c("WLS","FIML"), minMAF = .01, out = "out", runModel=TRUE) {
  fitfun <- match.arg(fitfun)
  minVar <- calcMinVar(minMAF)

  itemNames <- c(F1itemNames,F2itemNames)
  fac <- matrix(1, length(itemNames))
  thr <- matrix(1, length(itemNames))
  
  for(i in 1:length(itemNames)){
    fac[i] <-  is.factor(phenoData[,itemNames][,i])
    thr[i] <-  nlevels(phenoData[,itemNames][,i])-1
  }
  thr[thr< 0] <- 0
  maxThr <- max(thr)

  phenoData$snp <- rbinom(dim(phenoData)[1], 2, .5) # create placeholder
  latents   <- c("F1", "F2")
  lambda1    <- mxPath(from="F1", to=F1itemNames,values=1, labels = paste("lambda", F1itemNames, sep = "_")  )
  lambda2    <- mxPath(from="F2", to=F2itemNames,values=1, labels = paste("lambda", F2itemNames, sep = "_")  )
  facCor    <- mxPath(from="F1", to= "F2", arrows=2,free=T, values=.3, labels = c("corrs"))

  snpMu     <- mxPath(from = "one", to = "snp" , labels = "snpMean")
  snpBeta   <- mxPath(from = "snp", to = latents, labels = paste0("snp", 2, latents), values = 0, free = T)
  snpres    <- mxPath(from = "snp", arrows=2, values=1, free = T, labels = paste("snp", "res", sep = "_"))

  resid     <- mxPath(from = itemNames, arrows=2, values=1, free = c(fac==0), labels = paste(c(itemNames), "res", sep = "_"))
  facRes    <- mxPath(from=latents, arrows=2,free=F, values=1.0, labels = "facRes")
  if (length(covariates)) {
    covMean   <- mxPath(from = "one", to = covariates, free=FALSE, labels = paste0('data.',covariates)) 
    cov2item  <- mxPath(from = covariates, to = c(itemNames), connect = "all.pairs", labels = paste(rep(covariates, each = length(itemNames)), itemNames, sep = "_2_"))
  }
  itemMean  <- mxPath(from = 'one', to = itemNames, free= c(fac==0), values = 0, labels = paste0(itemNames, "Mean"))

  if(maxThr>0) thresh    <- mxThreshold(c(F1itemNames, F2itemNames)[c(fac==1)], nThresh=c(thr[fac==1]), free = T , labels = paste(rep(itemNames[c(fac==1)], each = maxThr), "_Thr_", 1:maxThr, sep = ""), values=mxNormalQuantiles(1))
  
  
  dat       <- mxData(observed=phenoData, type="raw", minVariance=minVar, warnNPDacov=FALSE)

  fun <- makeFitFunction(fitfun, maxThr)

  modelName <- "TwoFac"
  twoFacPre <- mxModel(model=modelName, type="RAM",
                       manifestVars = c("snp", F1itemNames, F2itemNames),
                       latentVars = c(latents, covariates),
                       lambda1, lambda2, facCor, snpMu, snpBeta, snpres, 
                       resid, facRes, covMean, cov2item,
                       itemMean, dat, fun  )


  if(maxThr>0) twoFacPre <- mxModel(twoFacPre, name = "TwoFac", thresh  )

  plan <- makeComputePlan(modelName, snpData, SNP, out)

  twoFac <- mxModel(twoFacPre, name = "TwoFac", plan)

  if (runModel) {
    twoFacFit <- mxRun(twoFac)
    summary(twoFacFit)
  } else { twoFac }
}
