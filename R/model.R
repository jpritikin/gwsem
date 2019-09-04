makeFitFunction <- function(fitfun)
{
  if(fitfun == "WLS")        mxFitFunctionWLS(allContinuousMethod= "marginals")
  else if(fitfun == "ML")  mxFitFunctionML()
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
    # TODO confirm this is the correct order
    # disagrees with https://www.cog-genomics.org/plink2/formats#bim
    onesnp <- c(
      onesnp,
      LC=mxComputeLoadContext(path=paste(stem, "bim", sep = "."),
                              column=c(1,2,4), sep='\t', header=FALSE,
                              col.names=c("CHROM", "SNP", "POS")))
  }

  onesnp <- c(
    onesnp,
    TC=mxComputeTryCatch(mxComputeSequence(list(
      GD=mxComputeGradientDescent(),
      SE=mxComputeStandardError()))),
    CK=mxComputeCheckpoint(path=paste(out, "log", sep = "."), standardErrors = TRUE))

  mxComputeLoop(onesnp, i=SNP)
}

setupThresholds <- function(model, fac)
{
  phenoData <- model$data$observed
  itemNames <- setdiff(model$manifestVars, 'snp')

  thr <- sapply(phenoData[,itemNames], nlevels)-1
  thr[thr< 0] <- 0

  if (max(thr) == 0) return(model)

  thresh <- mxThreshold(itemNames[fac], nThresh=thr[fac], free = T ,
                        labels = paste0(rep(itemNames[fac], each = max(thr)), "_Thr_", 1:max(thr)))
  mxModel(model, name = "OneFac", thresh)
}

setupCovariates <- function(model, covariates)
{
  if (length(covariates)==0) return(model)

  itemNames <- setdiff(model$manifestVars, 'snp')
  covMean   <- mxPath(from = "one", to = covariates, free=FALSE, labels = paste0('data.',covariates)) 
  cov2item  <- mxPath(from = covariates, to = itemNames, connect = "all.pairs",
                      labels = paste(rep(covariates, each = length(itemNames)), itemNames, sep = "_2_"))
  mxModel(model, covMean, cov2item)
}

#' Build a model suitable for a single factor genome-wide association study
#'
#' @template args-phenoData
#' @template args-snpData
#' @template args-snp
#' @template args-fitfun
#' @importFrom stats rbinom
#' @export
buildOneFac <- function(phenoData, snpData, itemNames, covariates=NULL,
                        SNP=NULL, fitfun = c("WLS","ML"), minMAF=0.01, out="out")
{
  fitfun <- match.arg(fitfun)
  minVar <- calcMinVar(minMAF)

  fac <- sapply(phenoData[,itemNames], is.factor)

  phenoData$snp <- rbinom(dim(phenoData)[1], 2, .5) # create placeholder
  latents   <- c("F")
  lambda    <- mxPath(from=latents, to=itemNames,values=1, free = T, labels = paste("lambda", itemNames, sep = "_")  )
  snpMu     <- mxPath(from = "one", to = "snp" , labels = "snpMean")
  snpBeta   <- mxPath(from = "snp", to = "F", labels = "snpReg", values = 0, free = T)
  snpres    <- mxPath(from = "snp", arrows=2, values=1, free = T, labels = paste("snp", "res", sep = "_"))
  resid     <- mxPath(from = c(itemNames), arrows=2, values=1, free = !fac, labels = paste(c(itemNames), "res", sep = "_"))
  facRes    <- mxPath(from=latents, arrows=2,free=F, values=1.0, labels = "facRes")
  itemMean  <- mxPath(from = 'one', to = itemNames, free= !fac, values = 0, labels = paste0(itemNames, "Mean"))

  dat       <- mxData(observed=phenoData, type="raw", minVariance=minVar, warnNPDacov=FALSE)

  modelName <- "OneFac"
  oneFacPre <- mxModel(model=modelName, type="RAM",
                       manifestVars = c("snp", itemNames),
                       latentVars = c(latents, covariates),
                       lambda, snpMu, snpBeta, snpres, resid, facRes,
                       itemMean, dat, makeFitFunction(fitfun))

  oneFacPre <- setupThresholds(oneFacPre, fac)
  oneFacPre <- setupCovariates(oneFacPre, covariates)

  plan <- makeComputePlan(modelName, snpData, SNP, out)

  mxModel(oneFacPre, name = "OneFac", plan)
}

#' Conduct a single factor genome-wide association study
#'
#' @template args-phenoData
#' @template args-snpData
#' @template args-snp
#' @template args-fitfun
#' @family GWAS
#' @export
oneFacGWAS <- function(phenoData, snpData, itemNames, covariates = NULL, SNP = NULL, fitfun = c("WLS","ML"), minMAF = .01, out = "out")
{
  fitfun <- match.arg(fitfun)
  if (!missing(minMAF) && fitfun != "WLS") warning("minMAF is ignored when fitfun != 'WLS'")
  oneFac <- buildOneFac(phenoData, snpData, itemNames, covariates, SNP, fitfun, minMAF, out)
  oneFacFit <- mxRun(oneFac)
  summary(oneFacFit)
}

#' Build a model suitable for a single factor residual genome-wide association study
#' 
#' @export
buildOneFacRes <- function(phenoData, snpData, itemNames , factor = F, res = itemNames, covariates = NULL, SNP = NULL, fitfun = c("WLS","ML"), minMAF = .01, out = "out")
{
  fitfun <- match.arg(fitfun)
  if (!missing(minMAF) && fitfun != "WLS") warning("minMAF is ignored when fitfun != 'WLS'")
  
  fac <- sapply(phenoData[,itemNames], is.factor)

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

  minVar <- calcMinVar(minMAF)
  dat       <- mxData(observed=phenoData, type="raw", minVariance=minVar, warnNPDacov=FALSE)

  modelName <- "OneFacRes"
  oneFacPre <- mxModel(model=modelName, type="RAM",
                       manifestVars = c("snp", itemNames),
                       latentVars = c(latents, covariates),
                       lambda, snpMu, snpFac, snpItemRes, snpres, resid, facRes,
                       itemMean, dat, makeFitFunction(fitfun))

  oneFacPre <- setupThresholds(oneFacPre, fac)
  oneFacPre <- setupCovariates(oneFacPre, covariates)

  plan <- makeComputePlan(modelName, snpData, SNP, out)

  mxModel(oneFacPre, name = "OneFacRes", plan)
}

#' Conduct a single factor genome-wide association study with a focus on residuals
#' 
#' @template args-phenoData
#' @template args-snpData
#' @template args-snp
#' @template args-fitfun
#' @family GWAS
#' @export
oneFacResGWAS <- function(phenoData, snpData, itemNames , factor = F, res = itemNames, covariates = NULL, SNP = NULL, fitfun = c("WLS","ML"), minMAF = .01, out = "out") {
  fitfun <- match.arg(fitfun)
  if (!missing(minMAF) && fitfun != "WLS") warning("minMAF is ignored when fitfun != 'WLS'")
  oneFac <- buildOneFacRes(phenoData, snpData, itemNames, factor, res, covariates, SNP, fitfun, minMAF, out)
  oneFacFit <- mxRun(oneFac)
  summary(oneFacFit)
}

#' Build a model suitable for a two factor genome-wide association study
#'
#' @export
buildTwoFac <- function(phenoData, snpData, F1itemNames, F2itemNames, covariates = NULL, SNP = NULL, fitfun = c("WLS","ML"), minMAF = .01, out = "out")
{
  fitfun <- match.arg(fitfun)
  if (!missing(minMAF) && fitfun != "WLS") warning("minMAF is ignored when fitfun != 'WLS'")

  minVar <- calcMinVar(minMAF)

  itemNames <- union(F1itemNames, F2itemNames)

  fac <- sapply(phenoData[,itemNames], is.factor)

  phenoData$snp <- rbinom(dim(phenoData)[1], 2, .5) # create placeholder
  latents   <- c("F1", "F2")
  lambda1    <- mxPath(from="F1", to=F1itemNames,values=1, labels = paste("lambda", F1itemNames, sep = "_")  )
  lambda2    <- mxPath(from="F2", to=F2itemNames,values=1, labels = paste("lambda", F2itemNames, sep = "_")  )
  facCor    <- mxPath(from="F1", to= "F2", arrows=2,free=T, values=.3)

  snpMu     <- mxPath(from = "one", to = "snp" , labels = "snpMean")
  snpBeta   <- mxPath(from = "snp", to = latents, labels = paste0("snp", 2, latents), values = 0, free = T)
  snpres    <- mxPath(from = "snp", arrows=2, values=1, free = T, labels = paste("snp", "res", sep = "_"))

  resid     <- mxPath(from = itemNames, arrows=2, values=1, free = c(fac==0), labels = paste(c(itemNames), "res", sep = "_"))
  facRes    <- mxPath(from=latents, arrows=2,free=F, values=1.0, labels = "facRes")
  itemMean  <- mxPath(from = 'one', to = itemNames, free= c(fac==0), values = 0, labels = paste0(itemNames, "Mean"))

  
  dat       <- mxData(observed=phenoData, type="raw", minVariance=minVar, warnNPDacov=FALSE)

  modelName <- "TwoFac"
  twoFacPre <- mxModel(model=modelName, type="RAM",
                       manifestVars = c("snp", itemNames),
                       latentVars = c(latents, covariates),
                       lambda1, lambda2, facCor, snpMu, snpBeta, snpres, 
                       resid, facRes,
                       itemMean, dat, makeFitFunction(fitfun))

  twoFacPre <- setupThresholds(twoFacPre, fac)
  twoFacPre <- setupCovariates(twoFacPre, covariates)

  plan <- makeComputePlan(modelName, snpData, SNP, out)

  mxModel(twoFacPre, name = "TwoFac", plan)
}

#' Conduct a two factor genome-wide association study
#' 
#' @template args-phenoData
#' @template args-snpData
#' @template args-snp
#' @template args-fitfun
#' @family GWAS
#' @export
twoFacGWAS <- function(phenoData, snpData, F1itemNames, F2itemNames, covariates = NULL, SNP = NULL, fitfun = c("WLS","ML"), minMAF = .01, out = "out") {
  fitfun <- match.arg(fitfun)
  if (!missing(minMAF) && fitfun != "WLS") warning("minMAF is ignored when fitfun != 'WLS'")
  twoFac <- buildTwoFac(phenoData, snpData, F1itemNames, F2itemNames, covariates, SNP, fitfun, minMAF, out)
  twoFacFit <- mxRun(twoFac)
  summary(twoFacFit)
}
