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

#' Run a genome-wide association study (GWAS) using the provided model
#'
#' @param model the MxModel object to be fit to each SNP
#' @template args-snpData
#' @template args-snp
#' @param out the filename where write shall be written
#' @export
GWAS <- function(model, snpData, SNP=NULL, out="out")
{
  model <- mxModel(model, makeComputePlan(model$name, snpData, SNP, out))
  model <- mxRun(model)
  message(paste("Done. See", omxQuotes(out), "for results"))
  invisible(model)
}

# TODO doc scale of thresholds
setupThresholds <- function(model, fac)
{
  phenoData <- model$data$observed
  itemNames <- setdiff(model$manifestVars, 'snp')

  thr <- sapply(phenoData[,itemNames], nlevels)-1
  thr[thr< 0] <- 0

  if (max(thr) == 0) return(model)

  thresh <- mxThreshold(itemNames[fac], nThresh=thr[fac], free = T ,
                        labels = paste0(rep(itemNames[fac], each = max(thr)), "_Thr_", 1:max(thr)))
  mxModel(model, thresh)
}

setupCovariates <- function(model, covariates)
{
  if (length(covariates)==0) return(model)

  itemNames <- setdiff(model$manifestVars, 'snp')
  covMean   <- mxPath(from = "one", to = covariates, free=FALSE, labels = paste0('data.',covariates)) 
  cov2item  <- mxPath(from = covariates, to = itemNames, connect = "all.pairs",
                      labels = paste(rep(covariates, each = length(itemNames)), itemNames, sep = "_"))
  mxModel(model, covMean, cov2item)
}

#' Build a model suitable for a single factor genome-wide association study
#'
#' @template args-phenoData
#' @param itemNames a vector of phenotypic item names that load on the latent factor
#' @template args-covariates
#' @template args-dots-barrier
#' @template args-fitfun
#' @template args-minmaf
#' @template args-modeltype
#' @family build
#' @importFrom stats rbinom
#' @export
buildOneFac <- function(phenoData, itemNames, covariates=NULL, ..., fitfun = c("WLS","ML"), minMAF=0.01,
			modelType=c('RAM','LISREL'))
{
  if (length(list(...)) > 0) stop("Rejected are any values passed in the '...' argument")
  fitfun <- match.arg(fitfun)
  minVar <- calcMinVar(minMAF)
  modelType <- match.arg(modelType)

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
  oneFacPre <- mxModel(model=modelName, type=modelType,
                       manifestVars = c("snp", itemNames),
                       latentVars = c(latents, covariates),
                       lambda, snpMu, snpBeta, snpres, resid, facRes,
                       itemMean, dat, makeFitFunction(fitfun))

  oneFacPre <- setupThresholds(oneFacPre, fac)
  setupCovariates(oneFacPre, covariates)
}

#' Build a model suitable for a single factor residual genome-wide association study
#' 
#' @param itemNames a vector of phenotypic item names that load on the latent factor
#' @param factor whether to estimate a regression from the SNP to the latent factor (default FALSE)
#' @param res character vector. Which indicators to estimate a regression to
#' @template args-phenoData
#' @template args-covariates
#' @template args-fitfun
#' @template args-minmaf
#' @template args-dots-barrier
#' @template args-modeltype
#' 
#' @family build
#' @export
buildOneFacRes <- function(phenoData, itemNames, factor = F, res = itemNames, covariates = NULL,
			   ..., fitfun = c("WLS","ML"), minMAF = .01, modelType=c('RAM','LISREL'))
{
  if (length(list(...)) > 0) stop("Rejected are any values passed in the '...' argument")
  fitfun <- match.arg(fitfun)
  if (!missing(minMAF) && fitfun != "WLS") warning("minMAF is ignored when fitfun != 'WLS'")
  modelType <- match.arg(modelType)
  
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
  oneFacPre <- mxModel(model=modelName, type=modelType,
                       manifestVars = c("snp", itemNames),
                       latentVars = c(latents, covariates),
                       lambda, snpMu, snpFac, snpItemRes, snpres, resid, facRes,
                       itemMean, dat, makeFitFunction(fitfun))

  oneFacPre <- setupThresholds(oneFacPre, fac)
  setupCovariates(oneFacPre, covariates)
}

#' Build a model suitable for a two factor genome-wide association study
#'
#' @param F1itemNames a vector of item names that load on the first latent factor
#' @param F2itemNames a vector of item names that load on the second latent factor
#' 
#' @template args-phenoData
#' @template args-covariates
#' @template args-fitfun
#' @template args-minmaf
#' @template args-dots-barrier
#' @template args-modeltype
#' @export
#' @family build
buildTwoFac <- function(phenoData, F1itemNames, F2itemNames, covariates = NULL, ...,
			fitfun = c("WLS","ML"), minMAF = .01, modelType=c('RAM','LISREL'))
{
  if (length(list(...)) > 0) stop("Rejected are any values passed in the '...' argument")
  fitfun <- match.arg(fitfun)
  if (!missing(minMAF) && fitfun != "WLS") warning("minMAF is ignored when fitfun != 'WLS'")
  modelType <- match.arg(modelType)

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
  twoFacPre <- mxModel(model=modelName, type=modelType,
                       manifestVars = c("snp", itemNames),
                       latentVars = c(latents, covariates),
                       lambda1, lambda2, facCor, snpMu, snpBeta, snpres, 
                       resid, facRes,
                       itemMean, dat, makeFitFunction(fitfun))

  twoFacPre <- setupThresholds(twoFacPre, fac)
  setupCovariates(twoFacPre, covariates)
}
