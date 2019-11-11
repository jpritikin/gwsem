makeFitFunction <- function(fitfun)
{
  if(fitfun == "WLS")        mxFitFunctionWLS(allContinuousMethod= "marginals")
  else if(fitfun == "ML")  mxFitFunctionML()
  else stop(paste("Unknown fitfun", omxQuotes(fitfun)))
}

calcMinVar <- function(minMAF) 2*minMAF*(1-minMAF)

#' Return a suitable compute plan for a genome-wide association study
#'
#' Instead of using OpenMx's default model processing sequence (i.e.,
#' \link[OpenMx]{omxDefaultComputePlan}), it is more efficient and
#' convienient to assemble a compute plan tailored for a genome-wide
#' association study.  This function returns a compute plan that loads
#' SNP data into model \code{modelName}, fits the model, outputs the
#' results to \code{out}, and repeats this procedure for all SNPs.
#'
#' @details
#'
#' You can request a specific list of SNPs using the \code{SNP}
#' argument. The numbers provided in \code{SNP} refer to offsets in
#' the \code{snpData} file. For example, \code{SNP=c(100,200)} will
#' process the 100th and 200th SNP. The first SNP in the
#' \code{snpData} file is at offset 1. When \code{SNP} is omitted then
#' all available SNPs are processed.
#' 
#' The suffix of \code{snpData} filename is interpreted to signal the
#' format of how the SNP data is stored on disk. Suffixes
#' \sQuote{pgen}, \sQuote{bed}, and \sQuote{bgen} are supported.
#' Per-SNP descriptions are found in different places depending on the
#' suffix. For \sQuote{bgen}, both the SNP data and description are
#' built into the same file. In the case of \sQuote{pgen}, an
#' associated file with suffix \sQuote{pvar} is expected to exist (see
#' the
#' \href{https://www.cog-genomics.org/plink/2.0/formats#pvar}{spec}
#' for details). In the case of \sQuote{bed}, an associated
#' \sQuote{bim} file is expected to exist (see the
#' \href{https://www.cog-genomics.org/plink2/formats#bim}{spec} for
#' details). The chromosome, base-pair coordinate, and variant ID are
#' added to each line of \code{out}.
#'
#' The code to implement method='pgen' is based on plink 2.0
#' alpha. plink's \sQuote{bed} file format is supported in addition
#' to \sQuote{pgen}. Data are coerced appropriately depending on the
#' type of the destination column. For a numeric column, data are
#' recorded as the values NA, 0, 1, or 2. An ordinal column must have
#' exactly 3 levels.
#'
#' For \code{method='bgen'}, the file \code{path+".bgi"} must also
#' exist. If not available, generate this index file with the
#' \href{https://bitbucket.org/gavinband/bgen/wiki/bgenix}{bgenix}
#' tool.
#'
#' For \sQuote{bgen} and \sQuote{pgen} formats, the numeric column can be
#' populated with a dosage (sum of probabilities multiplied by genotypes)
#' if these data are available.
#'
#' A compute plan does not do anything by itself. You'll need to combine
#' the compute plan with a model (such as returned by \link{buildOneFac})
#' to perform a GWAS.
#'
#' @template args-model
#' @template args-snpData
#' @template args-snp
#' @template args-out
#' @template args-dots-barrier
#' @template args-startfrom
#' @return
#' The given model with an appropriate compute plan.
#'
#' @export
#' @importFrom methods is
#' @seealso \link{GWAS}
#' @examples
#' m1 <- mxModel("test", mxFitFunctionWLS())
#' dir <- system.file("extdata", package = "gwsem")
#' m1 <- prepareComputePlan(m1, file.path(dir,"example.pgen"))
#' m1$compute
prepareComputePlan <- function(model, snpData, out="out.log", ...,
			       SNP=NULL, startFrom=1L)
{
  if (length(list(...)) > 0) stop("Rejected are any values passed in the '...' argument")
  modelName <- model$name
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
      LC=mxComputeLoadContext(path=paste(stem, "pvar", sep = "."), column=1:3, sep='\t',
			      col.names=c("CHR", "BP", "SNP")))
  } else if (snpFileExt == "bed") {
    onesnp <- c(
      onesnp,
      LC=mxComputeLoadContext(path=paste(stem, "bim", sep = "."),
                              column=c(1,2,4), sep='\t', header=FALSE,
                              col.names=c("CHR", "SNP", "BP")))
  }

  opt <- list(GD=mxComputeGradientDescent())
  if (is(model$fitfunction, "MxFitFunctionWLS")) {
	  opt <- c(opt, SE=mxComputeStandardError())
  } else {
	  opt <- c(opt,
		   ND=mxComputeNumericDeriv(),
		   SE=mxComputeStandardError(),
		   HQ=mxComputeHessianQuality())
  }

  onesnp <- c(
    onesnp,
    TC=mxComputeTryCatch(mxComputeSequence(opt)),
    CK=mxComputeCheckpoint(path=out, standardErrors = TRUE))

  mxModel(model, mxComputeLoop(onesnp, i=SNP, startFrom=startFrom))
}

#' Run a genome-wide association study (GWAS) using the provided model
#'
#' Adds a compute plan returned by \link{prepareComputePlan} to the
#' provided \code{model} and runs it. Once analyses are complete,
#' load your aggregated results with \link{loadResults}. 
#'
#' @template args-model
#' @template args-snpData
#' @template args-snp
#' @template args-out
#' @template args-dots-barrier
#' @template args-startfrom
#' @export
#' @return
#' The \link[OpenMx:MxModel-class]{MxModel} returned by \link[OpenMx]{mxRun}.
#' Data and estimates for the last SNP processed will be available for inspection.
#' @examples
#' dir <- system.file("extdata", package = "gwsem")
#' pheno <- data.frame(anxiety=rnorm(500))
#' m1 <- buildOneItem(pheno, 'anxiety')
#' GWAS(m1, file.path(dir,"example.pgen"),
#'      file.path(tempdir(),"out.log"))
GWAS <- function(model, snpData, out="out.log", ..., SNP=NULL, startFrom=1L)
{
	# verify model has a continuous 'snp' data column TODO
  if (length(list(...)) > 0) stop("Rejected are any values passed in the '...' argument")
  model <- prepareComputePlan(model, snpData, out=out,
			      SNP=SNP, startFrom=startFrom)
  model <- mxRun(model)
  message(paste("Done. See", omxQuotes(out), "for results"))
  invisible(model)
}

#' Set up thresholds for ordinal indicators
#'
#' Ordinal indicator thresholds are freely estimated with fixed means
#' and variance. This function adds thresholds to the given
#' \code{model}.  If no indicators are ordinal, the given \code{model}
#' is returned without changes.
#'
#' @details
#' 
#' Thresholds are added using \link[OpenMx]{mxThreshold}. Starting
#' values for thresholds use the defaults provided by this function
#' which assumes a mean of zero and variance of the square root of
#' two.  This variance is appropriate for \link{buildOneFac} where the
#' implied model variance of ordinal indicators is one plus the square
#' of the factor loading, and the loading's starting value is 1.0.
#'
#' @template args-model
#' @template detail-adv
#'
#' @return
#' The given \link[OpenMx:MxModel-class]{MxModel} with appropriate thresholds added.
#' @export
#' @examples
#' pheno <- data.frame(anxiety=cut(rnorm(500), c(-Inf, -.5, .5, Inf),
#'                     ordered_result = TRUE))
#' m1 <- buildOneItem(pheno, 'anxiety')
#' m1 <- setupThresholds(m1)
#' m1$Thresholds
setupThresholds <- function(model)
{
  phenoData <- model$data$observed
  itemNames <- setdiff(model$manifestVars, 'snp')

  thr <- sapply(phenoData[,itemNames,drop=FALSE], nlevels)-1
  fac <- thr >= 1
  thr[thr< 0] <- 0

  if (max(thr) == 0) return(model)

  thresh <- mxThreshold(itemNames[fac], nThresh=thr[fac], free = T ,
                        labels = paste0(rep(itemNames[fac], each = max(thr)), "_Thr_", 1:max(thr)),
			values=mxNormalQuantiles(thr[fac], sd=sqrt(2.0)))
  mxModel(model, thresh)
}

#' Set up model covariates
#'
#' In GWAS, including a number of the first principle components as
#' covariates helps reduce false positives caused by population
#' stratification. This function adds paths from covariates to
#' manifest indicators. Covariates are always treated as continuous
#' variables (not ordinal).
#'
#' @details
#' This is not the only way to adjust a model for
#' covariates. For example, in a single factor model (e.g., \link{buildOneFac}),
#' it would be more
#' appropriate to adjust the latent factor instead of the manifest
#' indicators. However, covariate adjustments to latent variables are only
#' possible with a maximum likelihood fit function
#' (\link[OpenMx]{mxFitFunctionML}).  For
#' \link[OpenMx]{mxFitFunctionWLS}, only manifest indicators can be
#' adjusted for covariates.
#' This function always adjusts manifest indicators regardless of the fit function.
#' 
#' @template args-model
#' @param covariates a character vector naming covariates available in the model data
#' @template detail-adv
#' @return
#' The given \link[OpenMx:MxModel-class]{MxModel} with paths
#' added from covariates to manifest indicators.
#' @export
#' @examples
#' m1 <- mxModel("test", type="RAM",
#'              latentVars = "sex", manifestVars = "anxiety")
#' m1 <- setupCovariates(m1, 'sex')
setupCovariates <- function(model, covariates)
{
  if (length(covariates)==0) return(model)

  itemNames <- setdiff(model$manifestVars, 'snp')
  covMean   <- mxPath(from = "one", to = covariates, free=FALSE, labels = paste0('data.',covariates)) 
  cov2item  <- mxPath(from = covariates, to = itemNames, connect = "all.pairs",
                      labels = paste(rep(covariates, each = length(itemNames)), itemNames, sep = "_"))
  mxModel(model, covMean, cov2item)
}

# export? TODO
setupData <- function(phenoData, gxe, customMinMAF, minMAF, fitfun)
{
  if (customMinMAF && fitfun != "WLS") warning("minMAF is ignored when fitfun != 'WLS'")
  minVar <- calcMinVar(minMAF)
  result <- list()
  aname <- c()
  if (length(gxe)) for (v1 in gxe) {
	  alg <- mxAlgebraFromString(paste0("data.snp * data.",v1),
				     name=paste0('snp_',v1,"Alg"),
				     dimnames=list(NULL,paste0('snp_',v1)))
	  result <- c(result, alg)
	  phenoData[[ paste0('snp_',v1) ]] <- 0.0  # placeholder
	  aname <- c(aname, paste0('snp_',v1))
  }
  c(mxData(observed=phenoData, type="raw", minVariance=minVar, warnNPDacov=FALSE,
	 algebra=aname), result)
}

#' @importFrom stats rbinom
addPlaceholderSNP <- function(phenoData) {
	# We use as.numeric because we currently only support dosages.
	phenoData$snp <- as.numeric(rbinom(dim(phenoData)[1], 2, .5))
	phenoData
}


#' Build a model suitable for a single item genome-wide association study
#'
#' @template detail-build
#'
#' @template args-phenoData
#' @param depVar the name of the single item to predict
#' @template args-covariates
#' @template args-dots-barrier
#' @template args-fitfun
#' @template args-minmaf
#' @template args-modeltype
#' @template args-gxe
#' @family model builder
#' @export
#' @return
#' A \link[OpenMx:MxModel-class]{MxModel}
#' @examples
#' pheno <- data.frame(anxiety=cut(rnorm(500), c(-Inf, -.5, .5, Inf),
#'                     ordered_result = TRUE))
#' m1 <- buildOneItem(pheno, 'anxiety')
buildOneItem <- function(phenoData, depVar, covariates=NULL, ..., fitfun = c("WLS","ML"), minMAF=0.01,
			 modelType=c('RAM','LISREL'), gxe=NULL)
{
  if (length(list(...)) > 0) stop("Rejected are any values passed in the '...' argument")
  fitfun <- match.arg(fitfun)
  modelType <- match.arg(modelType)

  if (length(depVar) != 1) {
    stop(paste("buildOneItem provided with", length(depVar), "dependent variables",
	       "instead of 1. Did you intend to use buildOneFac instead?"))
  }
  fac <- is.factor(phenoData[[depVar]])

  phenoData <- addPlaceholderSNP(phenoData)
  snpMu     <- mxPath(from = "one", to = "snp" , labels = "snpMean")
  snpBeta   <- mxPath(from = "snp", to = depVar, labels = "snpReg", values = 0, free = T)
  snpres    <- mxPath(from = "snp", arrows=2, values=1, free = T, labels = paste("snp", "res", sep = "_"))
  resid     <- mxPath(from = c(depVar), arrows=2, values=1, free = !fac, labels = paste(c(depVar), "res", sep = "_"))
  itemMean  <- mxPath(from = 'one', to = depVar, free= !fac, values = 0, labels = paste0(depVar, "Mean"))

  dat       <- setupData(phenoData, gxe, force(!missing(minMAF)), minMAF, fitfun)

  modelName <- "OneItem"
  oneFacPre <- mxModel(model=modelName, type=modelType,
                       manifestVars = c("snp", depVar),
                       latentVars = c(covariates),
                       snpMu, snpBeta, snpres, resid,
                       itemMean, dat, makeFitFunction(fitfun))

  oneFacPre <- setupThresholds(oneFacPre)
  setupCovariates(oneFacPre, covariates)
}

#' Build a model suitable for a single factor genome-wide association study
#'
#' @template detail-build
#'
#' @template args-phenoData
#' @param itemNames a vector of phenotypic item names that load on the latent factor
#' @template args-covariates
#' @template args-dots-barrier
#' @template args-fitfun
#' @template args-minmaf
#' @template args-modeltype
#' @template args-gxe
#' @family model builder
#' @export
#' @return
#' A \link[OpenMx:MxModel-class]{MxModel}
#' @examples
#' pheno <- list()
#' for (i in 1:5) pheno[[paste0('i',i)]] <- rnorm(500)
#' pheno <- as.data.frame(pheno)
#' buildOneFac(pheno, colnames(pheno))
buildOneFac <- function(phenoData, itemNames, covariates=NULL, ..., fitfun = c("WLS","ML"), minMAF=0.01,
			modelType=c('RAM','LISREL'), gxe=NULL)
{
  if (length(list(...)) > 0) stop("Rejected are any values passed in the '...' argument")
  fitfun <- match.arg(fitfun)
  modelType <- match.arg(modelType)

  fac <- sapply(phenoData[,itemNames,drop=FALSE], is.factor)

  phenoData <- addPlaceholderSNP(phenoData)
  latents   <- c("F")
  lambda    <- mxPath(from=latents, to=itemNames,values=1, free = T, labels = paste("lambda", itemNames, sep = "_")  )
  snpMu     <- mxPath(from = "one", to = "snp" , labels = "snpMean")
  snpBeta   <- mxPath(from = "snp", to = "F", labels = "snpReg", values = 0, free = T)
  snpres    <- mxPath(from = "snp", arrows=2, values=1, free = T, labels = paste("snp", "res", sep = "_"))
  resid     <- mxPath(from = c(itemNames), arrows=2, values=1, free = !fac, labels = paste(c(itemNames), "res", sep = "_"))
  facRes    <- mxPath(from=latents, arrows=2,free=F, values=1.0, labels = "facRes")
  itemMean  <- mxPath(from = 'one', to = itemNames, free= !fac, values = 0, labels = paste0(itemNames, "Mean"))

  dat       <- setupData(phenoData, gxe, force(!missing(minMAF)), minMAF, fitfun)

  modelName <- "OneFac"
  oneFacPre <- mxModel(model=modelName, type=modelType,
                       manifestVars = c("snp", itemNames),
                       latentVars = c(latents, covariates),
                       lambda, snpMu, snpBeta, snpres, resid, facRes,
                       itemMean, dat, makeFitFunction(fitfun))

  oneFacPre <- setupThresholds(oneFacPre)
  setupCovariates(oneFacPre, covariates)
}

#' Build a model suitable for a single factor residual genome-wide association study
#' 
#' @template detail-build
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
#' @template args-gxe
#' 
#' @family model builder
#' @export
#' @return
#' A \link[OpenMx:MxModel-class]{MxModel}
#' @examples
#' pheno <- list()
#' for (i in 1:5) pheno[[paste0('i',i)]] <- rnorm(500)
#' pheno <- as.data.frame(pheno)
#' buildOneFacRes(pheno, colnames(pheno))
buildOneFacRes <- function(phenoData, itemNames, factor = F, res = itemNames, covariates = NULL,
			   ..., fitfun = c("WLS","ML"), minMAF = .01, modelType=c('RAM','LISREL'), gxe=NULL)
{
  if (length(list(...)) > 0) stop("Rejected are any values passed in the '...' argument")
  fitfun <- match.arg(fitfun)
  modelType <- match.arg(modelType)
  
  fac <- sapply(phenoData[,itemNames,drop=FALSE], is.factor)

  phenoData <- addPlaceholderSNP(phenoData)
  latents   <- c("F")
  lambda    <- mxPath(from=latents, to=itemNames,values=1, free = T, labels = paste("lambda", itemNames, sep = "_")  )
  snpMu     <- mxPath(from = "one", to = "snp" , labels = "snpMean")
  snpFac  <- mxPath(from = "snp", to = "F", labels = "snp2fac", values = 0, free = factor)
  snpItemRes   <- mxPath(from = "snp", to = res, labels = paste("snp", res, sep = "2"), values = 0, free = T)
  snpres    <- mxPath(from = "snp", arrows=2, values=1, free = T, labels = paste("snp", "res", sep = "_"))
  resid     <- mxPath(from = c(itemNames), arrows=2, values=1, free = c(fac==0), labels = paste(c(itemNames), "res", sep = "_"))
  facRes    <- mxPath(from=latents, arrows=2,free=F, values=1.0, labels = "facRes")
  itemMean  <- mxPath(from = 'one', to = itemNames, free= c(fac==0), values = 0, labels = paste0(itemNames, "Mean"))

  dat       <- setupData(phenoData, gxe, force(!missing(minMAF)), minMAF, fitfun)

  modelName <- "OneFacRes"
  oneFacPre <- mxModel(model=modelName, type=modelType,
                       manifestVars = c("snp", itemNames),
                       latentVars = c(latents, covariates),
                       lambda, snpMu, snpFac, snpItemRes, snpres, resid, facRes,
                       itemMean, dat, makeFitFunction(fitfun))

  oneFacPre <- setupThresholds(oneFacPre)
  setupCovariates(oneFacPre, covariates)
}

#' Build a model suitable for a two factor genome-wide association study
#'
#' @template detail-build
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
#' @template args-gxe
#' @export
#' @family model builder
#' @return
#' A \link[OpenMx:MxModel-class]{MxModel}
#' @examples
#' pheno <- list()
#' for (i in 1:10) pheno[[paste0('i',i)]] <- rnorm(500)
#' pheno <- as.data.frame(pheno)
#' buildTwoFac(pheno, paste0('i',1:6), paste0('i',5:10))
buildTwoFac <- function(phenoData, F1itemNames, F2itemNames, covariates = NULL, ...,
			fitfun = c("WLS","ML"), minMAF = .01, modelType=c('RAM','LISREL'), gxe=NULL)
{
  if (length(list(...)) > 0) stop("Rejected are any values passed in the '...' argument")
  fitfun <- match.arg(fitfun)
  modelType <- match.arg(modelType)

  itemNames <- union(F1itemNames, F2itemNames)

  fac <- sapply(phenoData[,itemNames,drop=FALSE], is.factor)

  phenoData <- addPlaceholderSNP(phenoData)
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

  
  dat       <- setupData(phenoData, gxe, force(!missing(minMAF)), minMAF, fitfun)

  modelName <- "TwoFac"
  twoFacPre <- mxModel(model=modelName, type=modelType,
                       manifestVars = c("snp", itemNames),
                       latentVars = c(latents, covariates),
                       lambda1, lambda2, facCor, snpMu, snpBeta, snpres, 
                       resid, facRes,
                       itemMean, dat, makeFitFunction(fitfun))

  twoFacPre <- setupThresholds(twoFacPre)
  setupCovariates(twoFacPre, covariates)
}
