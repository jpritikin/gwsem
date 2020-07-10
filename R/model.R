defaultExogenous <- TRUE

makeFitFunction <- function(fitfun)
{
  if(fitfun == "WLS")        mxFitFunctionWLS(allContinuousMethod= "marginals")
  else if(fitfun == "ML")  mxFitFunctionML()
  else stop(paste("Unknown fitfun", omxQuotes(fitfun)))
}

calcMinVar <- function(minMAF) 2*minMAF*(1-minMAF)

forModels <- function(topModel, modelName, fn) {
  ret <- c()
  if (topModel$name %in% modelName) {
    ret <- c(fn(topModel))
    names(ret) <- topModel$name
  }
  if (length(modelName) > 1) {
    ret <- c(ret, sapply(setdiff(modelName, topModel$name), function(x) fn(topModel[[x]])))
  }
  ret
}

#' Return a suitable compute plan for a genome-wide association study
#'
#' \lifecycle{maturing}
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
#' pheno <- data.frame(anxiety=cut(rnorm(500), c(-Inf, -.5, .5, Inf),
#' ordered_result = TRUE))
#' m1 <- buildItem(pheno, 'anxiety')
#' dir <- system.file("extdata", package = "gwsem")
#' m1 <- prepareComputePlan(m1, file.path(dir,"example.pgen"))
#' m1$compute
prepareComputePlan <- function(model, snpData, out="out.log", ...,
			       SNP=NULL, startFrom=1L)
{
  if (length(list(...)) > 0) stop("Rejected are any values passed in the '...' argument")

  loc <- strsplit(snpData, ".", fixed=TRUE)
  locLen <- sapply(loc, length)
  if (any(locLen < 2)) {
    stop(paste("Please rename snpData", omxQuotes(snpData[locLen < 2]),
               "to the form file.ext where ext reflects the format of the data"))
  }
  snpFileExt <- mapply(function(pieces, l) pieces[l],
                       loc, locLen)
  stem <- mapply(function(pieces, l) paste(pieces[-l], collapse="."),
                       loc, locLen)
  method <- sapply(snpFileExt, function(ext1) {
    if (snpFileExt == 'pgen' || snpFileExt == 'bed') 'pgen'
    else if (snpFileExt == 'bgen') 'bgen'
    else NA
  })
  if (any(is.na(method))) {
    stop(paste("Unrecognized file extension", omxQuotes(snpFileExt[is.na(method)]),
               "inferred from snpData", omxQuotes(snpData[is.na(method)])))
  }

  if (length(snpData) > 1 && length(names(snpData)) == 0)
    stop(paste("Must provide model names for snpData. For example,\n",
               "c(",omxQuotes(model$name),"=",snpData[1],")"))
  modelName <- model$name
  if (length(names(snpData))) modelName <- names(snpData)

  onesnp <- mapply(function (mn1, snpData1, ext1, method1, stem1) {
    p1 <- list(LD=mxComputeLoadData(mn1, column='snp',
                              path=snpData1, method=method1))

    if (ext1 == "pgen") {
      p1 <- c(
        p1,
        LC=mxComputeLoadContext(path=paste(stem1, "pvar", sep = "."), column=1:5, sep='\t',
                                col.names=c("CHR", "BP", "SNP", "A1", "A2")))
    } else if (ext1 == "bed") {
      p1 <- c(
        p1,
        LC=mxComputeLoadContext(path=paste(stem1, "bim", sep = "."),
                                column=c(1,2,4:6), sep='\t', header=FALSE,
                                col.names=c("CHR", "SNP", "BP", "A2", "A1")))
    } else {
      p1 # built-in to BGEN already
    }
  }, modelName, snpData, snpFileExt, method, stem)

  opt <- list(GD=mxComputeGradientDescent())
  if (is(model$fitfunction, "MxFitFunctionWLS")) {
	  opt <- c(opt, SE=mxComputeStandardError())
  } else {
	  opt <- c(opt,
		   ND=mxComputeNumericDeriv(),
		   SE=mxComputeStandardError(),
		   HQ=mxComputeHessianQuality())
  }

  wantVcov <- any(forModels(model, modelName, function(m) {
    d1 <- m$data
    if (is.null(d1)) stop(paste("Model",omxQuotes(m$name),"contains no MxData"))
    obs <- d1$observed
    if (!('snp' %in% colnames(obs))) {
      stop(paste("No snp placeholder column in observed data of model", omxQuotes(m$name)))
    }
    length(d1$algebra) > 0
  }))

  onesnp <- c(
    ST=mxComputeSetOriginalStarts(),
    onesnp,
    TC=mxComputeTryCatch(mxComputeSequence(opt)),
    CK=mxComputeCheckpoint(path=out, standardErrors = TRUE, vcov = wantVcov))

  mxModel(model, mxComputeLoop(onesnp, i=SNP, startFrom=startFrom))
}

#' Run a genome-wide association study (GWAS) using the provided model
#'
#' \lifecycle{maturing}
#' The GWAS function is used to run a genome-wide association study based on the specified model. This function is design to take the output from \link{buildOneFac}, \link{buildOneFacRes}, and \link{buildTwoFac} as input, but can also take a similar user specified model. Users should be confident that the models they are running are statistically identified. It is advisable that the users empirically gauge time requirements by running a limited number of SNPs (e.g. 10) to ensure that all SNPs can be fit in a reasonable amount of time.
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
#' The results for each SNP are recorded in the specified log file (\code{out}).
#' In addition, data and estimates for the last SNP run are returned
#' as an \link[OpenMx:MxModel-class]{MxModel} object
#' (similar to the return value of \link[OpenMx]{mxRun}).
#' In this way, the last SNP processed is available for close inspection.
#' @examples
#' dir <- system.file("extdata", package = "gwsem")
#' pheno <- data.frame(anxiety=rnorm(500))
#' m1 <- buildItem(pheno, 'anxiety')
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
#' \lifecycle{experimental}
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
#' m1 <- buildItem(pheno, 'anxiety')
#' m1 <- setupThresholds(m1)
#' m1$Thresholds
setupThresholds <- function(model)
{
  phenoData <- model$data$observed
  manifestNames <- setdiff(model$manifestVars, 'snp')

  thr <- sapply(phenoData[,manifestNames,drop=FALSE], nlevels)-1
  fac <- thr >= 1
  thr[thr< 0] <- 0

  if (max(thr) == 0) return(model)

  thresh <- list()
  for (m1 in manifestNames[fac]) {
	  thresh <- c(thresh,
		      mxThreshold(m1, nThresh=thr[m1], free = T ,
				  labels = paste0(m1, "_Thr_", 1:max(thr[m1])),
				  values=mxNormalQuantiles(thr[m1], sd=sqrt(2.0))))
  }
  mxModel(model, thresh)
}

#' Set up exogenous model covariates
#'
#' \lifecycle{experimental}
#' In GWAS, including a number of the first principle components as
#' covariates helps reduce false positives caused by population
#' stratification. This function adds paths from covariates to
#' manifest indicators (\code{itemNames}). Covariates are always treated as continuous
#' variables (not ordinal).
#'
#' @details
#' This is not the only way to adjust a model for
#' covariates. For example, in a single factor model (e.g., \link{buildOneFac}),
#' it would be more
#' appropriate to adjust the latent factor instead of the manifest
#' indicators.
#' This is how endogenous covariates work.
#' However, exogenous covariate adjustments to latent variables are only
#' possible with a maximum likelihood fit function
#' (\link[OpenMx]{mxFitFunctionML}).  For
#' \link[OpenMx]{mxFitFunctionWLS}, only manifest indicators can be
#' adjusted for exogenous covariates.
#' This function always adjusts manifest indicators regardless of the fit function.
#'
#' @template args-model
#' @param covariates a character vector naming covariates available in the model data
#' @param itemNames a character vector of item names
#' @template detail-adv
#' @return
#' The given \link[OpenMx:MxModel-class]{MxModel} with paths
#' added from covariates to manifest indicators.
#' @export
#' @examples
#' m1 <- mxModel("test", type="RAM",
#'              latentVars = "sex", manifestVars = "anxiety",
#'              mxData(data.frame(sex=rbinom(10,1,.5)), 'raw'))
#' m1 <- setupExogenousCovariates(m1, 'sex', 'anxiety')
setupExogenousCovariates <- function(model, covariates, itemNames)
{
  if (length(covariates)==0) return(model)

  covMean   <- mxPath(from = "one", to = covariates, free=FALSE, labels = paste0('data.',covariates))
  cov2item  <- mxPath(from = covariates, to = itemNames, connect = "all.pairs",
                      labels = paste(rep(covariates, each = length(itemNames)), itemNames,
                                     sep = "_to_"))
  model <- mxModel(model, covMean, cov2item)
  emean <- mxGetExpected(model, "means")
  exoFree <- matrix(FALSE, length(emean), length(covariates),
		    dimnames=list(colnames(emean), covariates))
  exoFree[itemNames, covariates] <- TRUE
  model@data$exoFree <- exoFree
  model
}

# export? TODO
setupData <- function(phenoData, gxe, customMinMAF, minMAF, fitfun)
{
  phenoData <- as.data.frame(phenoData)
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
	  aname <- c(aname, paste0('snp_',v1,"Alg"))
  }
  c(mxData(observed=phenoData, type="raw", minVariance=minVar, warnNPDacov=FALSE,
	   algebra=aname, naAction='omit'), result)
}

#' @importFrom stats rbinom
addPlaceholderSNP <- function(phenoData) {
	if (!is.null(phenoData[[ 'snp' ]])) {
		warning("Data already contains placeholder data for the 'snp' column. This is okay for testing, but generally not recommended")
	} else {
		# We use as.numeric because we currently only support dosages.
		phenoData$snp <- as.numeric(rbinom(dim(phenoData)[1], 2, .5))
	}
	phenoData
}

endogenousSNPpath <- function(pred, depVar)
{
	paths <- list(mxPath(from = "one", to = pred, labels = paste0(pred, 'Mean')),
                mxPath(from = pred, to = depVar, values = 0,
                       labels=paste(pred, depVar, sep = "_to_")),
                mxPath(from = pred, arrows=2, values=1, lbound=0.001,
                       labels = paste(pred, "res", sep = "_")))
}

endogenousCovariatePaths <- function(phenoData, covariates, depVar)
{
	paths <- list()

	if (length(covariates)) for (cx in 1:length(covariates)) {
		c1 <- covariates[cx]
		if (is.factor(phenoData[[c1]])) {
			nth <- length(levels(phenoData[[c1]]))-1
			paths <- c(paths, list(
				mxPath(from=c1, arrows=2, values=1, free=FALSE)))
			# see setupThresholds
		} else {
			paths <- c(paths, list(
				mxPath(from='one',c1,labels=paste0(c1,"_mean")),
				mxPath(from=c1, arrows=2, values=1, lbound=0.001, labels=paste0(c1,"_var"))))
		}
		paths <- c(paths, list(
			mxPath(from=c1, to=depVar, labels=paste0(c1,'_to_',depVar))))
	}
	paths
}

#' @importFrom OpenMx omxQuotes
postprocessModel <- function(model, indicators, exogenous)
{
	model <- setupThresholds(model)

  isInt <- sapply(model$data$observed, function(x) typeof(x) == 'integer')

	if (length(exogenous) == 0 &&
		    is.na(model$expectation$thresholds) &&
		    is(model$fitfunction, 'MxFitFunctionWLS')) {
    if (any(isInt)) {
      warning(paste("Cumulants method prevented by integer observed columns:",
                    omxQuotes(names(isInt)[isInt]),
                    "; Remove these columns from your observed data for",
                    "improved estimation accuracy"))
    } else {
      # cumulants is substantially more precise than marginals
      model$fitfunction$continuousType <- 'cumulants'

      # discard model means
      if (is(model$expectation, "MxExpectationRAM")) {
        model$expectation$M <- as.character(NA)
        model <- mxModel(model, 'M', remove = TRUE)
      } else {
        stop("Only RAM models are supported")
      }
    }
	} else {
		model <- setupExogenousCovariates(model, exogenous, indicators)
	}
	model
}

#' Build a model suitable for a single item genome-wide association study
#'
#' \lifecycle{maturing}
#' @template detail-build
#'
#' @section WLS Technical Note:
#' When the \code{depVar} item is/are continuous,
#' covariates are endogenous (the default),
#' and the fit function is \code{WLS} then the
#' \code{cumulants} method is used to create observed summary
#' statistics (see \link[OpenMx]{mxFitFunctionWLS}). In other cases,
#' the \code{marginals} method is used. The \code{cumulants} method is
#' more accurate than \code{marginals}. The difference in accuracy
#' becomes vivid when comparing estimates against the \code{ML} fit
#' function.
#'
#' @template args-phenoData
#' @param depVar the name of items to predict
#' @template args-covariates
#' @template args-exogenous
#' @template args-dots-barrier
#' @template args-fitfun
#' @template args-minmaf
#' @template args-gxe
#' @template args-pred
#' @family model builder
#' @export
#' @return
#' A \link[OpenMx:MxModel-class]{MxModel}
#' @aliases buildOneItem
#' @examples
#' pheno <- data.frame(anxiety=cut(rnorm(500), c(-Inf, -.5, .5, Inf),
#'                     ordered_result = TRUE))
#' m1 <- buildItem(pheno, 'anxiety')
buildItem <- function(phenoData, depVar, covariates=NULL, ..., fitfun = c("WLS","ML"), minMAF=0.01,
			 gxe=NULL, exogenous=NA, pred = 'snp')
{
  if (length(list(...)) > 0) stop("Rejected are any values passed in the '...' argument")
  fitfun <- match.arg(fitfun)

  phenoData <- addPlaceholderSNP(phenoData)
  # Remove extraneous data columns that could prevent WLS cumulants
  phenoData <- phenoData[,intersect(colnames(phenoData),
				    c('snp', depVar, covariates, gxe))]
  fac <- sapply(phenoData[,depVar,drop=FALSE], is.factor)

  if (is.na(exogenous)) {
	  if (fitfun == 'WLS' && !any(fac)) exogenous <- FALSE
	  else exogenous <- TRUE
  }

  manifest <- depVar
  latents <- c()
  endoCovariates <- c()
  if (length(gxe)) pred <- c(pred, paste0('snp_', gxe))
  if (!exogenous) {
	  manifest <- c(manifest, pred, covariates)
	  endoCovariates <- covariates
  } else {
	  latents <- c(pred, covariates)
  }

  paths <- endogenousCovariatePaths(phenoData, endoCovariates, depVar)
  if (!exogenous) {
    paths <- c(paths, endogenousSNPpath(pred, depVar))
  }
  paths <- c(paths,
             mxPath(from = c(depVar), arrows=2, values=1, lbound=0.001,
                    free = !fac, labels = paste(c(depVar), "res", sep = "_")),
             mxPath(depVar, arrows=2, values=0, connect="unique.bivariate"),
             mxPath(from = 'one', to = depVar, free= !fac, values = 0,
                    labels = paste0(depVar, "Mean")))

  dat       <- setupData(phenoData, gxe, force(!missing(minMAF)), minMAF, fitfun)

  modelName <- "OneItem"
  model <- mxModel(model=modelName, type='RAM',
                       manifestVars = manifest,
                       latentVars = latents,
                       paths, dat, makeFitFunction(fitfun))

  postprocessModel(model, depVar, latents)
}

#' @export
#' @importFrom lifecycle deprecate_warn
buildOneItem <- function(phenoData, depVar, covariates=NULL, ..., fitfun = c("WLS","ML"),
                         minMAF=0.01, gxe=NULL, exogenous=NA)
{
  deprecate_warn("0.1.14", "buildOneItem()", "buildItem()")
  buildItem(phenoData, depVar, covariates, fitfun=fitfun,
            minMAF=minMAF, gxe=gxe, exogenous=exogenous)
}

#' Build a model suitable for a single factor genome-wide association study
#'
#' \lifecycle{maturing}
#' The \code{buildOneFac} function is used to specify a single factor latent variable model where the latent variable is predicted by a genomic variant such as a single nucleotide polymorphism, as well as range of covariates. \figure{singleFactor.jpg}{Single Factor Model}
#'
#' @template detail-build
#'
#' @template args-phenoData
#' @param itemNames a character list of the names of the items that load onto the latent variable. These names must match variable names in the phenoData file.
#' @template args-covariates
#' @template args-exogenous
#' @template args-dots-barrier
#' @template args-fitfun
#' @template args-minmaf
#' @template args-gxe
#' @template args-pred
#' @family model builder
#' @export
#' @return
#' \code{buildOneFac} returns an \link[OpenMx:MxModel-class]{MxModel} object that can serve as input for the \link{GWAS} function.
#' @examples
#' pheno <- list()
#' for (i in 1:5) pheno[[paste0('i',i)]] <- rnorm(500)
#' pheno <- as.data.frame(pheno)
#' buildOneFac(pheno, colnames(pheno))
buildOneFac <- function(phenoData, itemNames, covariates=NULL, ..., fitfun = c("WLS","ML"), minMAF=0.01,
			gxe=NULL, exogenous=NA, pred ='snp')
{
  if (length(list(...)) > 0) stop("Rejected are any values passed in the '...' argument")
  fitfun <- match.arg(fitfun)

  fac <- sapply(phenoData[,itemNames,drop=FALSE], is.factor)
  if (is.na(exogenous)) exogenous <- defaultExogenous

  phenoData <- addPlaceholderSNP(phenoData)
  manifest <- c(pred, itemNames)
  if (length(gxe)) manifest <- c(manifest, paste0('snp_', gxe))
  depVar   <- c("F")
  latents  <- depVar
  exoPred <- c()
  endoCovariates <- c()
  if (exogenous) {
	  latents <- c(latents, covariates)
	  exoPred <- covariates
  } else {
	  manifest <- c(manifest, covariates)
	  endoCovariates <- covariates
  }
  if (length(gxe)) pred <- c(pred, paste0('snp_', gxe))
  paths <- c(endogenousSNPpath(pred, depVar),
             endogenousCovariatePaths(phenoData, endoCovariates, depVar))
  paths <- c(paths,
	     mxPath(from=depVar, to=itemNames,values=1, free = T,
		    labels = paste("lambda", itemNames, sep = "_")  ),
	     mxPath(from = c(itemNames), arrows=2, values=1, lbound=0.001, free = !fac,
		    labels = paste(c(itemNames), "res", sep = "_")),
	     mxPath(from=depVar, arrows=2,free=F, values=1.0, labels = "facRes"),
	     mxPath(from = 'one', to = itemNames, free= !fac, values = 0,
		    labels = paste0(itemNames, "Mean")))

  dat       <- setupData(phenoData, gxe, force(!missing(minMAF)), minMAF, fitfun)

  modelName <- "OneFac"
  oneFacPre <- mxModel(model=modelName, type='RAM',
                       manifestVars = manifest,
                       latentVars = latents,
                       paths, dat, makeFitFunction(fitfun))

  postprocessModel(oneFacPre, itemNames, exoPred)
}

#' Build a model suitable for a single factor residual genome-wide association study
#'
#' \lifecycle{maturing}
#' The \code{buildOneFacRes} function is used to specify a single factor latent variable model where a combination of items as well as the latent variable may be predicted by a genomic variant such as a single nucleotide polymorphism, as well as range of covariates. \figure{resid.jpg}{Single Factor Model with a Focus on Residuals}
#'
#' Be aware that a latent variable model is not identified if all of the residuals as well as the latent variable are simultaneously predicted by the SNP.  Specifically, if users wish to use the SNP to predict the latent variable, they much choose at least one (and preferably more that one) item to not be predicted by the SNP.
#'
#' @template detail-build
#'
#' @param factor A logical expression (\code{FALSE} or \code{TRUE}) indicating whether to estimate a regression pathway from the SNP to the latent factor (default FALSE).
#' @param res A character vector of phenotypic item names that indicate which specific items the user wishes to regress on the SNP. The default is to regress all of the items on the SNP.
#' @template args-itemNames
#' @template args-phenoData
#' @template args-covariates
#' @template args-exogenous
#' @template args-fitfun
#' @template args-minmaf
#' @template args-dots-barrier
#' @template args-gxe
#' @template args-pred
#'
#' @family model builder
#' @export
#' @return
#' \code{buildOneFacRes} returns an \link[OpenMx:MxModel-class]{MxModel} object that can serve as input for the \link{GWAS} function.
#' @examples
#' pheno <- list()
#' for (i in 1:5) pheno[[paste0('i',i)]] <- rnorm(500)
#' pheno <- as.data.frame(pheno)
#' buildOneFacRes(pheno, colnames(pheno))
buildOneFacRes <- function(phenoData, itemNames, factor = F, res = itemNames, covariates = NULL,
			   ..., fitfun = c("WLS","ML"), minMAF = .01, gxe=NULL,
			 exogenous=NA,   pred ='snp')
{
  if (length(list(...)) > 0) stop("Rejected are any values passed in the '...' argument")
  fitfun <- match.arg(fitfun)

  fac <- sapply(phenoData[,itemNames,drop=FALSE], is.factor)
  if (is.na(exogenous)) exogenous <- defaultExogenous

  phenoData <- addPlaceholderSNP(phenoData)
  manifest <- c(pred, itemNames)
  if (length(gxe)) manifest <- c(manifest, paste0('snp_', gxe))
  latents   <- c("F")
  depVar <- res
  if (factor) depVar <- c(latents, res)
  exoPred <- c()
  endoCovariates <- c()
  if (exogenous) {
	  latents <- c(latents, covariates)
	  exoPred <- covariates
  } else {
	  manifest <- c(manifest, covariates)
	  endoCovariates <- covariates
  }

  if (length(gxe)) pred <- c(pred, paste0('snp_', gxe))
  paths <- c(endogenousSNPpath(pred, depVar),
             endogenousCovariatePaths(phenoData, endoCovariates, 'F'))
  paths <- c(paths,
	     mxPath(from="F", to=itemNames,values=1, free = T,
		    labels = paste("lambda", itemNames, sep = "_")  ),
	     mxPath(from = c(itemNames), arrows=2, values=1, lbound=0.001, free = c(fac==0),
		    labels = paste(c(itemNames), "res", sep = "_")),
	     mxPath(from="F", arrows=2,free=F, values=1.0, labels = "facRes"),
	     mxPath(from = 'one', to = itemNames, free= c(fac==0), values = 0,
		    labels = paste0(itemNames, "Mean")))

  dat       <- setupData(phenoData, gxe, force(!missing(minMAF)), minMAF, fitfun)

  modelName <- "OneFacRes"
  oneFacPre <- mxModel(model=modelName, type='RAM',
                       manifestVars = manifest,
                       latentVars = latents,
                       paths, dat, makeFitFunction(fitfun))

  postprocessModel(oneFacPre, itemNames, exoPred)
}

#' Build a model suitable for a two factor genome-wide association study
#'
#' \lifecycle{maturing}
#' The buildTwoFac function is used to specify a model with two latent variables where each latent variable is simultaneously predicted by a genomic variant such as a single nucleotide polymorphism, as well as range of covariates. The model allows the latent variables to correlate to accomodate comorbidity between latent traits. \figure{twoFactor.jpg}{Two Factor Model}
#'
#' @template detail-build
#'
#' @template args-F1itemNames
#' @template args-F2itemNames
#' @template args-phenoData
#' @template args-covariates
#' @template args-exogenous
#' @template args-fitfun
#' @template args-minmaf
#' @template args-dots-barrier
#' @template args-gxe
#' @template args-pred
#' @export
#' @family model builder
#' @return
#' \code{buildTwoFac} returns an \link[OpenMx:MxModel-class]{MxModel} object that can serve as input for the \link{GWAS} function.
#' @examples
#' pheno <- list()
#' for (i in 1:10) pheno[[paste0('i',i)]] <- rnorm(500)
#' pheno <- as.data.frame(pheno)
#' buildTwoFac(pheno, paste0('i',1:6), paste0('i',5:10))
buildTwoFac <- function(phenoData, F1itemNames, F2itemNames, covariates = NULL, ...,
			fitfun = c("WLS","ML"), minMAF = .01, gxe=NULL,
			exogenous=NA, pred= 'snp')
{
  if (length(list(...)) > 0) stop("Rejected are any values passed in the '...' argument")
  fitfun <- match.arg(fitfun)

  itemNames <- union(F1itemNames, F2itemNames)

  fac <- sapply(phenoData[,itemNames,drop=FALSE], is.factor)
  if (is.na(exogenous)) exogenous <- defaultExogenous

  phenoData <- addPlaceholderSNP(phenoData)
  manifest <- c(pred, itemNames)
  if (length(gxe)) manifest <- c(manifest, paste0('snp_', gxe))
  depVar <- c("F1", "F2")
  latents   <- depVar
  exoPred <- c()
  endoCovariates <- c()
  if (exogenous) {
	  latents <- c(latents, covariates)
	  exoPred <- covariates
  } else {
	  manifest <- c(manifest, covariates)
	  endoCovariates <- covariates
  }

  if (length(gxe)) pred <- c(pred, paste0('snp_', gxe))
  paths <- c(endogenousSNPpath(pred, depVar),
             endogenousCovariatePaths(phenoData, endoCovariates, depVar))
  paths <- c(paths,
	     mxPath(from="F1", to=F1itemNames,values=1, labels = paste("F1_lambda", F1itemNames, sep = "_")  ),
	     mxPath(from="F2", to=F2itemNames,values=1, labels = paste("F2_lambda", F2itemNames, sep = "_")  ),
	     mxPath(from="F1", to= "F2", arrows=2,free=T, values=.3, labels="facCov"),
	     mxPath(from = itemNames, arrows=2, values=1, lbound=0.001, free = c(fac==0),
		    labels = paste(c(itemNames), "res", sep = "_")),
	     mxPath(from=depVar, arrows=2,free=F, values=1.0, labels = "facRes"),
	     mxPath(from = 'one', to = itemNames, free= c(fac==0), values = 0, labels = paste0(itemNames, "Mean")))

  dat       <- setupData(phenoData, gxe, force(!missing(minMAF)), minMAF, fitfun)

  modelName <- "TwoFac"
  twoFacPre <- mxModel(model=modelName, type='RAM',
                       manifestVars = manifest,
                       latentVars = latents,
                       paths, dat, makeFitFunction(fitfun))

  postprocessModel(twoFacPre, itemNames, exoPred)
}
