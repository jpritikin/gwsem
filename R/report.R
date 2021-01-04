#' Determine which results are suspicious
#'
#' \lifecycle{maturing}
#' The \link{GWAS} function writes all results, both valid and
#' invalid, to a log file. This function uses heuristics to try to
#' classify rows as suspicious or unsuspicious.
#'
#' @details OpenMx reports exceptions in the \sQuote{catch1}
#'   column. Any error message in the \sQuote{catch1} column is
#'   suspicious. Any optimizer status code besides \sQuote{OK} is
#'   suspicious. It is suspicious if the focal parameter or its
#'   standard error is \code{NA}. If \sQuote{signAdj} was requested
#'   and it is \code{NA} then suspicion is also aroused.
#'
#' @template args-result
#' @param pars names of the parameters available in \code{result}
#' @return
#' a vector of logicals for each row of \code{result} indicating suspicion (if TRUE)
#'
#' @family reporting
#' @export
#' @importFrom stats sd
#' @examples
#' tdir <- tempdir()
#' dir <- system.file("extdata", package = "gwsem")
#' pheno <- data.frame(anxiety=rnorm(500))
#' m1 <- buildItem(pheno, 'anxiety')
#' GWAS(m1, file.path(dir,"example.bgen"),
#'     file.path(tdir,"out.log"))
#' r1 <- loadResults(file.path(tdir,"out.log"), "snp_to_anxiety")
#' r1[isSuspicious(r1, "snp_to_anxiety"),]
isSuspicious <- function(result, pars= attr(result, 'focus')) {
  # minMAF? TODO
  mask <- (!is.na(result[['catch1']]) & result[['catch1']] != '') |
    result[['statusCode']] != 'OK'
  for (p1 in pars) {
    if (is.null(result[[p1]])) {
      warning(paste(p1, "not found; ignoring"))
      next
    }
    mask <- mask | is.na(result[[p1]])
    pvar <- paste0('V',p1,':',p1)
    if (!is.null(result[[pvar]])) {
      mask <- mask | is.na(result[[pvar]])
    }
  }
  mask
}

#' Load GWAS results into a single data.frame
#'
#' \lifecycle{maturing}
#' A1 is the reference allele and A2 is the alternate allele.
#'
#' @param path vector of paths to result files created by \link{GWAS}
#' @param focus parameter name on which to calculate a Z score and p-value
#' @param extraColumns character vector of additional columns to load
#' @param .retainSE logical. Keep a column for the SE of the focus parameter
#' @param signAdj name of column. Value of focus parameter is multiplied by the sign of the named column
#' @param moderatorLevel \lifecycle{deprecated}
#' @template args-dots-barrier
#' @return a data.table with one row per SNP
#' @family reporting
#' @export
#' @importFrom utils read.delim
#' @importFrom data.table fread
#' @importFrom stats pnorm
#' @importFrom lifecycle deprecated deprecate_stop deprecate_soft
#' @examples
#' tdir <- tempdir()
#' dir <- system.file("extdata", package = "gwsem")
#' pheno <- data.frame(anxiety=rnorm(500))
#' m1 <- buildItem(pheno, 'anxiety')
#' GWAS(m1, file.path(dir,"example.bgen"),
#'     file.path(tdir,"out.log"))
#' loadResults(file.path(tdir,"out.log"), "snp2anxiety")
loadResults <- function(path, focus, ..., extraColumns=c(),
			.retainSE=deprecated(), signAdj=deprecated(), moderatorLevel=deprecated()) {

  if (lifecycle::is_present(.retainSE)) {
    deprecate_soft("2.0.6", "loadResults(.retainSE='TRUE')")
  }
  if (lifecycle::is_present(signAdj)) {
    deprecate_stop("2.0.6", "loadResults(signAdj='is moved to the signif() function')")
  }
  if (lifecycle::is_present(moderatorLevel)) {
    deprecate_stop("2.0.6", "loadResults(moderatorLevel='is moved to the signifGxE() function')")
  }

  if (length(list(...)) > 0) stop("Rejected are any values passed in the '...' argument")
  exampleRow <- read.delim(path[1], nrows = 1, check.names=FALSE)
  sel <- c('MxComputeLoop1', 'CHR','BP','SNP','A1','A2','statusCode','catch1',
	   focus, paste0('V',focus,':',focus), extraColumns)
  gxe <- regexpr('^snp_(\\w+)_to_(\\w+)$', focus, perl=TRUE)
  for (vx in 1:length(gxe)) {
    if (gxe[vx] == -1) next
    cstart <- attr(gxe,"capture.start")[vx,]
    clen <- attr(gxe,"capture.length")[vx,]
    mod <- substr(focus[vx], cstart[1], cstart[1]+clen[1]-1L)
    outcome <- substr(focus[vx], cstart[2], cstart[2]+clen[2]-1L)
    mainEffect <- paste0('snp_to_', outcome)
    covName <- c(paste0('V', focus[vx], ':', mainEffect),
                 paste0('V', mainEffect, ':', focus[vx]))
    covName <- covName[covName %in% colnames(exampleRow)]
    if (length(covName) != 1) stop(paste("Can't find parameter covariance between",
                                         mainEffect, 'and', focus[vx]))
    sel <- union(sel, c(mainEffect, paste0('V',mainEffect,':',mainEffect), covName))
 }
 got <- list()
  for (p1 in path) {
    d1 <- fread(p1, header=TRUE,
                sep="\t", check.names=FALSE, quote="", select = sel)
    got <- rbind(got, d1)
  }
  got
}

#' Load suspicious GWAS results into a single data.frame
#'
#' \lifecycle{deprecated}
#' See example for how to change your code.
#'
#' @param path vector of paths to result files created by \link{GWAS}
#' @template args-focus
#' @template args-dots-barrier
#' @param extraColumns character vector of additional columns to load
#' @param signAdj name of column. Value of focus parameter is multiplied by the sign of the named column
#' @param moderatorLevel see details
#' @export
#' @examples
#' tdir <- tempdir()
#' dir <- system.file("extdata", package = "gwsem")
#' pheno <- data.frame(anxiety=rnorm(500))
#' m1 <- buildItem(pheno, 'anxiety')
#' GWAS(m1, file.path(dir,"example.bgen"),
#'     file.path(tdir,"out.log"))
#' r1 <- loadResults(file.path(tdir,"out.log"), "snp_to_anxiety")
#' r1[isSuspicious(r1, "snp_to_anxiety"),]
loadSuspicious <- function(path, focus, ..., extraColumns=c(),
                           signAdj=NULL, moderatorLevel=NULL) {
  deprecate_stop("2.0.6", "gwsem::loadSuspicious()",
                 details="See example for a replacement.")
}

#' Compute Z score and p-value for parameter of focus
#'
#' The \code{signAdj} column is important and not optional for latent
#' factor models.  Loadings to factor indicators can take any sign. If
#' your focus is the regression from the SNP to the factor then this
#' regression estimate will need to be multiplied by the sign of one
#' of the factor loadings. Pick a loading associated with a strong
#' indicator of the factor.
#'
#' Two columns are added, \code{Z} and \code{P}. \code{Z} is the
#' focal parameter divded by its standard error. \code{P} is the
#' unadjusted two-sided normal CDF corresponding to the absolute
#' \code{Z} score.
#'
#' @template args-result
#' @template args-focus
#' @param signAdj name of column. Value of focus parameter is multiplied by the sign of the named column
#' @family reporting
#' @return result with new Z and P columns
#' @export
#' @examples
#' tdir <- tempdir()
#' dir <- system.file("extdata", package = "gwsem")
#' pheno <- data.frame(anxiety=rnorm(500))
#' m1 <- buildItem(pheno, 'anxiety')
#' GWAS(m1, file.path(dir,"example.bgen"),
#'     file.path(tdir,"out.log"))
#' r1 <- loadResults(file.path(tdir,"out.log"), "snp_to_anxiety")
#' r1 <- signif(r1, "snp_to_anxiety")
signif <- function(result, focus, signAdj = NULL) {
  if (length(focus) > 1) stop("Can only do one at a time")
  if (!is.null(signAdj)) {
    if (is.null(result[[signAdj]])) {
      stop(paste('signAdj=', signAdj, 'but no such column found',
                 'in result'))
    }
    result[[focus]] <- result[[focus]] * sign(result[[signAdj]])
  }
  result$Z <- suppressWarnings(
    result[[focus]] / sqrt(result[[paste0('V',focus,':',focus)]]))
  result$P <- 2*pnorm(-abs(result$Z))
  attr(result, 'focus') <- focus
  class(result) <- c("gwsemResult", class(result))
  result
}

#' Compute Z score and p-value for parameter of focus at particular
#' moderator level
#'
#' @template args-result
#' @template args-focus
#' @param level numeric level of the moderator
#' @family reporting
#' @return result with new Z and P columns
#' @export
#' @examples
#' # TODO
signifGxE <- function(result, focus, level) {
  if (length(focus) > 1) stop("Can only do one at a time")
  gxe <- regexpr('^snp_(\\w+)_to_(\\w+)$', focus, perl=TRUE)
  if (gxe == -1) stop("focus must be of the form snp_XXX_to_F")

  cstart <- attr(gxe,"capture.start")
  clen <- attr(gxe,"capture.length")
  mod <- substr(focus, cstart[1], cstart[1]+clen[1]-1L)
  outcome <- substr(focus, cstart[2], cstart[2]+clen[2]-1L)
  mainEffect <- paste0('snp_to_', outcome)
  covName <- c(paste0('V', focus, ':', mainEffect),
               paste0('V', mainEffect, ':', focus))
  covName <- covName[covName %in% colnames(result)]

  tmp <- result[[mainEffect]] + level * result[[focus]]
  tmpSE <- sqrt(result[[paste0('V',mainEffect,':',mainEffect)]] +
                  level^2 * result[[paste0('V',focus,':',focus)]] +
                  2*level * result[[covName]])
  result$Z <- tmp / tmpSE
  result$P <- 2*pnorm(-abs(result$Z))
  attr(result, 'focus') <- focus
  class(result) <- c("gwsemResult", class(result))
  result
}

#' Creates a Manhattan plot
#'
#' Uses the qqman package to create a Manhattan plot.
#'
#' @param x the result of \link{loadResults}
#' @param y an extra argument that should not be used
#' @param ... arguments forwarded to \link[qqman]{manhattan}
#' @export
#' @importFrom qqman manhattan
#' @family reporting
#' @return
#' A Manhattan plot.
#' @examples
#' tdir <- tempdir()
#' dir <- system.file("extdata", package = "gwsem")
#' pheno <- data.frame(anxiety=rnorm(500))
#' m1 <- buildItem(pheno, 'anxiety')
#' GWAS(m1, file.path(dir,"example.bgen"),
#'     file.path(tdir,"out.log"))
#' got <- loadResults(file.path(tdir,"out.log"), "snp_to_anxiety")
#' plot(got)
plot.gwsemResult <- function(x, y, ...) {
	if (!missing(y)) stop("plot does not accept a y= argument")
	x$P[is.na(x$P)] <- 1
	manhattan(x, ...)
}
