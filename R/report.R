#' Load GWAS results into a single data.frame
#'
#' A1 is the reference allele and A2 is the alternate allele.
#' 
#' Two columns are added, \code{Z} and \code{P}. \code{Z} is the
#' focal parameter divded by its standard error. \code{P} is the
#' unadjusted two-sided normal CDF corresponding to the absolute
#' \code{Z} score.
#'
#' @param path vector of paths to result files created by \link{GWAS}
#' @param focus parameter name on which to calculate a Z score and p-value
#' @param extraColumns character vector of additional columns to load
#' @param .retainSE logical. Keep a column for the SE of the focus parameter
#' @template args-dots-barrier
#' @export
#' @importFrom data.table fread
#' @importFrom stats pnorm
#' @examples
#' tdir <- tempdir()
#' dir <- system.file("extdata", package = "gwsem")
#' pheno <- data.frame(anxiety=rnorm(500))
#' m1 <- buildItem(pheno, 'anxiety')
#' GWAS(m1, file.path(dir,"example.pgen"),
#'     file.path(tdir,"out.log"))
#' loadResults(file.path(tdir,"out.log"), "snp2anxiety")
loadResults <- function(path, focus, ..., extraColumns=c(),
			.retainSE=FALSE) {
  sel <- c('MxComputeLoop1', 'CHR','BP','SNP','A1','A2','statusCode','catch1',
	   focus,paste0(focus,'SE'), extraColumns)
  got <- list()
  for (p1 in path) {
    got <- rbind(got, fread(p1, stringsAsFactors = FALSE, header=TRUE,
			    sep="\t", check.names=FALSE, quote="", select = sel))
  }
  got$Z <- got[[focus]] / got[[paste0(focus,'SE')]]
  if (!.retainSE) got[[paste0(focus,'SE')]] <- NULL # redundent; save RAM
  got$P <- 2*pnorm(-abs(got$Z))
  attr(got, 'focus') <- focus
  class(got) <- c("gwsemResult", class(got))
  got
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
#' @return
#' A Manhattan plot.
#' @examples
#' tdir <- tempdir()
#' dir <- system.file("extdata", package = "gwsem")
#' pheno <- data.frame(anxiety=rnorm(500))
#' m1 <- buildItem(pheno, 'anxiety')
#' GWAS(m1, file.path(dir,"example.pgen"),
#'     file.path(tdir,"out.log"))
#' got <- loadResults(file.path(tdir,"out.log"), "snp2anxiety")
#' plot(got)
plot.gwsemResult <- function(x, y, ...) {
	if (!missing(y)) stop("plot does not accept a y= argument")
	x$P[is.na(x$P)] <- 1
	manhattan(x, ...)
}
