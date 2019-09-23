#' Load GWAS results into a single data.frame
#'
#' Two columns are added, \code{Z} and \code{P}. \code{Z} is the
#' focal parameter divded by its standard error. \code{P} is the
#' unadjusted two-sided normal CDF corresponding to the absolute
#' \code{Z} score.
#' 
#' @param path vector of paths to result files created by \link{GWAS}
#' @param focus parameter name on which to calculate a Z score and p-value
#' @param extraColumns character vector of additional columns to load
#' @template args-dots-barrier
#' @export
#' @importFrom data.table fread
#' @importFrom stats pnorm
loadResults <- function(path, focus="snpReg", ..., extraColumns=c()) {
  sel <- c('MxComputeLoop1', 'CHR','BP','SNP','statusCode','catch1',
	   focus,paste0(focus,'SE'), extraColumns)
  got <- list()
  for (p1 in path) {
    got <- rbind(got, fread(p1, stringsAsFactors = FALSE, header=TRUE,
			    sep="\t", check.names=FALSE, quote="", select = sel))
  }
  got$Z <- got[[focus]] / got[[paste0(focus,'SE')]]
  got[[paste0(focus,'SE')]] <- NULL # redundent; save RAM
  got$P <- 2*pnorm(-abs(got$Z))
  got
}
