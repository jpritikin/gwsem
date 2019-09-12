#' Load GWAS results into a single data.frame
#' 
#' @param path vector of paths to result files created by \link{GWAS}
#' @export
#' @importFrom utils read.table
loadResults <- function(path) {
  got <- list()
  for (p1 in path) {
    got <- rbind(got, read.table(p1, stringsAsFactors = FALSE, header=TRUE,
				 sep="\t", check.names=FALSE, quote="", comment.char=""))
  }
  got
}
