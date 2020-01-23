#' Genome-wide Structural Equation Modeling
#'
#' @docType package
#' @name gwsem-package
#' @aliases gwsem
#'
#' @description
#'
#' Psychometricians have long known that little information can be
#' gleened from a single item. Hence, there is a long-standing
#' tradition in education, psychology, and many other fields to use
#' more than one item to measure a latent trait. For example, a math
#' test will always consist of more than one problem (or the single
#' problem with consist of many parts).
#'
#' Phenotypic data gathered at the same time as genetic data sometimes
#' contains multiple items that measure different aspects of the same
#' latent construct. However, due to the astronomic number of single
#' nucleotide polymorphisms (SNPs) to test, fast analysis methods are
#' generally preferred with much of the prior research employing
#' regression. Regression is fast, but can only predict a single item
#' at time. Hence, associations with rich phenotypic data cannot be
#' properly investigated.
#'
#' \pkg{gwsem} contains low-level C/C++ code to permit OpenMx to
#' rapidly read genetic data encoded in U.K. Biobank or plink formats.
#' The association between SNPs and a factor model can be explored.
#'
#' @useDynLib gwsem, .registration = TRUE
#' @import OpenMx
#' @importFrom Rcpp evalCpp
#' 
NULL

.onAttach <- function(libname, pkgname) {
	if (.Platform$OS.type == "windows" && .Machine$sizeof.pointer == 4) {
		packageStartupMessage("C++ exceptions do not work correctly on 32-bit Windows. You should use GNU/Linux for any serious work.")
	}
}
