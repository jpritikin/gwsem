#' @param fitfun The fitfun argument specifies which fit function should 
#'         be used in evaluating the GWAS model. Users may choose between 
#'         the relatively rapid "WLS", or the slower but asymptotically 
#'         optimal "ML". In many cases the the differences between the 
#'         fit functions is trivial and the faster "WLS" option should be 
#'         used, but in some situations the differences can be quite 
#'         meaningful (such as when data are Missing at Random - MAR).
