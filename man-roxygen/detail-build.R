#' @details
#'
#' This function is designed to be passed to the \link{GWAS} function. The model that is returned, however, is a valid OpenMx model and can be fitted using \link[OpenMx]{mxRun} or \link[OpenMx]{mxTryHard}. Models should be tested to ensure that they are identified and fit the data as expected to avoid unnecessary computation of an invalid model.
#'
#' There is no limit on the number of items that can be included, but more items will exponentially increase computation time.  To address this, we suggest that users use the \sQuote{WLS} fit function. The \sQuote{WLS} fit function is dramatically faster than the \sQuote{ML} fit function, especially for ordinal items.
#' 
#' Ordinal indicator thresholds are setup by
#' \link{setupThresholds}. Exogenous covariates adjustments are setup by
#' \link{setupExogenousCovariates}.
#' You can plot the model using \link[OpenMx]{omxGraphviz}.
