#' @param exogenous This argument specifies how you would like to integrate 
#'        the covariates into the analysis. If exogenous = T, each items will 
#'        be directly regressed on each covariate. If exogenous = F, the latent 
#'        factor(s) will be directly regressed on each covariate. Setting 
#'        exogenous = T does not assume that the items are related to the 
#'        covariates proportional to their factor loadings (which is probably 
#'        preferable in most cases).
