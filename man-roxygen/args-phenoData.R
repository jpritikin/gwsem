#' @param phenoData the file pathway for the phenotypic data (e.g. "myData.txt" or 
#'         "phenotype/myData.txt"). This data file can include more variables than 
#'         those included in the analysis, but GW-SEM will only use the items/covariates 
#'         that are specified. (The dangers of very large dataset is that they can take 
#'         a long time to load and can take up space in the R environment. This should 
#'         not affect processing speed for the GWAS analysis, but can create headaches 
#'         for pre-processing).
