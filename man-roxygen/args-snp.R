#' @param SNP a numerical range that specifies the number of SNPs to be evaluated 
#'        from the snpData file. This argument can be used to evaluate a subset of 
#'        snps for model testing. e.g. 1:10 will run the first 10 snps to make sure 
#'        that the model is functioning the way the users intends, that the files 
#'        exist pathways are correct. This option is also very useful to specify a 
#'        range of snps to be evaluated that is smaller than the complete file. For 
#'        example, users may wish to run several discrete batches of analyses for 
#'        chromosome 1, by running 1:10000, 100001:200000, etc. This prevents users 
#'        from constructing numerous snap files for each chromosome. The default 
#'        value of the SNP argument is NULL, which will run all snps in the file.
