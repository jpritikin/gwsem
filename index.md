# Welcome to the alternate GW SEM documentation site!

***
Here we demonstrate how to use the relevant functions to conduct GWAS analyses.  Because GWAS data cannot be easily shared on the internet, the demonstrations use simulated data.  In some cases, the data are loosely based on the analyses that we conducted for publication. To keep files sizes manageable, demonstrations consist of 6,000 "individuals" and 2,000 SNPs.  While we know that this is a laughably small sample for a real GWAS study, it works well to demonstrate how the functions work.  Suggestions for analyzing real data are made throughout the tutorial to help bridge the gap between didactic and practical applications.

***

## Installation

You can install the released version of gwsem from [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("gwsem")
```

GW-SEM is currently under intense development. Therefore, the CRAN version of the software will not include the latest modifactions, enhancements, enhancements or functionality.

If you want to use a development snapshot, clone the source code and do

```
git clone https://github.com/jpritikin/gwsem
cd gwsem 
./tools/rox
cd ..
R CMD INSTALL gwsem R_LIBS_USER
Rscript tools/test.R
```

You cannot use **devtools** `install_github` because it does not run **roxygen2**.

GW-SEM utilize the optimization function of OpenMx. Therefore to use GW-SEM you must have an up-to-date version of OpenMx (Newer than version X.X-XX)

You can install the released version of OpenMx from [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("OpenMx")
```

If you want a new version of OpenMx, you can follow the instruction to build it from source [HERE](https://openmx.ssri.psu.edu).

***

The goal of GW-SEM is to provide users with the opportunity to to analyze the complex, interconnected array of risk factors, biomarkers, environmental antecedents, comorbid disorders and other health outcomes on a genome-wide basis using structural equation modeling techniques.

GW-SEM is a computationally efficient, flexible and accessible algorithm to conduct multivariate genome-wide association studies within a structural equation modeling framework.  In the current release of GW-SEM, we have enhanced the computational efficiency of the algorithm, allowing users to select a fit function that is appropriate for their analyses, expanded compatibility with standard genomic data formats, and formatted the output to be read seamlessly into other standard post-GWAS processing software, and streamlining the requirements for user's to specify custom GWAS models. 

GW-SEM analyses provide substantially deeper insights into the underlying genomic pathway than is possible with standard GWAS methods.

Tutorials:

- [Standard GWAS](1.-Standard-GWAS.html)
- [One factor model](2.-One-Factor-Model.html)
- [One factor residuals model](3.-Residuals-Model.html)
- [Two factor model](4.-Two-Factor-Model.html)
- [User specified custom GWAS](5.-User-Specified-GWAS-Models.html)
- [Gene environment interaction](6.-Gene-Environment-Interaction--GxE--Models.html)
- [Post GWAS processing functions](7.-Post-GWAS-processing-functions.html)
- [Latent growth curve](growth.html)
