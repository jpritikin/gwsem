# gwsem

<!-- badges: start -->
[![Travis build status](https://travis-ci.org/jpritikin/gwsem.svg?branch=master)](https://travis-ci.org/jpritikin/gwsem)
[![Codecov test coverage](https://codecov.io/gh/jpritikin/gwsem/branch/master/graph/badge.svg)](https://codecov.io/gh/jpritikin/gwsem?branch=master)
[![cran version](http://www.r-pkg.org/badges/version/gwsem)](https://cran.r-project.org/package=gwsem)
[![Monthly Downloads](http://cranlogs.r-pkg.org/badges/gwsem)](http://cranlogs.r-pkg.org/badges/gwsem)
[![Total Downloads](http://cranlogs.r-pkg.org/badges/grand-total/gwsem)](http://cranlogs.r-pkg.org/badges/grand-total/gwsem)
<!-- badges: end -->

The goal of gwsem is to provide users with the opportunity to to analyze the complex, interconnected array of risk factors, biomarkers, environmental antecedents, comorbid disorders and other health outcomes on a genome-wide basis using structural equation modeling techniques.

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
