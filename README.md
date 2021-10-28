# gwsem

<!-- badges: start -->
[![Travis build status](https://api.travis-ci.com/jpritikin/gwsem.svg?branch=master)](https://app.travis-ci.com/github/jpritikin/gwsem)
[![Codecov test coverage](https://codecov.io/gh/jpritikin/gwsem/branch/master/graph/badge.svg)](https://app.codecov.io/gh/jpritikin/gwsem?branch=master)
[![cran version](http://www.r-pkg.org/badges/version/gwsem)](https://cran.r-project.org/package=gwsem)
[![Monthly Downloads](https://cranlogs.r-pkg.org/badges/gwsem)](https://cranlogs.r-pkg.org/badges/gwsem)
[![Total Downloads](https://cranlogs.r-pkg.org/badges/grand-total/gwsem)](https://cranlogs.r-pkg.org/badges/grand-total/gwsem)
<!-- badges: end -->

The goal of gwsem is to provide users with the opportunity to analyze the complex, interconnected array of risk factors, biomarkers, environmental antecedents, comorbid disorders, and other health outcomes on a genome-wide basis using structural equation modeling techniques.

At the moment, plink formats are not supported on the ARM64 architecture. This shortcoming is likely easy to cure once we get access to ARM64 hardware for testing.

## Installation

GW-SEM utilizes the optimization function of OpenMx.

You can install the released version of OpenMx from [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("OpenMx")
```

If you want a new version of OpenMx, you can follow the instruction to build it from source [HERE](https://openmx.ssri.psu.edu).

You can install the released version of gwsem from [CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("gwsem")
```

GW-SEM is currently under development. Therefore, the CRAN version of the software may not include the latest modifactions, enhancements, or functionality.

If you want to use a development snapshot, clone the source code and do

```
git clone https://github.com/jpritikin/gwsem
cd gwsem
./tools/rox
R CMD INSTALL .
Rscript tools/test.R
```

You cannot use **devtools** `install_github` because it does not run **roxygen2**.
