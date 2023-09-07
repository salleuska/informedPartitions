### Info

The scripts in this folder reproduces plots presented in Section 3.2 of the paper and in the Appendix. T

### Preliminaries


In addition to libraries cited in the main folder, make sure to have the following `R` packages installed.

```r
install.packages("devtools")  ## to install packages from github

install.packages("salso")

## packages for output processing and visualization
install.packages("reshape2")
install.packages("ggplot2")
install.packages("grid")
install.packages("gridExtra")
install.packages("latex2exp")

## to evaluate LSP prior
install.packages("Rcpp")
install.packages("RcppArmadillo")
## to evaluate iCRP prior
install.packages("betafunctions")

install.packages("mclust") ## adjustedRandIndex
install.packages("mcclust") ## for variation of information

#install.packages("bayestestR")
#install.packages("cowplot")


```
The scripts make use of two `R` packages available on Github.

```r
library(devtools)

devtools::install_github("gpage2990/drpm") ## package for dependent random partition models
devtools::install_github("salleuska/CPLogit", subdir = "CPLogit") ## package for centered partition process regression model

```
