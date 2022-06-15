
<!-- README.md is generated from README.Rmd. Please edit that file -->

# converTree

<!-- badges: start -->
<!-- badges: end -->

The goal of converTree is to convert tree format among between newick,
parent vector and ancestor matrix.

## Installation

You can install the development version of converTree like so:

``` r
devtools::install_github("xiayh17/converTree")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(converTree)
## basic example code
```

convert like so:

``` r
## parent vector > child list > nwk
p = c(11, 2, 3, 14, 14, 16, 8, 6, 9, 1, 15, 8, 10, 14, 5, 13, 17)
cl = getChildListFromParentVector(p,16)
childList2nwk(list = cl,n = 16)
## parent vector > nwk
parentVector2nwk(parents=p,n = 16)
## parent vector > anc
parentVector2anc(parents=p,n = 16)
## anc > parent vector
anc <- parentVector2anc(parents=p,n = 16)
anc2parentVector(anc=anc,n=16)
## anc > nwk
anc2nwk(anc=anc,n=16)
```

You’ll still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this. You could also
use GitHub Actions to re-render `README.Rmd` every time you push. An
example workflow can be found here:
<https://github.com/r-lib/actions/tree/v1/examples>.
