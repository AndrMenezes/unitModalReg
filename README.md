
<!-- README.md is generated from README.Rmd. Please edit that file -->

# `unitModalReg`: A Collection of Parametric Modal Regressions Models for Bounded Data

<!-- badges: start -->

<!-- badges: end -->

The goal of `unitModalReg` is to provide tools for fitting parametric
modal regression models for bounded response variables.

## Installation

You can install the development version from
[GitHub](https://github.com/) with:

``` r
# install.packages("remotes")
remotes::install_github("AndrMenezes/unitModalReg")
```

## Example

The packages follows the structure of `glm` objects. The main function
is `unitModalReg`.

``` r
library(unitModalReg)

# Simulate data
n <- 500
betas <- c(beta0 = 2, beta1 = -1, beta2 = 0.6)
x1  <- runif(n)
x2  <- rbinom(n, size = 1, prob = 0.65)
eta <- betas[1] + betas[2] * x1 + betas[3] * x2
mu  <- exp(eta) / exp(1 + eta)
phi <- 4
y   <- rugamma(n = n, mu = mu, phi = phi)
data_sim <- data.frame(y = y, x1 = x1, x2 = factor(x2))

# Families of distributions
distr <- list(
  Be = "betaMode",
  UGz = "unitGompertz",
  UGa = "unitGamma",
  Kum = "Kumaraswamy"
)

# Fit the models
fits <- lapply(distr, function(m) {
  unitModalReg(y ~ x1 + x2, data = data_sim, family = m, link = "logit")
})

# Get the coefficients
vapply(fits, coef, numeric(4))
#>                      Be         UGz         UGa        Kum
#> (Intercept) -0.72578488 -1.54060142 -0.76066727 -0.6221908
#> x1           0.19229626  0.37335286  0.21208397  0.1398764
#> x21          0.06202745  0.04657786  0.04251816  0.0850151
#> phi          1.52605812  0.98750372  3.74129743  2.1089950
```

We also implemented the `gof` function to compare the models based on
the information criterion (AIC, BIC and HQIC) and [Voung
test](https://en.wikipedia.org/wiki/Vuong's_closeness_test) for
non-nested models.

``` r
# Compare the models fitting
gof(lt = fits)
#> 
#> Goodness-of-Fit for modal regression models 
#> 
#> Call:  gof(lt = fits)
#> 
#> Voung tests for all model combinations: 
#>                                 statistic p.value
#> betaMode versus unitGompertz    3.508     0.000  
#> betaMode versus unitGamma       -0.357    0.360  
#> betaMode versus Kumaraswamy     1.385     0.083  
#> unitGompertz versus unitGamma   -6.267    0.000  
#> unitGompertz versus Kumaraswamy -5.313    0.000  
#> unitGamma versus Kumaraswamy    1.103     0.135  
#> 
#> 
#> Information criterions of modal regression models: 
#>      betaMode unitGompertz unitGamma Kumaraswamy
#> aic  -328.821 -203.590     -329.241  -324.864   
#> bic  -311.963 -186.732     -312.382  -308.006   
#> hqic -322.206 -196.975     -322.625  -318.249
```

The currently methods implemented are

``` r
methods(class = "unitModalReg")
#> [1] coef      confint   fitted    hnp       logLik    print     residuals
#> [8] summary   vcov     
#> see '?methods' for accessing help and source code
```
