# R-package **`mvs`**

Methods for high-dimensional multi-view learning based on the multi-view
stacking (MVS) framework. Data have a multi-view structure when features
comprise different ‘views’ of the same observations. For example, the
different views may comprise omics, imaging or electronic health
records. Package `mvs` provides functions to fit stacked penalized
logistic regression (StaPLR) models, which are a special case of
multi-view stacking (MVS). Additionally, `mvs` generalizes the StaPLR
model to settings with a Gaussian or Poisson outcome distribution, and
to hierarchical multi-view structures with more than two levels. For
more information about the StaPLR and MVS methods, see Van Loon,
Fokkema, Szabo, & De Rooij (2020) and Van Loon et al. (2022).

## Installation

The current stable release can be installed directly from CRAN:

``` r
utils::install.packages("mvs")
```

The current development version can be installed from GitLab using
package **`devtools`**:

``` r
devtools::install_gitlab("wsvanloon/mvs@develop")
```

## Using **`mvs`**

The two main functions are `StaPLR()` (alias `staplr`), which fits
penalized and stacked penalized regression models models with up to two
levels, and `MVS()` (alias `mvs`), which fits multi-view stacking models
with \>= 2 levels. Objects returned by either function have associated
`coef` and `predict` methods.

### Example: `StaPLR`

``` r
library("mvs")
```

Generate 1000 observations with four two-feature views with varying
within- and between-view correlation:

``` r
set.seed(012)
n <- 1000
cors <- seq(0.1, 0.7, 0.1)
X <- matrix(NA, nrow=n, ncol=length(cors)+1)
X[ , 1] <- rnorm(n)
for (i in 1:length(cors)) {
  X[ , i+1] <- X[ , 1]*cors[i] + rnorm(n, 0, sqrt(1-cors[i]^2))
}
beta <- c(1, 0, 0, 0, 0, 0, 0, 0)
eta <- X %*% beta
p <- exp(eta)/(1+exp(eta))
y <- rbinom(n, 1, p)
```

Fit StaPLR:

``` r
view_index <- rep(1:(ncol(X)/2), each=2)
set.seed(012)
fit <- StaPLR(X, y, view_index)
```

Extract coefficients at the view level:

``` r
coefs <- coef(fit)
coefs$meta
```

    ## 5 x 1 sparse Matrix of class "dgCMatrix"
    ##                    s1
    ## (Intercept) -2.345398
    ## V1           4.693861
    ## V2           .       
    ## V3           .       
    ## V4           .

We see that the only the first view has been selected. The data was
generated so that only the first feature (from the first view) was a
true predictor, but it was also substantially correlated with features
from other views (see `cor(X)`), most strongly with the features from
the fourth view.

Extract coefficients at the base level:

``` r
coefs$base
```

    ## [[1]]
    ## 3 x 1 sparse Matrix of class "dgCMatrix"
    ##                      s1
    ## (Intercept) -0.05351035
    ## V1           0.86273113
    ## V2           0.09756006
    ## 
    ## [[2]]
    ## 3 x 1 sparse Matrix of class "dgCMatrix"
    ##                        s1
    ## (Intercept) -6.402186e-02
    ## V1           1.114585e-38
    ## V2           1.156060e-38
    ## 
    ## [[3]]
    ## 3 x 1 sparse Matrix of class "dgCMatrix"
    ##                      s1
    ## (Intercept) -0.06875322
    ## V1           0.26176566
    ## V2           0.35602028
    ## 
    ## [[4]]
    ## 3 x 1 sparse Matrix of class "dgCMatrix"
    ##                      s1
    ## (Intercept) -0.03101978
    ## V1           0.27605205
    ## V2           0.39234018

We see that the first feature has the strongest effect on the predicted
outcome, with a base-level regression coefficient of 0.86. The features
in views two, three and four all have zero effect, since the meta-level
coefficients for these views are zero.

Compute predictions:

``` r
new_X <- matrix(rnorm(16), nrow=2)
predict(fit, new_X)
```

    ##      lambda.min
    ## [1,]  0.8698197
    ## [2,]  0.1819153

By default, the predictions are made using the values of the penalty
parameters which minimize the cross-validation error (lambda.min).

## Generalizations

-   As StaPLR was developed in the context of binary classification
    problems, the default outcome distribution is `family = "binomial"`.
    Other outcome distributions (e.g., Gaussian, Poisson) can be modeled
    by specifying, e.g., `family = "gaussian"` or `family = "poisson"`.

-   A generalization of stacked penalized (logistic) regression to three
    or more hierarchical levels is implemented in function `MVS` (alias
    `mvs`).

## References

Van Loon, W., De Vos, F., Fokkema, M., Szabo, B., Koini, M., Schmidt,
R., & De Rooij, M. (2022). Analyzing hierarchical multi-view MRI data
with StaPLR: An application to Alzheimer’s disease classification.
*Frontiers in Neuroscience*, *16*, 830630.
<https://doi.org/10.3389/fnins.2022.830630>

Van Loon, W., Fokkema, M., Szabo, B., & De Rooij, M. (2020). Stacked
penalized logistic regression for selecting views in multi-view
learning. *Information Fusion*, *61*, 113–123.
<https://doi.org/10.1016/j.inffus.2020.03.007>
