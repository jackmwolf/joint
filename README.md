
<!-- README.md is generated from README.Rmd. Please edit that file -->

# joint

<!-- badges: start -->

<!-- badges: end -->

`joint` is an R package designed to efficiently estimate treatment
effects in randomized trials by jointly modeling primary and secondary
efficacy endpoints.

## Installation

You can install the development version of joint from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("jackmwolf/joint")
```

## Example

``` r
library(joint)
```

The example dataset includes $`n = 200`$ observations. Participants were
randomly assigned to the treatment (`A == 1`) or control (`A == 0`)
conditions and four endpoints were measured (`Y1`, `Y2`, `Y3`, and
`Y4`). Categorical versions of `Y3` (`Y3_cat`; binary) and `Y4`
(`Y4_cat`; ordinal with three levels) are also provided.

``` r
data("joint_example")
head(joint_example)
#>   A          Y1         Y2         Y3         Y4 Y3_cat Y4_cat
#> 1 0 -0.41669863 -1.2090084  1.2310544  0.3310951      1      1
#> 2 0  1.51359605 -1.7221864  1.8497439  1.7118287      1      2
#> 3 0  0.41971048  0.7547703 -0.3940902  1.4821349      0      2
#> 4 0  0.73734891  0.3873351 -0.7896880 -0.1314977      0      0
#> 5 0 -0.09441292 -0.4982108 -0.4985990 -2.2440471      0      0
#> 6 0 -1.23037670 -0.3774145 -0.1705403  2.3951030      0      2
```

To estimate the average treatment effect (ATE) on `Y1` while leveraging
data from `Y2` and `Y3`, we could assume that the following one-factor
structural equation model is correctly specified:

``` math
 \begin{pmatrix} Y_1 \\ Y_2 \\ Y_3 \end{pmatrix} |A \sim N\left\{\vec\nu + \gamma\vec\lambda A, \text{diag}(\vec\theta)+\vec\lambda\vec\lambda^T \right\}
```
and estimate the ATE via the maximum likelihood estimator,
$`\widehat\tau_\text{SEM}=\widehat\gamma\widehat{\vec\lambda}`$. This is
accomplished by `joint_sem()`:

``` r
fit_sem <- joint_sem(
  data0 = joint_example,
  endpoints = c("Y1", "Y2", "Y3"),
  treatment = "A"
)
```

ATE estimates along with variance estimates can be extracted using
`estimate_effects()`:

``` r
estimate_effects(fit_sem)
#> $estimate
#>         Y1         Y2         Y3 
#>  0.4505280  0.4511448 -0.4552091 
#> 
#> $vcov
#>              Y1           Y2           Y3
#> Y1  0.021643736  0.001822643 -0.001686521
#> Y2  0.001822643  0.022190595 -0.001502853
#> Y3 -0.001686521 -0.001502853  0.023348443
```

However, we suspect that the model may be misspecified, we may wish to
use model averaging and additionally consider a saturated means model
(which would result in unbiased and consistent treatment effect
estimation) and the estimator
$`\widehat\tau_\text{MA}(\omega)=\omega\widehat\tau_\text{SEM}+(1-\omega)\widehat\tau_\text{Saturated}`$
for $`\omega \in [0,1]`$. We could then select $`\omega`$ using super
learning (`joint_sl()`) or Bayesian model averaging (`joint_bma()`):

``` r
set.seed(1)
fit_sl <- joint_sl(
  data0 = joint_example,
  endpoints = c("Y1", "Y2", "Y3"),
  treatment = "A",
  primary = "Y1",
  n_boot = 5 # note the choice of n_boot = 5; more are recommended in practice
)
fit_sl$summary
#>      Method  Estimate        se   se_boot   lb_boot   ub_boot Omega
#> 1    SL-SEM 0.4505280        NA 0.1022004 0.3598210 0.5896298     1
#> 2       SEM 0.4505280 0.1471181 0.1143030 0.3649402 0.6265119     1
#> 3 Saturated 0.4299304 0.1549339 0.1176454 0.3063521 0.5845191     0
```

## Future Work

- Support parallel computation for `joint_sl` and `joint_bma`

- Develop `summary` methods for `joint_sl` and `joint_bma` output

## References

- Wolf, J. M., Koopmeiners, J. S., & Vock, D. M. (2025). Jointly
  modeling multiple endpoints for efficient treatment effect estimation
  in randomized controlled trials. arXiv preprint arXiv:2506.03393
