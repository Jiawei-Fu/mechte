## R package for Causal Mediation Analysis

The package contains a suite of functions designed to estimate average causal mediation effects. In general, the method transforms the challenging mediation problem into a simple linear regression problem without compromising the non-parametric nature. More details can be found:

Fu, Jiawei. 2023. "Extract Mechanisms from Heterogeneous Effects: A New Identification Strategy for Mediation Analysis." Working paper, New York University. Available at https://jiaweifu.org/pdf/mechanism.pdf

**(In particular, it is helpful to read section 6.2 and 6.3, two applications.)**

To install and use the latest version of the package, use the following code:
```r
install.packages("devtools")
devtools::install_github("Jiawei-Fu/mechte")
library(mechte)
```

To estimate the average mediation effects using the package, it is assumed that researchers have access to multiple estimated average treatment effects (ATE) for both the mediator and the outcome, along with the corresponding standard errors. These ATEs can be derived using causal trees/forests or through meta-analysis techniques, as discussed in the working paper Fu (2023). For those interested in employing causal forests to find subgroups with significant treatment effects, a useful resource is the tutorial available at https://ml-in-econ.appspot.com/lab3.html. 

Within our package, we provide a pseudo-dataset containing such estimates to facilitate the estimation process. This allows users to familiarize themselves with the package's functionality and to practice the estimation process before applying it to their own data. 

```r
data(example_dat)
example_dat <- as.data.frame(example_dat)
```

The following code snippet demonstrates the implementation of the SIMEX estimator in R. `gamma_hat` is the vector of average treatment effects on the mediator, `tau_hat` is the vector of average treatment effects on the outcome, and `sd_u` is the vector of the standard error of `gamma_hat`. (SIMEX does not need the standard error of `tau_hat`, which is the `sd_v`.) 

```r
simest(example_dat$gamma_hat,example_dat$tau_hat,example_dat$sd_u)
```

We can also visualize the result with following code. The resulting plot will display the observations of average treatment effects on both the mediator and the outcome, their corresponding standard errors, and the estimated regression line. The slope of this regression line represents the key parameter `beta`, which indicates the rate of change in the pure indirect effect relative to a one-unit change in the average treatment effects on the mediator. This parameter is crucial for constructing the average mediation effects.

```r
plot_mech(example_dat$gamma_hat,example_dat$tau_hat,example_dat$sd_u,example_dat$sd_v)
```


