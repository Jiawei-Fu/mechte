% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simest.R
\name{simest}
\alias{simest}
\title{Causal Mediation with SIMEX estimator}
\usage{
simest(gamma_hat, tau_hat, sd_u, prop = 1, alpha = 0.05, b = 1000)
}
\arguments{
\item{gamma_hat}{a vector of the treatment effect on the mediator.}

\item{tau_hat}{a vector of the total treatment effect on the outcome.}

\item{sd_u}{a vector of the standard error of \code{gamma_hat}, the treatment effect on the mediator.}

\item{prop}{the proportion of each subgroups. The default number 1 means equal proportion.}

\item{alpha}{significant level. The default number is 0.05.}

\item{b}{the number of bootstrap replicates. The default number is 1000.}
}
\value{
\code{acme} the estimated average causal mediation effects.

\code{beta} estimated value of beta.

\code{alpha} estimated value of alpha.

\code{se_alpha} standard error of the estimated alpha.

\code{se_beta} standard error of the estimated beta.

\code{p_alpha} P-value of alpha.

\code{p_beta} P-value of beta.

\code{ci_up_eta} the lower limit of the at least (1-\code{alpha})\% CI

\code{ci_low_eta} the lower limit of the at least (1-\code{alpha})\% CI

\code{Q} Cochran’s Q for heterogeneity test

\code{pvalue_q} pvalue for the Cochran’s Q for heterogeneity test

\code{I_2} Higgins & Thompson’s I^2 for heterogeneity test

\code{dat_sub} the hypotheis test and at least (1-\code{alpha})\% CI for each subgroup
}
\description{
This function estimates the average mediation effects with the SIMEX estimator, following the methodology outlined in Fu (2023). In general, the method transforms the challenging mediation problem into a simple linear regression problem without compromising the non-parametric nature. To implement the SIMEX estimator, we assume that researchers have obtained multiple estimated average treatment effects on the mediator and outcome, and corresponding standard errors. To get those average treatment effects, one can use causal trees/forests or meta-analysis, as described by Fu (2023). With those estimates as inputs, this function returns average mediation effects.
}
\details{
Total treatment effect on the outcome, denoted by \code{tau}, is decomposed into \verb{tau = alpha + beta gamma + epsilon}, where \code{alpha} is the expectation of the average direct effects (can be regarded as the intercept), \code{epsilon} is the error term, and \verb{beta gamma} is the average mediation effect where \code{gamma} is the average treatment effect on the mediator and \code{beta} is a parameter we want to estimate. Function \code{simest} estimates \code{beta} and \code{alpha}.
}
\examples{
data(example_dat)  # input data
example_dat <- as.data.frame(example_dat)
simest(example_dat$gamma_hat,example_dat$tau_hat,example_dat$sd_u) # simex estimator

# assign the result to the tmp
tmp <- simest(example_dat$gamma_hat,example_dat$tau_hat,example_dat$sd_u)
tmp$ave_med  # extract the value
tmp$dat_sub  # extract the subgroup information

}
\references{
Jiawei Fu. 2024. "Extract Mechanisms from Heterogeneous Effects: A New Identification Strategy for Mediation Analysis" \emph{Working Paper}.
}
