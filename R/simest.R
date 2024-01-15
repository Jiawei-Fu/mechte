#' @title Causal Mediation with SIMEX estimator
#' @description This function estimates the average mediation effects with the SIMEX estimator, following the methodology outlined in Fu (2023). In general, the method transforms the challenging mediation problem into a simple linear regression problem without compromising the non-parametric nature. To implement the SIMEX estimator, we assume that researchers have obtained multiple estimated average treatment effects on the mediator and outcome, and corresponding standard errors. To get those average treatment effects, one can use causal trees or meta-analysis, as described by Fu (2023). With those estimates as inputs, this function returns average mediation effects.
#'
#' @param gamma_hat a vector of the effect of the treatment on the mediator.
#' @param tau_hat a vector of the total effect of the treatment on the outcome.
#' @param sd_u a vector of the standard error of gamma, the effect of the treatment on the mediator.
#' @param b the number of bootstrap replicates. The default number is 1000.
#'
#' @returns `ave_med` the estimated average mediation effects.
#' @returns `beta` estimated value of beta.
#' @returns `alpha` estimated value of alpha.
#' @returns `se_alpha` standard error of the estimated alpha.
#' @returns `se_beta` standard error of the estimated beta.
#' @returns `p_alpha` P-value of alpha.
#' @returns `p_beta` P-value of beta.
#'
#'
#' @export
#'
#' @examples
#' data(example_dat)  # input data
#' example_dat <- as.data.frame(example_dat)
#' simest(example_dat$gamma_hat,example_dat$tau_hat,example_dat$sd_u) # simex estimator
#'
#' # assign the result to the tmp
#' tmp <- simest(example_dat$gamma_hat,example_dat$tau_hat,example_dat$sd_u)
#' tmp$ave_med  # extract the value
#'
#' @references Jiawei Fu. 2013. "Extract Mechanisms from Heterogeneous Effects: A New Identification Strategy for Mediation Analysis" \emph{Working Paper}.
simest <- function(gamma_hat,tau_hat,sd_u,b=1000){

  n1 <- length(gamma_hat)
  n2 <- length(tau_hat)
  n3 <- length(sd_u)

  if ( (n1 < 2) |(n2 < 2)|(n3 < 2))
    stop("The length of inputs should be larger than 1")

  if(length(unique(c(n1,n2,n3)))>1)
    stop("The length of inputs should be equal")

  mod <- lm(tau_hat~gamma_hat,x=TRUE,y=TRUE)

  mod_sim <- simex::simex(mod,B=b,
                   measurement.error = sd_u,
                   SIMEXvariable="gamma_hat",fitting.method ="quad",asymptotic="FALSE")

  ### output

  cat("The regression outcomes with SIMEX:", "\n")
  printCoefmat(summary(mod_sim)$coefficients$jackknife, P.values = TRUE, has.Pvalue = TRUE)

  tmp_beta <- mod_sim$coefficients[2]
  ave_med <- mean(tmp_beta*gamma_hat)

  cat("\n")
  cat("The average mediation effect is",ave_med, "\n")

  # for extract

  output2 <- list()

  output2$beta <- mod_sim$coefficients[2]
  output2$alpha <- mod_sim$coefficients[1]

  output2$se_alpha <- summary(mod_sim)$coefficients$jackknife[1,2]
  output2$se_beta <- summary(mod_sim)$coefficients$jackknife[2,2]

  output2$p_beta <- summary(mod_sim)$coefficients$jackknife[2,4]
  output2$p_alpha <- summary(mod_sim)$coefficients$jackknife[1,4]

  output2$ave_med <- ave_med

  invisible(output2)
}
