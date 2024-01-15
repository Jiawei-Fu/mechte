#' @title Causal Mediation with SIMEX estimator
#' @description Estimation of the beta (key part of the causal mediation effect) with the SIMEX estimator.
#'
#' @param gamma_hat a vector of the effect of the treatment on the mediator.
#' @param tau_hat a vector of the total effect of the treatment on the outcome.
#' @param sd_u a vector of the standard error of gamma, the effect of the treatment on the mediator.
#' @param b the number of bootstrap replicates. The default number is 1000.
#'
#' @return beta estimated value of beta.
#' @return alpha estimated value of alpha.
#' @return se_alpha standard error of the estimated alpha.
#' @return se_beta standard error of the estimated beta.
#' @return p_alpha P-value of alpha.
#' @return p_beta P-value of beta.
#' @return ave_med the estimated average mediation effects.
#'
#' @export
#'
#' @examples
#' data(example_dat)  # input data
#' example_dat <- as.data.frame(example_dat)
#' simest(example_dat$gamma_hat,example_dat$tau_hat,sd_u) # simex estimator
#'
#' tmp <- simest(example_dat$gamma_hat,example_dat$tau_hat,example_dat$sd_u) # save the result
#' tmp$beta  # extract the value
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
