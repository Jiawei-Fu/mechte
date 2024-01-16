#' @title Causal Mediation with BCES estimator
#' @description This function estimates the average mediation effects with the BCES estimator, following the methodology outlined in Fu (2023).. However, it's important to note that we recommend against using this function for your analyses. Instead, we suggest employing the function `simest`, which is likely to yield more reliable and accurate results.
#'
#'@details Total treatment effect on the outcome, denoted by `tau`, is decomposed into `tau = alpha + beta gamma + epsilon`, where `alpha` is the expectation of the average direct effects (can be regarded as the intercept), `epsilon` is the error term, and `beta gamma` is the average mediation effect where `gamma` is the average treatment effect on the mediator and `beta` is a parameter we want to estimate. Function `bces` estimates `beta` and `alpha`.
#'
#' @param gamma_hat a vector of the treatment effect on the mediator.
#' @param tau_hat a vector of the total treatment effect on the outcome.
#' @param sd_u a vector of the standard error of `gamma_hat`, the treatment effect on the mediator.
#' @param sd_v a vector of the standard error of `tau_hat`, the total treatment effect on the outcome.
#' @param confid the size of the test. The default is 0.05.
#' @param bootstrap logical value. Whether to use bootstrap.
#' @param b the number of bootstrap replicates. The default number is 1000.
#' @param seed a single value passed to \code{set.seed}.
#'
#' @return ave_med the estimated average mediation effect.
#' @return beta estimated value of beta.
#' @return alpha estimated value of alpha.
#' @return sd_alpha standard error of the estimated alpha.
#' @return sd_beta standard error of the estimated beta.
#' @return ci_alpha confidence Interval of alpha based on analytic formula.
#' @return ci_beta confidence Interval of beta based on analytic formula.
#' @return boot_ci_alpha confidence Interval of alpha based on bootstrap.
#' @return boot_ci_beta confidence Interval of beta based on bootstrap.
#' @return p_alpha P-value of alpha based on analytic formula.
#' @return p_beta P-value of beta based on analytic formula.
#' @return boot_p_alpha P-value of alpha based on bootstrap.
#' @return boot_p_beta P-value of beta based on bootstrap.
#'
#'
#' @export
#'
#' @examples
#' \dontrun{
#' data(example_dat)  # input data
#' example_dat <- as.data.frame(example_dat)
#' bces(example_dat$gamma_hat,example_dat$tau_hat,example_dat$sd_u,example_dat$sd_v) # simex estimator
#'
#' tmp <- bces(example_dat$gamma_hat,example_dat$tau_hat,sd_u) # save the result
#' tmp$beta # extract the value
#'}
#'
bces <- function(gamma_hat,tau_hat,sd_u,sd_v,confid=0.05, bootstrap=TRUE, b=1000, seed=320){

  n <- length(gamma_hat)

  n1 <- length(gamma_hat)
  n2 <- length(tau_hat)
  n3 <- length(sd_u)
  n4 <- length(sd_v)

  if ( (n1 < 2) |(n2 < 2)|(n3 < 2)|(n4 < 2))
    stop("The length of inputs should be larger than 1")

  if(length(unique(c(n1,n2,n3,n4)))>1)
    stop("The length of inputs should be equal")

  if (!is.null(seed))
    set.seed(seed)

  tau_bar <- mean(tau_hat)
  gamma_bar <- mean(gamma_hat)

  sigma_u <- sd_u^2
  sigma_v <- sd_v^2

  nume <- sum((gamma_hat-gamma_bar)*tau_hat)
  deno <- sum((gamma_hat-gamma_bar)^2) - sum(sigma_u)

  beta_bces <- nume/deno # beta
  alpha_bces <- tau_bar - beta_bces*gamma_bar # alpha

  ##### variance

  xi_1 <- ( (gamma_hat-gamma_bar)*(tau_hat-beta_bces*gamma_hat-alpha_bces) + beta_bces*sigma_u )/(sigma_u- mean(sigma_u))
  zeta_1 <- tau_hat - beta_bces*gamma_hat - gamma_bar*xi_1

  sigma_beta <- mean( (xi_1-mean(xi_1))^2 )
  sigma_alpha <- mean( (zeta_1-mean(zeta_1))^2 )

  sd_beta <- sqrt(sigma_beta/n)
  sd_alpha <- sqrt(sigma_alpha/n)

  pvalue_beta <- 2*pnorm(q=abs(beta_bces)/sd_beta,lower.tail=FALSE)
  pvalue_alpha <- 2*pnorm(q=abs(alpha_bces)/sd_alpha,lower.tail=FALSE)

  ci_beta <- c(beta_bces-qnorm(confid/2,lower.tail=FALSE)*sd_beta,beta_bces+qnorm(confid/2,lower.tail=FALSE)*sd_beta)
  ci_alpha <- c(alpha_bces-qnorm(confid/2,lower.tail=FALSE)*sd_alpha,alpha_bces+qnorm(confid/2,lower.tail=FALSE)*sd_alpha)

  ### bootstrap

  ### percentile-t method

  if(bootstrap==TRUE){

    b_beta_t <- rep(NA,b)
    b_alpha_t <- rep(NA,b)

    ### locally set seed
    old <- .Random.seed
    set.seed(seed)

    for (j in 1:b) {

      b_ind <- sample(1:n,n,TRUE)
      b_gamma_hat <- gamma_hat[b_ind]
      b_tau_hat <- tau_hat[b_ind]

      b_sigma_u <-sigma_u[b_ind]
      b_sigma_v <- sigma_v[b_ind]

      b_tau_bar <- mean(b_tau_hat)
      b_gamma_bar <- mean(b_gamma_hat)

      b_nume <- sum((b_gamma_hat-b_gamma_bar)*b_tau_hat)
      b_deno <- sum((b_gamma_hat-b_gamma_bar)^2) - sum(b_sigma_u)

      b_beta_bces <- b_nume/b_deno  # beta
      b_alpha_bces <- b_tau_bar - b_beta_bces*b_gamma_bar # alpha

      b_xi_1 <- ( (b_gamma_hat-b_gamma_bar)*(b_tau_hat-b_beta_bces*b_gamma_hat-b_alpha_bces) + b_beta_bces*b_sigma_u )/(b_sigma_u- mean(b_sigma_u))
      b_zeta_1 <- b_tau_hat - b_beta_bces*b_gamma_hat - b_gamma_bar*b_xi_1

      b_sigma_beta <- mean( (b_xi_1-mean(b_xi_1))^2 )
      b_sigma_alpha <- mean( (b_zeta_1-mean(b_zeta_1))^2 )

      b_sd_beta <- sqrt(b_sigma_beta/n)
      b_sd_alpha <- sqrt(b_sigma_alpha/n)

      b_beta_t[j] <- (b_beta_bces - beta_bces)/b_sd_beta
      b_alpha_t[j] <- (b_alpha_bces - alpha_bces)/b_sd_alpha
    }

    .Random.seed <<- old  # locally seed

    ### remove NaN

    check_na_beta <- sum(is.na(b_beta_t))
    check_na_alpha <- sum(is.na(b_alpha_t))

    if(check_na_beta+check_na_alpha>0){
      ind_na_beta <- c(1:b)[is.na(b_beta_t)]
      ind_na_alpha <- c(1:b)[is.na(b_alpha_t)]

      ind_na <- union(ind_na_beta,ind_na_alpha)

      b_beta_t <- b_beta_t[-ind_na]
      b_alpha_t <- b_alpha_t[-ind_na]
    }

    #### cont

    quant_beta <- (quantile(b_beta_t, c(confid/2, 1 - confid/2)))
    quant_alpha <- (quantile(b_alpha_t, c(confid/2, 1 - confid/2)))

    ci_beta_low <- beta_bces-sd_beta*quant_beta[[2]]
    ci_beta_upp <- beta_bces-sd_beta*quant_beta[[1]]

    ci_alpha_low <- alpha_bces-sd_alpha*quant_alpha[[2]]
    ci_alpha_upp <- alpha_bces-sd_alpha*quant_alpha[[1]]

    boot_p_beta <-  2 * min(sum((beta_bces/sd_beta) >= b_beta_t)/b, sum((beta_bces/sd_beta) < b_beta_t)/b)

    boot_p_alpha <-  2 * min(sum((alpha_bces/sd_alpha) >= b_alpha_t)/b, sum((alpha_bces/sd_alpha) < b_alpha_t)/b)
  }


  ### output

  output <- matrix(NA, nrow = 2, ncol = 5)
  rownames(output) <- c("(Intercept)","gamma_hat")
  colnames(output) <- c("Estimate","Std. Error" , "CI.low", "CI.up", "P-value")

  output[1,1] <- alpha_bces
  output[2,1] <- beta_bces

  output[1,2] <- sd_alpha
  output[2,2] <- sd_beta

  output[1,3] <- ci_alpha[1]
  output[2,3] <- ci_beta[1]

  output[1,4] <- ci_alpha[2]
  output[2,4] <- ci_beta[2]

  output[1,5] <- pvalue_alpha
  output[2,5] <- pvalue_beta

  cat("The regression outcomes without bootstrap:", "\n")
  printCoefmat(output, P.values = TRUE, has.Pvalue = TRUE)

  if(bootstrap==TRUE){
    output1 <- matrix(NA, nrow = 2, ncol = 4)
    rownames(output1) <- c("(Intercept)","gamma_hat")
    colnames(output1) <- c("Estimate" , "CI.low", "CI.up", "P-value")

    output1[1,1] <- alpha_bces
    output1[2,1] <- beta_bces

    output1[1,2] <- ci_alpha_low
    output1[2,2] <- ci_beta_low

    output1[1,3] <- ci_alpha_upp
    output1[2,3] <- ci_beta_upp

    output1[1,4] <- boot_p_alpha
    output1[2,4] <- boot_p_beta
    cat("\n")
    cat("The regression outcomes with bootstrap:", "\n")
    printCoefmat(output1, P.values = TRUE, has.Pvalue = TRUE)
  }

  ave_med <- mean(beta_bces*gamma_hat)

  cat("\n")
  cat("The average mediation effect is",ave_med,".", "\n")

  # for extract
  output2 <- list()

  output2$beta <- beta_bces
  output2$alpha <- alpha_bces

  output2$sd_alpha <- sd_alpha
  output2$sd_beta <- sd_beta

  output2$ci_alpha <- ci_alpha
  output2$ci_beta <- ci_beta

  output2$p_beta <- pvalue_beta
  output2$p_alpha <- pvalue_alpha

  output2$ave_med <- ave_med

  if(bootstrap==TRUE){
    output2$boot_ci_beta <- c(ci_beta_low,ci_beta_upp)
    output2$boot_ci_alpha <- c(ci_alpha_low,ci_alpha_upp)
    output2$boot_p_beta <-  boot_p_beta
    output2$boot_p_alpha <-  boot_p_alpha
  }

  invisible(output2)
}

