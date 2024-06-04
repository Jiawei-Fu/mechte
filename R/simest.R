#' @title Causal Mediation with SIMEX estimator
#' @description This function estimates the average mediation effects with the SIMEX estimator, following the methodology outlined in Fu (2023). In general, the method transforms the challenging mediation problem into a simple linear regression problem without compromising the non-parametric nature. To implement the SIMEX estimator, we assume that researchers have obtained multiple estimated average treatment effects on the mediator and outcome, and corresponding standard errors. To get those average treatment effects, one can use causal trees/forests or meta-analysis, as described by Fu (2023). With those estimates as inputs, this function returns average mediation effects.
#'
#' @details Total treatment effect on the outcome, denoted by `tau`, is decomposed into `tau = alpha + beta gamma + epsilon`, where `alpha` is the expectation of the average direct effects (can be regarded as the intercept), `epsilon` is the error term, and `beta gamma` is the average mediation effect where `gamma` is the average treatment effect on the mediator and `beta` is a parameter we want to estimate. Function `simest` estimates `beta` and `alpha`.
#'
#' @param gamma_hat a vector of the treatment effect on the mediator.
#' @param tau_hat a vector of the total treatment effect on the outcome.
#' @param sd_u a vector of the standard error of `gamma_hat`, the treatment effect on the mediator.
#' @param prop the proportion of each subgroups. The default number 1 means equal proportion.
#' @param b the number of bootstrap replicates. The default number is 1000.
#' @param alpha significant level. The default number is 0.05.
#'
#' @returns `acme` the estimated average causal mediation effects.
#' @returns `beta` estimated value of beta.
#' @returns `alpha` estimated value of alpha.
#' @returns `se_alpha` standard error of the estimated alpha.
#' @returns `se_beta` standard error of the estimated beta.
#' @returns `p_alpha` P-value of alpha.
#' @returns `p_beta` P-value of beta.
#' @returns `ci_up_eta` the lower limit of the at least (1-`alpha`)% CI
#' @returns `ci_low_eta` the lower limit of the at least (1-`alpha`)% CI
#' @returns `Q` Cochran’s Q for heterogeneity test
#' @returns `pvalue_q` pvalue for the Cochran's Q for heterogeneity test
#' @returns `I_2` Higgins & Thompson's I^2 for heterogeneity test
#' @returns `dat_sub` the hypotheis test and at least (1-`alpha`)% CI for each subgroup
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
#' tmp$dat_sub  # extract the subgroup information
#'
#' @references Jiawei Fu. 2024. "Extract Mechanisms from Heterogeneous Effects: A New Identification Strategy for Mediation Analysis" \emph{Working Paper}.
simest <- function(gamma_hat,tau_hat,sd_u,X=NULL,prop=1,alpha=0.05,b=1000){

  if(sum(is.na(gamma_hat))>0){stop("gamma has NA")}
  if(sum(is.na(sd_u))>0){stop("sd_u has NA")}
  if(sum(is.na(tau_hat))>0){stop("tau_hat has NA")}


  ### remove na

  n1 <- length(gamma_hat)
  n2 <- length(tau_hat)
  n3 <- length(sd_u)

  if ( (n1 < 2) |(n2 < 2)|(n3 < 2))
    stop("The length of inputs should be larger than 1")

  if(length(unique(c(n1,n2,n3)))>1)
    stop("The length of inputs should be equal")


  ###### HETEROGENEITY for gamma

  w <- 1/(sd_u)^2
  theta_fix <- sum(w*gamma_hat)/sum(w)
  Q <- sum((w)*(gamma_hat-theta_fix)^2)

  pvalue_q <- pchisq(Q, n1-1, ncp = 0, lower.tail=FALSE)

  tmp_i <- NA
  if(Q-n1+1<0){tmp_i <- 0}else{tmp_i <- Q-n1+1}

  I_2 <- tmp_i/Q

  #### ESTIMATION

  if(!is.null(X)){
    lm_formula <- as.formula(paste("tau_hat", paste(c("gamma_hat",names(X)), collapse = ' + '), sep = " ~ "))
    mod <- lm(lm_formula,data=X, x=TRUE,y=TRUE)
  }else{
    mod <- lm(tau_hat~gamma_hat,x=TRUE,y=TRUE)
  }


  mod_sim <- simex::simex(mod,B=b,
                   measurement.error = sd_u,
                   SIMEXvariable="gamma_hat",fitting.method ="quad",asymptotic="FALSE")

  beta_p <- summary(mod_sim)$coefficients$jackknife["gamma_hat",4]

  beta_sd <- summary(mod_sim)$coefficients$jackknife["gamma_hat",2]
  beta <- mod_sim$coefficients[2]

  ### inference

  if(length(prop)==1 & prop[1]==1){
    prop_new <- rep(1/n1,n1)
    gamma_new <- prop_new*gamma_hat
    gamma_sd_new <- prop_new*sd_u}
  if(length(prop)==1 & prop[1]!=1){stop("The proportion has one value and is not 1; it should be a vetor or 1 ")}
  if(length(prop)!=1 & length(prop)!=n1){stop("Length of prop is not equal to the length of gamma.")}
  if(length(prop)!=1 & sum(prop)!=n1){cat("Sum of the prop is not equal to 1. \n")}
  if(length(prop)!=1 & sum(prop<=0)>0){stop("Prop has non positive terms.")}
  if(length(prop)!=1 ){
    gamma_new <- prop*gamma_hat
    gamma_sd_new <- prop*sd_u}

  z_gamma_mean <- sum(gamma_new,na.rm = T)
  z_gamma_var <- sum(gamma_sd_new^2,na.rm = T)
  z_gamma_sd <- sqrt(z_gamma_var)

  # p value two sided

  p_value_gamma <- 2*pnorm(-abs(z_gamma_mean/z_gamma_sd))

  null_test <- NA

  cat(p_value_gamma,beta_p)

  if( (p_value_gamma<=alpha) & (beta_p<=alpha) ){
    null_test <- "rejected"
  }else{null_test <-  "not rejected"}

  p_value <- max(p_value_gamma,beta_p)  ### our definition of p value under H0

  # 1-alpha interval
  # we need sqrt(1-alpha) interval for gamma and beta

  ci_up_gamma <- z_gamma_mean + qnorm((1+sqrt(1-alpha))/2)*z_gamma_sd
  ci_low_gamma <- z_gamma_mean - qnorm((1+sqrt(1-alpha))/2)*z_gamma_sd

  ci_up_beta <- beta + qnorm((1+sqrt(1-alpha))/2)*beta_sd
  ci_low_beta <- beta - qnorm((1+sqrt(1-alpha))/2)*beta_sd

  ### CI for eta

  ci_up_eta <- ci_low_eta <-  NA

  tmp_a <- ci_up_gamma*ci_up_beta
  tmp_b <- ci_up_gamma*ci_low_beta
  tmp_c <- ci_low_gamma*ci_up_beta
  tmp_d <- ci_low_gamma*ci_low_beta

  ci_up_eta <- max(tmp_a,tmp_b,tmp_c,tmp_d)
  ci_low_eta <- min(tmp_a,tmp_b,tmp_c,tmp_d)

  ave_med <- z_gamma_mean*beta


  ##### CI and hypothesis test for each group

  rec_sub_test <- rep(NA,n1)
  p_value_subgroup <- 2*pnorm(-abs(gamma_hat/sd_u))

  for (i in 1:n1) {
    if( (p_value_subgroup[i]<=alpha) & (beta_p<=alpha) ){
      rec_sub_test[i] <- "rejected"
    }else{rec_sub_test[i] <-  "not rejected"}
  }


  dat_sub <- data.frame("group" = 1:n1,
                        "gamma" = gamma_hat,
                        "test" = rec_sub_test,
                        "p_value" = NA,
                        "ci_up" = NA,
                        "ci_low" =NA)

  ci_up_subgamma <- gamma_hat + qnorm((1+sqrt(1-alpha))/2)*sd_u
  ci_low_subgamma <- gamma_hat - qnorm((1+sqrt(1-alpha))/2)*sd_u

  tmp_sub_a <- ci_up_subgamma*ci_up_beta
  tmp_sub_b <- ci_up_subgamma*ci_low_beta
  tmp_sub_c <- ci_low_subgamma*ci_up_beta
  tmp_sub_d <- ci_low_subgamma*ci_low_beta

  for (i in 1:n1) {
    dat_sub$ci_up[i] <- max(tmp_sub_a[i],tmp_sub_b[i],tmp_sub_c[i],tmp_sub_d[i])
    dat_sub$ci_low[i] <- min(tmp_sub_a[i],tmp_sub_b[i],tmp_sub_c[i],tmp_sub_d[i])
    dat_sub$p_value[i] <- max(p_value_subgroup[i],beta_p)
  }



  ### output

  cat("The regression outcomes with SIMEX:", "\n")
  printCoefmat(summary(mod_sim)$coefficients$jackknife, P.values = TRUE, has.Pvalue = TRUE)

  cat("\n")
  cat("Heterogenity Test: Cochran’s Q",Q, "p-value is",pvalue_q, "; \n")
  cat("Heterogenity Test: Higgins & Thompson’s I^2 is",I_2*100, "%; \n")

  cat("The average causal mediation effect (ACME) is",ave_med, "; \n")
  cat("The at least",(1-alpha)*100,"% Confidence Interval of ACME is",c(ci_low_eta,ci_up_eta), "; \n")
  cat("The test of Null Hypothesis ACME=0 at level",alpha,"is",null_test,", p-value is",p_value, " \n")


  # for extract

  output2 <- list()

  output2$beta <- mod_sim$coefficients[2]
  output2$alpha <- mod_sim$coefficients[1]

  output2$se_alpha <- summary(mod_sim)$coefficients$jackknife[1,2]
  output2$se_beta <- summary(mod_sim)$coefficients$jackknife[2,2]

  output2$p_beta <- summary(mod_sim)$coefficients$jackknife[2,4]
  output2$p_alpha <- summary(mod_sim)$coefficients$jackknife[1,4]

  output2$acme <- ave_med

  output2$p_value_gamma <- p_value_gamma
  output2$p_value <- p_value #eta
  output2$ci_up_eta <- ci_up_eta
  output2$ci_low_eta <- ci_low_eta

  output2$Q <- Q
  output2$pvalue_q <- pvalue_q
  output2$I_2 <- I_2


  output2$dat_sub <- dat_sub

  invisible(output2)
}
