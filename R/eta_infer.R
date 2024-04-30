eta_infer <- function(gamma,gamma_sd,beta,beta_sd,beta_p,alpha=0.05,prop=1){  
  
  # H0: eta=0 <-> sum gamma_k =0 or beta =0 
  
  if(sum(is.na(gamma))>0){cat("gamma has NA")}
  
  l <- length(gamma)
  if(prop==1){
    prop_new <- rep(1/l,l)
    gamma_new <- prop_new*gamma
    gamma_sd_new <- prop_new*gamma_sd}
  if(prop!=1 & length(prop)!=l){stop("Length of prop is not equal to the length of gamma.")}
  if(prop!=1 & sum(prop)!=l){stop("Sum of the prop is not equal to 1.")}
  if(prop!=1 & sum(prop<=0)>0){stop("Prop has non positive terms.")}
  if(prop!=1){
    gamma_new <- prop*gamma
    gamma_sd_new <- prop*gamma_sd}
  
  z_gamma_mean <- sum(gamma_new,na.rm = T)
  z_gamma_var <- sum(gamma_sd_new^2,na.rm = T)
  
  z_gamma_sd <- sqrt(z_gamma_var)
  
  # p value two sided
  
  p_value_gamma <- 2*pnorm(-abs(z_gamma_mean/z_gamma_sd))
  
  null_test <- NA
  
  if( (p_value_gamma<=alpha) & (beta_p<=alpha) ){
    null_test <- "reject"
  }else{null_test <-  "accept"}
  
  p_value <- min(p_value_gamma,beta_p)  ### our definition of p value under H0 
  
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
  
  
  # output
  eta_tab <- list()
  eta_tab$p_value_gamma <- p_value_gamma
  eta_tab$p_value <- p_value #eta
  eta_tab$ci_up_eta <- ci_up_eta
  eta_tab$ci_low_eta <- ci_low_eta
  eta_tab$eta <- z_gamma_mean*beta
  
  invisible(eta_tab)
  
}

