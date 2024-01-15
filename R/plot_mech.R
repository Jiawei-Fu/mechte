#' Plot the estimates by
#'
#' @param gamma_hat a vector of the effect of the treatment on the mediator.
#' @param tau_hat a vector of the total effect of the treatment on the outcome.
#' @param sd_u a vector of the standard error of gamma, the effect of the treatment on the mediator.
#' @param sd_v a vector of the standard error of tau, the total effect of the treatment on the outcome.
#' @param confid the size of the test. The default is 0.05.
#' @param b the number of bootstrap replicates. The default number is 1000.
#' @param fit the estimator to use. Default is "SIMEX"; additional options are "BCES" or "both".
#' @param seed a single value passed to \code{set.seed}.
#' @param legend logical value; whether to include legend in the figure. Default is "TRUE".
#' @param xlab the label on the x-axis. Default is "ATE on the Mediator".
#' @param ylab the label on the y-axis. Default is "ATE on the Outcome".
#' @param col_error_bar the color of the error bar. Default isn "#219ebc".
#' @param col_point the color of the point. Default isn "#219ebc".
#' @param col_line the color of the fitted line. Default is "#8ac926" and "#fe7f2d".
#' @param size_title the size of the title. Default is 14.
#' @param size_text the size of the text on the axis. Default is 12.
#' @param size_legend the size of the legend. Default is 12.
#'
#' @return Figure.
#' @export
#'
#' @examples
#'\dontrun{
#' data(example_dat)  # input data
#' example_dat <- as.data.frame(example_dat)
#' plot_mech(example_dat$gamma_hat,example_dat$tau_hat,example_dat$sd_u,example_dat$sd_v)
#'}
plot_mech <- function(gamma_hat,tau_hat,sd_u,sd_v,confid=0.05,b=1000,fit="SIMEX",seed=320, legend="TRUE",xlab="ATE on the Mediator",ylab="ATE on the Outcome",col_error_bar="#219ebc",col_point="#219ebc",col_line=c("#8ac926","#fe7f2d"),size_title=14,size_text=12,size_legend=12){

  gamma_hat <- gamma_hat
  tau_hat <- tau_hat
  sd_u <- sd_u
  sd_v <- sd_v

  dat <- data.frame(gamma_hat,tau_hat,sd_u,sd_v)

  mod <- bces(gamma_hat,tau_hat,sd_u,sd_v,confid=confid, bootstrap=TRUE, b=b, seed=seed)

  mod_sim <- simex::simex(lm(tau_hat~gamma_hat,x=TRUE,y=TRUE),B=b,
                   measurement.error = sd_u,
                   SIMEXvariable="gamma_hat",fitting.method ="quad",asymptotic="FALSE")

  beta_bces <- mod$beta_bces
  alpha_bces <- mod$alpha_bces

  beta_sim <- summary(mod_sim)[[1]][[1]][2,1]
  alpha_sim <- summary(mod_sim)[[1]][[1]][1,1]


  if(fit=="both"){
    dat_est <- data.frame(alpha=c(alpha_bces,alpha_sim),beta=c(beta_bces,beta_sim),Estimators=c("BCES","SIMEX"))
  }

  if(fit=="BCES"){
    dat_est <- data.frame(alpha=alpha_bces,beta=beta_bces,Estimators="BCES")
  }

  if(fit=="SIMEX"){
    dat_est <- data.frame(alpha=alpha_sim,beta=beta_sim,Estimators="SIMEX")
  }

  ggplot2::ggplot(data=dat, ggplot2::aes(x=gamma_hat, y=tau_hat)) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin=tau_hat-qnorm(confid/2,lower.tail=FALSE)*sd_v, ymax=tau_hat+qnorm(confid/2,lower.tail=FALSE)*sd_v), colour=col_error_bar, width=0,alpha=0.5) +
    ggplot2::geom_errorbarh(ggplot2::aes(xmin=gamma_hat-qnorm(confid/2,lower.tail=FALSE)*sd_u, xmax=gamma_hat+qnorm(confid/2,lower.tail=FALSE)*sd_u), colour=col_error_bar, height=0,alpha=0.5) +
    ggplot2::geom_point(alpha=1,colour=col_point) +
    ggplot2::geom_abline(data=dat_est,ggplot2::aes(intercept=alpha, slope=beta,colour=Estimators),size=1, show.legend=legend)+
    #    ggplot2::geom_abline(dat_est,ggplot2::aes(intercept=alpha_sim, slope=beta_sim),colour=group,size=1, show.legend=TRUE)+
    ggplot2::scale_color_manual(values = col_line )+
    ggplot2::labs(x=xlab, y=ylab)+
    ggplot2::theme_light()+
    ggplot2::theme(axis.text=element_text(size=size_text),
                   axis.title=element_text(size=size_title,face="bold"))+
    ggplot2::theme(legend.text=element_text(size=size_legend))

}
