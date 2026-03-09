#' Summary Methods for mechte Results
#'
#' @param object A `mechte_bces` or `mechte_simest` object.
#' @param ... Additional arguments, currently unused.
#'
#' @return The input object, invisibly.
#' @export
summary.mechte_bces <- function(object, ...) {
  cat("The regression outcomes without bootstrap:\n")
  stats::printCoefmat(object$reg_table, P.values = TRUE, has.Pvalue = TRUE)

  if (isTRUE(object$bootstrap) && !is.null(object$boot_table)) {
    cat("\n")
    cat("The regression outcomes with bootstrap:\n")
    stats::printCoefmat(object$boot_table, P.values = TRUE, has.Pvalue = TRUE)
  }

  cat("\n")
  cat("The average mediation effect is", object$ave_med, ".\n")

  invisible(object)
}

#' @export
summary.mechte_simest <- function(object, ...) {
  cat("The regression outcomes with SIMEX:\n")
  stats::printCoefmat(object$reg_table, P.values = TRUE, has.Pvalue = TRUE)

  cat("\n")
  cat("Heterogenity Test: Cochran’s Q", object$Q, "p-value is", object$pvalue_q, ";\n")
  cat("Heterogenity Test: Higgins & Thompson’s I^2 is", object$I_2 * 100, "%;\n")
  cat("The average causal mediation effect (ACME) is", object$acme, ";\n")
  cat(
    "The at least",
    (1 - object$alpha_level) * 100,
    "% Confidence Interval of ACME is",
    c(object$ci_low_eta, object$ci_up_eta),
    ";\n"
  )
  cat(
    "The test of Null Hypothesis ACME=0 at level",
    object$alpha_level,
    "is",
    object$null_test,
    ", p-value is",
    object$p_value,
    "\n"
  )

  invisible(object)
}
