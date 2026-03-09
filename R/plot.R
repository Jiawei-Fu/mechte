#' Plot Method for SIMEX Mediation Results
#'
#' @param x A `mechte_simest` object returned by `simest()`.
#' @param sd_v Optional standard errors for `tau_hat`. If supplied, vertical
#'   error bars are added to match `plot_mech()`.
#' @param confid The size of the test. Defaults to the `alpha` used in `simest()`.
#' @param legend Logical; whether to include the legend.
#' @param xlab The x-axis label.
#' @param ylab The y-axis label.
#' @param col_error_bar The color of the error bars.
#' @param col_point The color of the points.
#' @param col_line The color of the fitted SIMEX line.
#' @param size_title The size of the axis titles.
#' @param size_text The size of the axis text.
#' @param size_legend The size of the legend text.
#' @param ... Additional arguments, currently unused.
#'
#' @return A `ggplot2` object.
#' @export
plot.mechte_simest <- function(
  x,
  sd_v = NULL,
  confid = x$alpha_level,
  legend = TRUE,
  xlab = "ATE on the Mediator",
  ylab = "ATE on the Outcome",
  col_error_bar = "#219ebc",
  col_point = "#219ebc",
  col_line = "#8ac926",
  size_title = 14,
  size_text = 12,
  size_legend = 12,
  ...
) {
  dat <- x$plot_data

  if (is.null(dat)) {
    stop("Plot data are not stored in this mechte_simest object.")
  }

  if (!is.null(sd_v)) {
    if (length(sd_v) != nrow(dat)) {
      stop("The length of sd_v should be equal to the number of rows in the fitted object.")
    }
    dat$sd_v <- sd_v
  }

  dat_est <- data.frame(
    alpha = unname(x$alpha),
    beta = unname(x$beta),
    Estimators = "SIMEX"
  )

  p <- ggplot2::ggplot(data = dat, ggplot2::aes(x = gamma_hat, y = tau_hat))

  if (!is.null(sd_v)) {
    p <- p +
      ggplot2::geom_errorbar(
        ggplot2::aes(
          ymin = tau_hat - stats::qnorm(confid / 2, lower.tail = FALSE) * sd_v,
          ymax = tau_hat + stats::qnorm(confid / 2, lower.tail = FALSE) * sd_v
        ),
        colour = col_error_bar,
        width = 0,
        alpha = 0.5
      )
  }

  p +
    ggplot2::geom_errorbarh(
      ggplot2::aes(
        xmin = gamma_hat - stats::qnorm(confid / 2, lower.tail = FALSE) * sd_u,
        xmax = gamma_hat + stats::qnorm(confid / 2, lower.tail = FALSE) * sd_u
      ),
      colour = col_error_bar,
      height = 0,
      alpha = 0.5
    ) +
    ggplot2::geom_point(alpha = 1, colour = col_point) +
    ggplot2::geom_abline(
      data = dat_est,
      ggplot2::aes(intercept = alpha, slope = beta, colour = Estimators),
      linewidth = 1,
      show.legend = legend
    ) +
    ggplot2::scale_color_manual(values = col_line) +
    ggplot2::labs(x = xlab, y = ylab) +
    ggplot2::theme_light() +
    ggplot2::theme(
      axis.text = ggplot2::element_text(size = size_text),
      axis.title = ggplot2::element_text(size = size_title, face = "bold"),
      legend.text = ggplot2::element_text(size = size_legend)
    )
}
