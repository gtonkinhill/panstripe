#' plot_residuals
#'
#' @importFrom rlang .data
#'
#' @description Plots the randomized quantile residuals of the fitted regression model as described in Dunn and Smyth (1996). This plot can be used to assess the fit of the regression.
#'
#' @param fit the result of running the `panstripe` function.
#' @param plot whether to generate the plot (default) or return a data.frame
#' @param text_size the base text size of the plot (default=14)
#' 
#'
#' @return either a ggplot2 object or a `data.frame` with the data needed to recreate the plot
#'
#' @examples
#'
#' sim <- simulate_pan(rate=1e-4)
#' fA <- panstripe(sim$pa, sim$tree, nboot=0)
#' plot_residuals(fA)
#' 
#' @export
plot_residuals <- function(fit,
                           plot=TRUE, 
                           text_size=14){
  # check inputs
  if (class(fit) != 'panfit') {
    stop('fit is not of class `panfit`!')
  } else {
    validate_panfit(fit)
  }
  
  plot_data <- fit$model$data %>%
    tibble::add_column(residuals = statmod::qresiduals(fit$model)) %>%
    tibble::add_column(predicted = predict(fit$model, type = 'response', data=fit$data))
  
  if (!plot){
    return(plot_data)
  }
  
  # bound <- max(c(plot_data$acc, plot_data$predicted))
  # gg1 <- ggplot2::ggplot(plot_data, ggplot2::aes(x=acc, y=predicted)) +
  #   ggplot2::geom_point() +
  #   ggplot2::scale_y_continuous(limits = c(0, bound)) +
  #   ggplot2::scale_x_continuous(limits = c(0, bound)) +
  #   ggplot2::geom_abline(slope=1, intercept = 0, col='red')
  # 
  # gg1 <- gg1 + 
  #   ggplot2::theme_bw(base_size = text_size) +
  #   ggplot2::xlab('accessory distance') +
  #   ggplot2::ylab('predicted')
  # 
  # gg1
  
  gg <- ggplot2::ggplot(plot_data, ggplot2::aes(x=core, y=residuals)) +
    ggplot2::geom_point() +
    ggplot2::geom_hline(yintercept = 0)
  
  gg <- gg + 
    ggplot2::theme_bw(base_size = text_size) +
    ggplot2::xlab('core distance') +
    ggplot2::ylab('Randomized quantile residuals')
  
  gg
  
  return(gg)
  
}
