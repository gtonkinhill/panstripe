#' plot_dist_params
#' 
#' @importFrom rlang .data
#'
#' @description Plots the fitted pangenome tweedie regression model.
#'
#' @param fit the result of running the `panstripe` function. Multiple results can be passed as a named list.
#' @param plot (default=FALSE)
#' @param text_size the base text size of the plot (default=14)
#' @param color_pallete the pallete number passed to `scale_fill_brewer`
#'
#' @return either a ggplot2 object or a `data.frame` with the data needed to recreate the plot
#'
#' @examples
#'
#' sim <- simulate_pan(rate=1e-3)
#' fA <- panstripe(sim$pa, sim$tree, nboot=10)
#' plot_dist_params(fA)
#' plot_dist_params(list(a=fA,b=fA))
#'
#' @export
plot_dist_params <- function(fit, plot=TRUE, text_size=14, color_pallete=6){
  
  # check inputs
  if (class(fit)!='panfit'){
    purrr::map(fit, ~{
      if (class(.x) != 'panfit') stop('fit is not of class `panfit`!')
      validate_panfit(.x)
    })
  } else {
    validate_panfit(fit)
    fit <- list(fit)
  }
  
  dist_dat <- purrr::imap_dfr(fit, ~{
    rbind(.x$bootrap_replicates %>% 
            dplyr::group_by(core) %>%
            dplyr::summarise(
              "value"=quantile(tpoisson.lambda, 0.5),
              "lower"=quantile(tpoisson.lambda, 0.025),
              "upper"=quantile(tpoisson.lambda, 0.975)
            ) %>%
            tibble::add_column(parameter='Rate (Poisson mean)'),
          .x$bootrap_replicates %>% 
            dplyr::group_by(core) %>%
            dplyr::summarise(
              "value"=quantile(tgamma.mean, 0.5),
              "lower"=quantile(tgamma.mean, 0.025),
              "upper"=quantile(tgamma.mean, 0.975)
            ) %>%
            tibble::add_column(parameter='Size (Gamma mean)')) %>%
      tibble::add_column(pangenome=.y)
  })
  
  if (!plot){
    return(dist_dat)
  }
  
  if (length(fit)>1) {
    gg <- ggplot2::ggplot(dist_dat, ggplot2::aes(x=.data$core, y=.data$value, col=.data$pangenome)) +
      ggplot2::geom_line() +
      ggplot2::geom_ribbon(ggplot2::aes(ymin=.data$lower, ymax=.data$upper, fill=.data$pangenome), alpha=0.3) +
      ggplot2::facet_wrap(~parameter, ncol = 2, scales = 'free_y')
  } else {
    gg <- ggplot2::ggplot(dist_dat, ggplot2::aes(x=.data$core, y=.data$value)) +
      ggplot2::geom_line() +
      ggplot2::geom_ribbon(ggplot2::aes(ymin=.data$lower, ymax=.data$upper), alpha=0.3) +
      ggplot2::facet_wrap(~parameter, ncol = 2, scales = 'free_y')
  }
  
  gg <- gg + 
    ggplot2::xlab('core phylogentic branch distance') +
    ggplot2::ylab('estimated parameter value') +
    ggplot2::theme_bw(base_size = text_size) +
    ggplot2::scale_colour_brewer(type = 'qual', palette = color_pallete) +
    ggplot2::scale_fill_brewer(type = 'qual', palette = color_pallete)
  
  return(gg)
  
}
