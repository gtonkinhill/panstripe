#' plot_dist_params
#'
#' Plots the fitted pangenome tweedie regression model.
#' 
#' @description 
#'
#' @param fit the result of running the `panstripe` function. Multiple results can be passed as a named list.
#' @param plot (default=FALSE)
#'
#' @return result
#'
#' @examples
#'
#' sim <- simulate_pan(rate=1e-3)
#' fA <- panstripe(sim$pa, sim$tree, nboot=10)
#' plot_dist_params(fA)
#' sim <- simulate_pan(rate=1e-2)
#' fB <- panstripe(sim$pa, sim$tree, nboot=10)
#' plot_dist_params(list(a=fA,b=fB))
#'
#' @export
plot_dist_params <- function(fit, plot=TRUE){
  
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
            dplyr::filter(converged) %>%
            dplyr::group_by(core) %>%
            dplyr::summarise(
              "value"=quantile(tpoisson.lambda, 0.5),
              "lower"=quantile(tpoisson.lambda, 0.025),
              "upper"=quantile(tpoisson.lambda, 0.975)
            ) %>%
            tibble::add_column(parameter='Poisson mean'),
          .x$bootrap_replicates %>% 
            dplyr::filter(converged) %>%
            dplyr::group_by(core) %>%
            dplyr::summarise(
              "value"=quantile(tgamma.mean, 0.5),
              "lower"=quantile(tgamma.mean, 0.025),
              "upper"=quantile(tgamma.mean, 0.975)
            ) %>%
            tibble::add_column(parameter='Gamma mean')) %>%
      tibble::add_column(pangenome=.y)
  })
  
  if (!plot){
    return(dist_dat)
  }
  
  if (length(fit)>1) {
    gg <- ggplot2::ggplot(dist_dat, ggplot2::aes(x=core, y=value, col=pangenome)) +
      ggplot2::geom_line() +
      ggplot2::geom_ribbon(ggplot2::aes(ymin=lower, ymax=upper, fill=pangenome), alpha=0.3) +
      ggplot2::facet_wrap(~parameter, ncol = 2, scales = 'free_y')
  } else {
    gg <- ggplot2::ggplot(dist_dat, ggplot2::aes(x=core, y=value)) +
      ggplot2::geom_line() +
      ggplot2::geom_ribbon(ggplot2::aes(ymin=lower, ymax=upper), alpha=0.3) +
      ggplot2::facet_wrap(~parameter, ncol = 2, scales = 'free_y')
  }
  
  gg <- gg + 
    ggplot2::xlab('core phylogentic branch distance') +
    ggplot2::ylab('estimated parameter value') +
    ggplot2::theme_bw(base_size = 14)
  
  return(gg)
  
}
