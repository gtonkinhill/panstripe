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
#' @param boot_ci whether to use bootstrap CI or standard error
#'
#' @return either a ggplot2 object or a `data.frame` with the data needed to recreate the plot
#'
#' @examples
#'
#' sim <- simulate_pan(rate=1e-3, mean_trans_size = 1)
#' fit <- panstripe(sim$pa, sim$tree, nboot=100, ci_type='perc')
#' plot_dist_params(fit)
#' simB <- simulate_pan(rate=1e-3, mean_trans_size = 10)
#' fitB <- panstripe(simB$pa, simB$tree, nboot=100, ci_type='perc')
#' plot_dist_params(list(a=fit,b=fitB))
#'
#' @export
plot_dist_params <- function(fit, plot=TRUE, text_size=14, color_pallete=6, boot_ci=FALSE){
  
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
  
  dist_dat <- purrr::imap_dfr(fit, function(f, y) {
    
    coef_names <- c('Intercept', strsplit(as.character(f$model$formula), ' \\+ ')[[3]])
    if (ncol(f$ci_samples$t)>length(coef_names)){
      coef_names <- c(coef_names, 'p', 'phi')
    }
    
    if (boot_ci){
      df <- tibble::tibble(
        param=f$summary$term,
        estimate = f$summary$estimate,
        lower = f$summary$`bootstrap CI 2.5%`,
        upper = f$summary$`bootstrap CI 97.5%`
      )
    } else {
      df <- tibble::tibble(
        param=f$summary$term,
        estimate = f$summary$estimate,
        lower = f$summary$estimate-2*f$summary$std.error,
        upper = f$summary$estimate+2*f$summary$std.error
      )
    }
    
    df$pangenome=y
    df <- df[df$param!='Intercept',]
    df$param <- factor(df$param, levels = coef_names[-1])
    return(df)
  })
  
  
  if (!plot){
    return(dist_dat)
  }
  
  if (length(fit)>1) {
    gg <- ggplot2::ggplot(dist_dat, ggplot2::aes(x=.data$pangenome, y=.data$estimate, col=.data$pangenome, group=.data$pangenome)) +
      ggplot2::geom_point(size=2) +
      ggplot2::geom_errorbar(ggplot2::aes(ymin=.data$lower, ymax=.data$upper), width=.2) +
      ggplot2::facet_wrap(~param, nrow = 1, scales = 'free_y')
  } else {
    gg <- ggplot2::ggplot(dist_dat, ggplot2::aes(x=.data$param, y=.data$estimate, group=.data$param)) +
      ggplot2::geom_point(size=2) +
      ggplot2::geom_errorbar(ggplot2::aes(ymin=.data$lower, ymax=.data$upper), width=.2)
  }
  
  gg <- gg + 
    ggplot2::xlab('') +
    ggplot2::ylab('estimated parameter value') +
    ggplot2::theme_bw(base_size = text_size) +
    ggplot2::scale_colour_brewer(type = 'qual', palette = color_pallete) +
    ggplot2::scale_fill_brewer(type = 'qual', palette = color_pallete)
  # gg
  return(gg)
  
}

convert_tweedie <- function(xi, mu, phi){
  return(list(
    poisson.lambda=(mu^(2-xi))/(phi*(2-xi)),
    gamma.mean=(2-xi)*phi*(mu^(xi-1)),
    gamma.phi=(2-xi)*(xi-1)*(phi^2)*(mu^(2*(xi-1)))
  ))
}
