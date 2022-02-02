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
#' fit <- panstripe(sim$pa, sim$tree, nboot=10, ci_type='perc')
#' plot_dist_params(fit)
#' plot_dist_params(list(a=fit,b=fit))
#'
#' @export
plot_dist_params <- function(fit, plot=TRUE, text_size=14, color_pallete=6){
  
  # check inputs
  if (class(fit)!='panfit'){
    purrr::map(fit, ~{
      if (class(.x) != 'panfit') stop('fit is not of class `panfit`!')
      panstripe:::validate_panfit(.x)
    })
  } else {
    panstripe:::validate_panfit(fit)
    fit <- list(fit)
  }
  
  
  dist_dat <- purrr::imap_dfr(fit, ~{
    
    p_ci <- boot_ci_pval(.x$ci_samples, 5, 'norm', transformation='logit')
    phi_ci <- boot_ci_pval(.x$ci_samples, 6, 'norm', transformation='inverse')
    
    estimate <- convert_tweedie(xi = .x$model$p, 
                    phi = .x$model$phi, 
                    mu = c(0:max(.x$data$acc)))
    lower <- convert_tweedie(xi = p_ci[[1]], 
                             phi = phi_ci[[1]], 
                             mu = c(0:max(.x$data$acc)))
    upper <- convert_tweedie(xi = p_ci[[2]], 
                             phi = phi_ci[[2]], 
                             mu = c(0:max(.x$data$acc)))
    
    tibble::tibble(
      pangenome=.y,
      parameter=gsub('[0-9]+','',names(unlist(estimate))),
      tweedie_mean=rep(c(0:max(.x$data$acc)), 3),
      estimate=unlist(estimate),
      lower=unlist(lower),
      upper=unlist(upper),
    )
    
  })
  
  new_names <- c(poisson.lambda='Poisson mean',gamma.mean='Gamma mean',gamma.phi='Gamma dispersion')
  dist_dat$parameter <- factor(new_names[dist_dat$parameter], levels=new_names)
  
  if (!plot){
    return(dist_dat)
  }
  
  if (length(fit)>1) {
    gg <- ggplot2::ggplot(dist_dat, ggplot2::aes(x=.data$tweedie_mean, y=.data$estimate, col=.data$pangenome)) +
      ggplot2::geom_line() +
      ggplot2::geom_ribbon(ggplot2::aes(ymin=.data$lower, ymax=.data$upper, fill=.data$pangenome), alpha=0.3) +
      ggplot2::facet_wrap(~parameter, nrow = 1, scales = 'free_y')
  } else {
    gg <- ggplot2::ggplot(dist_dat, ggplot2::aes(x=.data$tweedie_mean, y=.data$estimate)) +
      ggplot2::geom_line() +
      ggplot2::geom_ribbon(ggplot2::aes(ymin=.data$lower, ymax=.data$upper), alpha=0.3) +
      ggplot2::facet_wrap(~parameter, nrow = 1, scales = 'free_y')
  }
  
  gg <- gg + 
    ggplot2::xlab('mean of Tweedie distribution') +
    ggplot2::ylab('estimated parameter value') +
    ggplot2::theme_bw(base_size = text_size) +
    ggplot2::scale_colour_brewer(type = 'qual', palette = color_pallete) +
    ggplot2::scale_fill_brewer(type = 'qual', palette = color_pallete)
  
  return(gg)
  
}

convert_tweedie <- function(xi, mu, phi){
  return(list(
    poisson.lambda=(mu^(2-xi))/(phi*(2-xi)),
    gamma.mean=(2-xi)*phi*(mu^(xi-1)),
    gamma.phi=(2-xi)*(xi-1)*(phi^2)*(mu^(2*(xi-1)))
  ))
}
