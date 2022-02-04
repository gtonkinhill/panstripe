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
  
  
  dist_dat <- purrr::imap_dfr(fit, function(f, y) {

    sims <- sim_glm(f$model, boot_samples=f$ci_samples)
    sims_ci <- apply(sims, 2, function(x) quantile(x, c(0.025,0.5,0.975)))
    
    p_ci <- quantile(f$ci_samples$t[,6], c(0.025,0.5,0.975))
    phi_ci <- quantile(f$ci_samples$t[,7], c(0.025,0.5,0.975))
    
    sim_params <- purrr::map_dfr(1:nrow(sims), function(i){
      conv <- convert_tweedie(xi = f$ci_samples$t[i,6], 
                      phi = f$ci_samples$t[i,7], 
                      mu = sims[i,])
      tibble::tibble(
        index=1:ncol(sims),
        core=f$model$data$core,
        poisson.lambda=conv$poisson.lambda,
        gamma.mean=conv$gamma.mean,
        gamma.phi=conv$gamma.phi
      )
    }) %>%
      dplyr::group_by(index,core) %>%
      dplyr::summarise(
        parameter=c('Poisson mean','Gamma mean','Gamma dispersion'),
        estimate=c(quantile(poisson.lambda, 0.5),quantile(gamma.mean, 0.5),quantile(gamma.phi, 0.5)),
        lower=c(quantile(poisson.lambda, 0.025),quantile(gamma.mean, 0.025),quantile(gamma.phi, 0.025)),
        upper=c(quantile(poisson.lambda, 0.975),quantile(gamma.mean, 0.975),quantile(gamma.phi, 0.975))
      ) %>% 
      tibble::add_column(pangenome=y, .before=1)
    
  return(sim_params)
    
  }) %>% dplyr::arrange(parameter, core)
  
  dist_dat <- purrr::map_dfr(split(dist_dat, paste(dist_dat$pangenome, dist_dat$parameter)), ~{
    .x <- .x[!duplicated(cut(.x$core, breaks = 1000)),]
    pacc <- stats::smooth.spline(x = .x$core, y = .x$estimate, cv=TRUE, all.knots=TRUE)
    
    tibble::tibble(
      pangenome=unique(.x$pangenome),
      parameter=unique(.x$parameter),
      core=pacc$x,
      estimate=pmax(0, pacc$y),
      lower=pmax(0, stats::smooth.spline(x = .x$core, y = .x$lower, all.knots=TRUE, cv = TRUE)$y),
      upper=pmax(0, stats::smooth.spline(x = .x$core, y = .x$upper, all.knots=TRUE, cv = TRUE)$y)
    )
  })

  if (!plot){
    return(dist_dat)
  }
  
  if (length(fit)>1) {
    gg <- ggplot2::ggplot(dist_dat, ggplot2::aes(x=.data$core, y=.data$estimate, col=.data$pangenome)) +
      ggplot2::geom_line() +
      ggplot2::geom_ribbon(ggplot2::aes(ymin=.data$lower, ymax=.data$upper, fill=.data$pangenome), alpha=0.3) +
      ggplot2::facet_wrap(~parameter, nrow = 1, scales = 'free_y')
  } else {
    gg <- ggplot2::ggplot(dist_dat, ggplot2::aes(x=.data$core, y=.data$estimate)) +
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
