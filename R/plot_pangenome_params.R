#' plot_pangenome_params
#'
#' @importFrom rlang .data
#'
#' @description Plots the estimated coefficients of the pangenome regression model(s). 
#' It is recommended that the \link{compare_pangenomes} function is used to formally compare the slopes of different pangenome datasets.
#'
#' @param fit the result of running the `panstripe` function. Multiple fits can be passed as a named list.
#' @param boot_ci whether to use the estimated Bootstrap confidence intervals for the 'core' and 'tip' parameters (default=TRUE)
#' @param plot whether to generate the plot (default) or return a data.frame
#' @param legend toggles the display of the legend on and off
#' @param text_size the base text size of the plot (default=14)
#' @param color_pallete the pallete number passed to `scale_fill_brewer`
#' 
#'
#' @return either a ggplot2 object or a `data.frame` with the data needed to recreate the plot
#'
#' @examples
#'
#' sim <- simulate_pan(rate=1e-3)
#' fA <- panstripe(sim$pa, sim$tree, ci_type='perc')
#' plot_pangenome_params(fA, color_pallete=6)
#' sim <- simulate_pan(rate=1e-2)
#' fB <-panstripe(sim$pa, sim$tree, ci_type='perc')
#' plot_pangenome_params(list(a=fA,b=fB), color_pallete=6, boot_ci=TRUE)
#' 
#' @export
plot_pangenome_params <- function(fit,
                                boot_ci=TRUE,
                                plot=TRUE, 
                                legend=TRUE,
                                text_size=14,
                                color_pallete=6){
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
  
  # Transform coefficients of some models prior to plotting
  fit <- purrr::map(fit, ~{
    if (any(grepl('glmmTMB', .x$model$call))){
      family <- .x$model$modelInfo$family$family
    } else {
      family <- .x$model$family$family
    }
    
    if (family %in% c("tweedie", "Tweedie", "quasipoisson", "poisson")){
      k <- .x$summary$term %in% c('core', 'istip', 'depth')
      .x$summary$estimate[k] <- exp(.x$summary$estimate[k])
      .x$summary$`bootstrap CI 2.5%`[k] <- exp(.x$summary$`bootstrap CI 2.5%`[k])
      .x$summary$`bootstrap CI 97.5%`[k] <- exp(.x$summary$`bootstrap CI 97.5%`[k])
    }
    return(.x)
  })

  plot_data <- purrr::imap_dfr(fit, ~{
    .x$summary %>%
      tibble::add_column(pangenome=.y, .before=1)
  }) %>%
    dplyr::filter(.data$term %in% c('core', 'istip', 'phi'))
  plot_data$term[plot_data$term=='istip'] <- 'tip'
  plot_data$term <- factor(plot_data$term, levels = c('core', 'tip', 'phi'))
  
  if (boot_ci){
    plot_data$lower <- plot_data$`bootstrap CI 2.5%`
    plot_data$upper <- plot_data$`bootstrap CI 97.5%`
  } else {
    isphi <- plot_data$term == 'phi'
    plot_data$lower <- plot_data$estimate-plot_data$std.error
    plot_data$upper <- plot_data$estimate+plot_data$std.error
    plot_data$lower[isphi] <- plot_data$`bootstrap CI 2.5%`[isphi]
    plot_data$upper[isphi] <- plot_data$`bootstrap CI 97.5%`[isphi]
  }
  
  if (!plot){
    return(plot_data)
  }
  
  if (length(fit)>1) {
    gg <- ggplot2::ggplot(plot_data, ggplot2::aes(x=.data$pangenome, y=.data$estimate, colour=.data$pangenome))
  } else {
    gg <- ggplot2::ggplot(plot_data, ggplot2::aes(x=.data$term, y=.data$estimate))
  }
    
  gg <- gg + ggplot2::geom_point()
  gg <- gg + ggplot2::geom_errorbar(ggplot2::aes(ymin=lower, ymax=upper))

  if (length(fit)>1) {
    gg <- gg + ggplot2::facet_wrap(~ .data$term, nrow = 1, scales = 'free_y')
    gg <- gg + ggplot2::xlab('')
  } else {
    gg <- gg + ggplot2::xlab('parameter')
  }
  
  gg <- gg + 
    ggplot2::scale_colour_brewer(type = 'qual', palette = color_pallete) +
    ggplot2::scale_fill_brewer(type = 'qual', palette = color_pallete) +
    ggplot2::theme_bw(base_size = text_size) +
    ggplot2::ylab('estimate')
  
  if (!legend){
    gg <- gg + ggplot2::theme(legend.position = 'none')
  }
  
  gg
  
  return(gg)
  
}
