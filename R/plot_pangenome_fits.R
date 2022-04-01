#' plot_pangenome_fits
#'
#' @importFrom rlang .data
#'
#' @description Plots the fitted pangenome regression model. 
#' Optionally the data points representing the branches at the tips of the phylogeny can be included. 
#' Internal branches are not plotted as this would require a third dimension to account for the depth of the branch. 
#' The fit of the model can be further assessed using the \link{plot_residuals} function.
#'
#' @param fit the result of running the `panstripe` function. Multiple fits can be passed as a named list.
#' @param ci whether or not to include the confidence interval ribbon (default=TRUE)
#' @param plot whether to generate the plot (default) or return a data.frame
#' @param legend toggles the display of the legend on and off
#' @param text_size the base text size of the plot (default=14)
#' @param color_pallete the pallete number passed to `scale_fill_brewer`
#' @param include_data whether or not to include raw data points in plot (default=FALSE) Note: The ribbon indicates the uncertainty in the estimated mean and thus not all points will lie within these margins.
#' @param trim whether or not to trim the plots to cover the same range for each pangenome (default=TRUE)
#' @param facet whether or not to generate separate plots for each pangenome (default=FALSE).
#' 
#'
#' @return either a ggplot2 object or a `data.frame` with the data needed to recreate the plot
#'
#' @examples
#'
#' sim <- simulate_pan(rate=1e-4)
#' fA <- panstripe(sim$pa, sim$tree, nboot=10, ci_type='perc')
#' plot_pangenome_fits(fA, color_pallete=6, include_data=TRUE)
#' sim <- simulate_pan(rate=1e-3)
#' fB <-panstripe(sim$pa, sim$tree, nboot=10, ci_type='perc')
#' plot_pangenome_fits(list(a=fA,b=fB), color_pallete=6, ci=TRUE)
#' 
#' @export
plot_pangenome_fits <- function(fit,
                                ci=TRUE,
                                plot=TRUE, 
                                legend=TRUE,
                                text_size=14,
                                color_pallete=6,
                                include_data=FALSE,
                                trim=TRUE,
                                facet=FALSE){
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
  
  plot_data <- purrr::imap(fit, ~{

    point_data <- .x$data %>% 
      tibble::add_column(pangenome=.y, .before = 1) %>% 
      dplyr::filter(istip)
    
    
    ilink <- family(.x$model)$linkinv
    
    p <- predict(.x$model, type = 'link', se.fit = TRUE,
                 newdata=data.frame(
                   core = seq(0,max(point_data$core), length.out=100),
                   istip = TRUE,
                   depth = max(.x$model$data$depth)-seq(0,max(point_data$core), length.out=100)
                 ))

    fit_data <- tibble::tibble(
      pangenome=.y,
      core = seq(0,max(point_data$core), length.out=100),
      acc = ilink(p$fit),
      lower = ilink(p$fit - 2*p$se.fit),
      upper = ilink(p$fit + 2*p$se.fit),
    )
    
    return(list(point_data=point_data, fit_data=fit_data))
  })
  
  if (trim){
    max_core <- min(purrr::map_dbl(plot_data, ~ max(.x$point_data$core)))
  } else {
    max_core <- Inf
  }
  
  plot_data <- list(
    point_data=purrr::map_dfr(plot_data, ~ .x$point_data) %>%
      dplyr::filter(core<max_core),
    fit_data=purrr::map_dfr(plot_data, ~ .x$fit_data)%>%
      dplyr::filter(core<max_core)
  )
  
  if (!plot){
    return(plot_data)
  }
  
  if (length(fit)>1) {
    gg <- ggplot2::ggplot(plot_data$fit_data, ggplot2::aes(x=.data$core, y=.data$acc, colour=.data$pangenome))
    if (include_data) {
      gg <- gg + ggplot2::geom_point(data = plot_data$point_data, ggplot2::aes(y=.data$acc, colour=.data$pangenome), alpha=1)
    }
    gg <- gg + ggplot2::geom_line() 
    if (ci && !any(is.na(plot_data$fit_data$upper))){
      gg <- gg + ggplot2::geom_ribbon(ggplot2::aes(ymin=.data$lower, ymax=.data$upper, fill=.data$pangenome), alpha=0.3)
    }
    if (facet){
      gg <- gg + ggplot2::facet_wrap(~pangenome, ncol = 1)
    }
  } else {
    gg <- ggplot2::ggplot(plot_data$fit_data, ggplot2::aes(x=.data$core, y=.data$acc))
    if (include_data) {
      gg <- gg + ggplot2::geom_point(data = plot_data$point_data, ggplot2::aes(y=.data$acc))
    }
    gg <- gg + ggplot2::geom_line()
    if (ci && !any(is.na(plot_data$fit_data$upper))){
      gg <- gg + ggplot2::geom_ribbon(ggplot2::aes(ymin=.data$lower, ymax=.data$upper), alpha=0.3)
    }
  }
  
  gg <- gg + 
    ggplot2::scale_colour_brewer(type = 'qual', palette = color_pallete) +
    ggplot2::scale_fill_brewer(type = 'qual', palette = color_pallete) +
    ggplot2::theme_bw(base_size = text_size) +
    ggplot2::xlab('core phylogentic branch distance') #+
  
  
  gg <- gg + ggplot2::ylab('accessory distance')
  
  gg
  
  return(gg)
  
}
