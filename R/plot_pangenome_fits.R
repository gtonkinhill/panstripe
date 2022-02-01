#' plot_pangenome_fits
#'
#' @importFrom rlang .data
#'
#' @description Plots the fitted pangenome tweedie regression model.
#'
#' @param fit the result of running the `panstripe` function. Multiple fits can be passed as a named list.
#' @param ci whether or not to include the confidence interval ribbon (default=TRUE)
#' @param boot whether to use estimated SE for confidence intervals (default) or bootstrap estimates.
#' @param plot whether to generate the plot (default) or return a data.frame
#' @param legend toggles the display of the legend on and off
#' @param text_size the base text size of the plot (default=14)
#' @param color_pallete the pallete number passed to `scale_fill_brewer`
#'
#' @return either a ggplot2 object or a `data.frame` with the data needed to recreate the plot
#'
#' @examples
#'
#' sim <- simulate_pan(rate=0)
#' fA <- panstripe(sim$pa, sim$tree, nboot=10, ci_type='perc')
#' plot_pangenome_fits(fA, color_pallete=6)
#' sim <- simulate_pan(rate=1e-3)
#' fB <- panstripe(sim$pa, sim$tree, nboot=10, ci_type='perc')
#' plot_pangenome_fits(list(a=fA,b=fB), color_pallete=6, boot=FALSE, ci=FALSE)
#'
#' @export
plot_pangenome_fits <- function(fit,
                                ci=TRUE,
                                boot=FALSE,
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
  
  fit_dat <- purrr::imap_dfr(fit, ~{
    p <- predict(.x$model,
                 data.frame(core=rep(max(.x$data$depth)/100, 101), 
                            depth = seq(0, max(.x$data$depth), max(.x$data$depth)/100), 
                            istip=FALSE), 
                 se.fit = TRUE)
    ilink <- family(.x$model)$linkinv
    
    temp_tree <- .x$tree
    temp_tree$edge.length <- .x$data$acc
    
    df <- tibble(
      core = ape::node.depth.edgelength(tree),
      acc = ape::node.depth.edgelength(temp_tree)
    )
    
    
    tibble::tibble(
      pangenome=.y,
      core=seq(0, max(.x$data$depth), max(.x$data$depth)/100),
      val=ilink(p$fit),
      lower=ilink(p$fit - (2 * p$se.fit)),
      upper=ilink(p$fit + (2 * p$se.fit))
    )
  })
  
  
  point_dat <- purrr::imap_dfr(fit, ~{
    tibble::add_column(.x$data, pangenome=.y, .before=1)
  })
  
  if (!plot){
    return(list(
      fit_dat=fit_dat,
      point_dat=point_dat
    ))
  }
  
  if (length(fit)>1) {
    gg <- ggplot2::ggplot(fit_dat, ggplot2::aes(x=.data$core, y=.data$val, colour=.data$pangenome)) +
      ggplot2::geom_line() +
      ggplot2::geom_point(data = point_dat, ggplot2::aes(y=.data$acc, colour=.data$pangenome))
    if (ci){
      gg <- gg + ggplot2::geom_ribbon(ggplot2::aes(ymin=.data$lower, ymax=.data$upper, fill=.data$pangenome), alpha=0.3)
    }
  } else {
    gg <- ggplot2::ggplot(fit_dat, ggplot2::aes(x=.data$core, y=.data$val)) +
      ggplot2::geom_line() +
      ggplot2::geom_point(data = point_dat, ggplot2::aes(y=.data$acc))
    if (ci){
      gg <- gg + ggplot2::geom_ribbon(ggplot2::aes(ymin=.data$lower, ymax=.data$upper), alpha=0.3)
    }
  }
  
  gg <- gg + 
    ggplot2::scale_colour_brewer(type = 'qual', palette = color_pallete) +
    ggplot2::scale_fill_brewer(type = 'qual', palette = color_pallete) +
    ggplot2::theme_bw(base_size = text_size) +
    ggplot2::xlab('core phylogentic branch distance') +
    ggplot2::ylab('estimated accessory distance on branch')

  return(gg)
  
}
