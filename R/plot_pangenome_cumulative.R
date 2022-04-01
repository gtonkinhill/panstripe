#' plot_pangenome_cumulative
#'
#' @importFrom rlang .data
#'
#' @description Plots the cumulative branch length versus the cumulative number of gene gain and loss events. 
#' This is similar to the style of plot produced by TempEst. This function plots only the data points as the fit 
#' can not be included as Panstripe models the individual branches rather than the cumulative total.
#'
#' @param fit the result of running the `panstripe` function. Multiple fits can be passed as a named list.
#' @param plot whether to generate the plot (default) or return a data.frame
#' @param legend toggles the display of the legend on and off
#' @param text_size the base text size of the plot (default=14)
#' @param color_pallete the pallete number passed to `scale_fill_brewer`
#' @param facet whether or not to generate separate plots for each pangenome (default=FALSE).
#'
#' @return either a ggplot2 object or a `data.frame` with the data needed to recreate the plot
#'
#' @examples
#'
#' sim <- simulate_pan(rate=0)
#' fA <- panstripe(sim$pa, sim$tree, nboot=0)
#' plot_pangenome_cumulative(fA, color_pallete=6)
#' sim <- simulate_pan(rate=1e-4)
#' fB <-panstripe(sim$pa, sim$tree, nboot=0)
#' plot_pangenome_cumulative(list(a=fA,b=fB), color_pallete=6)
#' 
#' @export
plot_pangenome_cumulative <- function(fit,
                                plot=TRUE, 
                                legend=TRUE,
                                text_size=14,
                                color_pallete=6,
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
  
  plot_data <- purrr::imap_dfr(fit, ~{
    
    temp_tree <- .x$tree
    temp_tree$edge.length <- .x$data$acc
    
    return(tibble::tibble(
      pangenome=.y,
      acc = ape::node.depth.edgelength(temp_tree),
      core = ape::node.depth.edgelength(.x$tree),
      tip = c(1:(.x$tree$Nnode + length(.x$tree$tip.label)))<length(.x$tree$tip.label)
    ))
  })
  

  if (!plot){
    return(plot_data)
  }
  
  if (length(fit)>1) {
    gg <- ggplot2::ggplot(plot_data, ggplot2::aes(x=.data$core, y=.data$acc, shape=tip, colour=.data$pangenome)) + 
      ggplot2::geom_point()
    if (facet){
      gg <- gg + ggplot2::facet_wrap(~pangenome, ncol = 1)
    }
  } else {
    gg <- ggplot2::ggplot(plot_data, ggplot2::aes(x=.data$core, y=.data$acc, shape=tip))+ 
      ggplot2::geom_point()
  }
  
  gg <- gg + 
    ggplot2::scale_colour_brewer(type = 'qual', palette = color_pallete) +
    ggplot2::scale_fill_brewer(type = 'qual', palette = color_pallete) +
    ggplot2::theme_bw(base_size = text_size) +
    ggplot2::xlab('core phylogentic branch distance') + 
    ggplot2::ylab('accessory distance')
  
  gg
  
  return(gg)
  
}
