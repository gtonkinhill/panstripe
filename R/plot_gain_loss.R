#' plot_gain_loss
#' 
#' @importFrom rlang .data
#'
#' @description Plots a expected gain and loss events onto the branches of the phylogeny.
#'
#' @param fit a fitted pangenome model output by running 'panstripe'
#' @param tip_label whether or not to show tree tip labels (default=TRUE)
#' @param text_size adjusts the size of text in the plot
#' @param color_pallete the colour pallete to use. A number between 1 & 9. See 'scale_colour_brewer' for more details
#'
#' @return a plot of the phylogeny coloured by the inferred total gene gain/loss events per branch 
#'
#' @examples
#'
#' sim <- simulate_pan(rate=1e-3)
#' fA <- panstripe(sim$pa, sim$tree, nboot=0)
#' plot_gain_loss(fA, color_pallete=7)
#'
#' @export
plot_gain_loss <- function(fit,
                           tip_label=TRUE,
                           text_size=14,
                           color_pallete=7){
  
  #check inputs
  if (class(fit) != 'panfit') stop('fit is not of class `panfit`!')
  validate_panfit(fit)
  
  stopifnot(all(fit$data$core==fit$tree$edge.length))
  
  gt <- dplyr::full_join(ggtree::fortify(fit$tree), 
                            data.frame(node = fit$tree$edge[,2],
                                       trait = fit$data$acc), 
                            by = 'node')
  
  gg <- ggtree::ggtree(gt, ggplot2::aes(color=.data$trait), size=1) +
    ggplot2::labs(colour='Total genes\ngained & lost') +
    ggplot2::scale_color_binned(type = 'viridis') +
    ggtree::theme_tree2()
  
  if (tip_label){
    gg <- gg + ggtree::geom_tiplab(align=TRUE, colour='black')
  }
    
  return(gg)
}
