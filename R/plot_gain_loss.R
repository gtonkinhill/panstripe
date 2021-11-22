#' plot_gain_loss
#'
#' Plots a expected gain and loss events onto the branches of the phylogeny.
#' 
#' @description 
#'
#' @param res the output of running 'fit_tweedie'
#' @param adjusted whether of not to adjust the plot for errors (see description)
#' @param text_size adjusts the size of text in the plot
#' @param label_genes whether to print the gene labels in the presence/absence matrix
#' @param color_pallete the colour pallete to use. A number between 1 & 9. See 'scale_colour_brewer' for more details
#'
#' @return a plot of the phylogeny coloured by the expected total gene gain/loss events per branch 
#'
#' @examples
#'
#' sim <- simulate_pan(rate=1e-3)
#' res <- fit_tweedie(sim$pa, sim$tree)
#' plot_gain_loss(res, 'both', color_pallete=7)
#'
#' @export
plot_gain_loss <- function(res,
                           adjusted='both',
                           text_size=14,
                           color_pallete=7){
  
  stopifnot(all(res$points$core==res$anc$tree$edge.length))
  invisible(capture.output({ef <- tweedie::tweedie.profile(acc ~ core, data = res$points, p.vec = seq(1.1,1.9,0.1),
                                                           do.smooth = FALSE, method="series")}))
  m <- glm(acc ~ core, data = res$points, family = statmod::tweedie(var.power = ef$p.max, link.power = 0))
  
  data <- ggtree::fortify(res$anc)
  
  if (adjusted=='no'){
    data$gl <- res$points$acc[match(data$parent, res$anc$edge[,1])]
  } else if (adjusted=='yes'){
    data$gl <- res$points$acc[match(data$parent, res$anc$edge[,1])] - 
      residuals(m, type = "response")[match(data$parent, res$anc$edge[,1])]
  } else {
    gl <- c(res$points$acc[match(data$parent, res$anc$edge[,1])],
            res$points$acc[match(data$parent, res$anc$edge[,1])] - 
              residuals(m, type = "response")[match(data$parent, res$anc$edge[,1])])
    data <- rbind(data, data)
    data$gl <- gl
    data$adjusted <- rep(c('No Adjustment', 'Adjusted'), each=nrow(data)/2)
  }
  
  gg <- ggtree::ggtree(data, ggplot2::aes(col=gl)) +
    ggplot2::scale_color_fermenter(type = 'div', palette = color_pallete) +
    ggplot2::theme(legend.position="right",
                   legend.text=ggplot2::element_text(size=text_size),
                   legend.title=ggplot2::element_text(size=text_size),
                   strip.text.x = ggplot2::element_text(size = text_size)) +
    ggplot2::labs(col='Expected gene\ngain/loss events')
  
  if (adjusted=='both'){
    gg <- gg +
      ggplot2::facet_wrap(~adjusted)
  }
  
  return(gg)
  
}
