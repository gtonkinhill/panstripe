#' plot_umap
#'
#' @importFrom rlang .data
#'
#' @description Generates a UMAP plot from a gene presence/absence matrix
#'
#' @param pa binary presence/absence matrix
#' @param category a factor vector which can be used to colour the points
#' @param plot whether to return the plot or just the data.frame used for plotting (default=TRUE)
#' @param n_neighbours number of neighbours with which the UMAP embedding is constructed
#' @param dim dimensionality of the target output space
#' @param n_epochs number of iterations performed during layout optimisation
#' @param metric distance metric to determine distance between data points (see umap.defaults documentation for options)
#'
#' @return either a ggplot2 object or a `data.frame` with the data needed to recreate the plot
#'
#' @examples
#'
#' sim <- simulate_pan(rate=1e-3)
#' plot_umap(sim$pa)
#'
#' @export
plot_umap <- function(pa, category=NULL, plot=TRUE, n_neighbours=15, dim=2, n_epochs=200, metric="euclidean"){
  pa.umap <- umap::umap(pa,
                       n_neighbours = n_neighbours,
                       n_components = dim,
                       metric = metric,
                       n_epochs = n_epochs)
  
  plotdf = data.frame(pa.umap$layout)
  
  if (!plot){
    return(plotdf)
  }
  
  if (is.null(category)){
    gg <- ggplot2::ggplot(plotdf, ggplot2::aes(x=.data$X1, y=.data$X2)) +
      ggplot2::geom_point() +
      ggplot2::theme_bw(base_size = 14) +
      ggplot2::labs(x = 'dim1', y = 'dim2')
  } else {
    gg <- ggplot2::ggplot(plotdf, ggplot2::aes(x=.data$X1, y=.data$X2, colour=category)) +
      ggplot2::geom_point() +
      ggplot2::theme_bw(base_size = 14) +
      ggplot2::labs(x = 'dim1', y = 'dim2')
  }
  
  return(gg)
  
}