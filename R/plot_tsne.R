#' plot_tsne
#'
#' @importFrom rlang .data
#'
#' @description Generates a tSNE plot from a gene presence/absence matrix
#'
#' @param pa binary presence/absence matrix
#' @param category a factor vector which can be used to colour the points
#' @param plot whether to return the plot or just the data.frame used for plotting (default=TRUE)
#' @param perplexity the `perplexity` argument passed to tSNE
#' @param pcadims the number of principle components passed from the inital stage used in the main tSNE algorithm
#'
#' @return either a ggplot2 object or a `data.frame` with the data needed to recreate the plot
#'
#' @examples
#'
#' sim <- simulate_pan(rate=1e-3)
#' plot_tsne(sim$pa)
#'
#' @export
plot_tsne <- function(pa, category=NULL, plot=TRUE, perplexity=20, pcadims=50){
  
  pcadims <- min(pcadims, nrow(pa))
  if (floor((nrow(pa) - 1) / 3) < perplexity) {
    perplexity <- floor((nrow(pa) - 1) / 3)
    message(paste("Small number of samples, reducing perplexity to", perplexity))
  }
  
  pca <- stats::prcomp(pa)
  
  result <- Rtsne::Rtsne(X = pca$x[,1:pcadims], 
                         check_duplicates=FALSE,
                         pca = FALSE,
                         perplexity = perplexity)
  
  if (is.null(category)){
    stopifnot(all(names(category)==rownames(pa)))
    plotdf <- tibble::tibble(
      dim1 = result$Y[,1],
      dim2 = result$Y[,2]
    )
  } else {
    plotdf <- tibble::tibble(
      category = category,
      dim1 = result$Y[,1],
      dim2 = result$Y[,2]
    )
  }
  
  if (!plot){
    return(plotdf)
  }
  
  if (is.null(category)){
    gg <- ggplot2::ggplot(plotdf, ggplot2::aes(x=.data$dim1, y=.data$dim2)) +
      ggplot2::geom_point() +
      ggplot2::theme_bw(base_size = 14)
  } else {
    gg <- ggplot2::ggplot(plotdf, ggplot2::aes(x=.data$dim1, y=.data$dim2, colour=category)) +
      ggplot2::geom_point() +
      ggplot2::theme_bw(base_size = 14)
  }
  
  return(gg)
  
}
