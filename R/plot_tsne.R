#' plot_tsne
#'
#' Generates a tSNE plot from a gene presence/absence matrix
#' 
#' @description 
#'
#' @param pa binary presence/absence matrix
#' @param plot whether to return the plot or just the data.frame used for plotting (default=FALSE)
#'
#' @return result
#'
#' @examples
#'
#' sim <- simulate_pan(rate=1e-3)
#' plot_tsne(sim$pa)
#'
#' @export
plot_tsne <- function(pa, category=NULL, plot=TRUE, perplexity=20, pcadims=50){
  
  pcadims <- min(pcadims, nrow(pa))
  perplexity <- min(perplexity, nrow(pa))
  
  pca <- stats::prcomp(pa)
  
  result <- Rtsne::Rtsne(X = pca$x[,1:pcadims], 
                  check_duplicates=FALSE,
                  pca = FALSE,
                  perplexity = 10)
  
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
    gg <- ggplot2::ggplot(plotdf, aes(x=dim1, y=dim2)) +
      ggplot2::geom_point() +
      ggplot2::theme_bw(base_size = 14)
  } else {
    gg <- ggplot2::ggplot(plotdf, aes(x=dim1, y=dim2, colour=category)) +
      ggplot2::geom_point() +
      ggplot2::theme_bw(base_size = 14)
  }
  
  return(gg)
  
}
