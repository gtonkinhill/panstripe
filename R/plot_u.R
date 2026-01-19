#' plot_u
#'
#' @importFrom rlang .data
#'
#' @description Plots a pangenome U-plot of gene frequencies with optional weights for phylogeny.
#'
#' @param pa binary presence/absence matrix  
#' @param tree a phylogeny in 'phylo' format
#' @param plot logical indicating whether to display the plot (default) or just return a `data.frame`
#' @param text_size adjusts the size of text in the plot
#' @param bin_width the width of the histogram bins
#' @param color_pallete a number indicating which RColorBrewer palette to use
#'
#' @return a ggplot2 object
#'
#' @examples
#'
#' sim <- simulate_pan(rate=1e-3)
#' plot_u(sim$pa)
#' plot_u(sim$pa, sim$tree)
#'
#' @export
plot_u <- function(pa, 
                   tree=NULL,
                   plot = TRUE,
                   text_size=14,
                   bin_width = 5,
                   color_pallete=6){
  
  
  if (!is.null(tree)){
    #Calculate weights using method described in Gerstein et al., 1994
    norm_weights <- get_gsc_weights(tree)
  }
  
  # d <- data.frame(label = names(norm_weights), weight = get_gsc_weights(tree))
  # d <- as.data.frame(picante::evol.distinct(tree))
  # colnames(d) <- c('label', 'weight')
  # ggtree(tree) %<+% d +                  # Use %<+% to attach the data frame to the tree
  #   geom_tiplab(offset = 0.2) +          # Add species names
  #   geom_tippoint(aes(size = weight,     # Scale size by weight
  #                     color = weight)) + # Scale color by weight
  #   scale_color_viridis_c() +            # Use a nice color gradient
  #   theme_tree2() +                      # Add a scale bar for branch lengths
  #   labs(title = "Phylogenetic Normalization Weights",
  #        subtitle = "Larger/brighter points indicate higher weight (unique lineages)",
  #        size = "Weight",
  #        color = "Weight")
  
  if (!is.null(tree)){
    freqs <- colSums(pa * norm_weights)
  } else {
    freqs <- colSums(pa)
  }
  
  # plot gene frequencies
  pdf <- data.frame(
    gene = colnames(pa),
    freq = freqs
  )
  
  # filter out zero frequency genes
  pdf <- pdf[pdf$freq > 0, ]
  
  if (!plot){
    return(pdf)
  }
  
  gg <- ggplot2::ggplot(pdf, ggplot2::aes(x=freq)) +
    ggplot2::geom_histogram(binwidth = bin_width) +
    ggplot2::scale_colour_brewer(type = 'qual', palette = color_pallete) +
    ggplot2::scale_fill_brewer(type = 'qual', palette = color_pallete) +
    ggplot2::theme_bw(base_size = text_size) +
    ggplot2::labs(
      x = "Gene Frequency",
      y = "Number of Genes",
      title = "Pangenome U-Plot"
    )
  
  return(gg)
  
}
