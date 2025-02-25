#' plot_tree_pa
#'
#' @importFrom rlang .data
#'
#' @description Plots a phylogeny alongside a presence/absence matrix of user selected genes. 
#'
#' @param tree a phylogeny in 'phylo' format
#' @param pa binary presence/absence matrix
#' @param genes a vector of gene to include in the plot. Must be a subset of colnames(pa) 
#' @param align whether or not to align the tips of the phylogeny with the presence/absence matrix using dotted lines
#' @param order whether to order the genes based upon their presence/absence pattern
#' @param plot_titles the subplot titles (default=c('phylogeny', 'gene presence/absence'))
#' @param text_size adjusts the size of text in the plot
#' @param label_genes whether to print the gene labels in the presence/absence matrix
#' @param label_tips whether to print the tree tip labels
#' @param cols colours for the gene presence/absence matrix. These will be recycled if #genes > #colours.
#'
#' @return a patchwork object plotting the phylogeny alongside the gene presence/absence matrix
#'
#' @examples
#'
#' sim <- simulate_pan(rate=1e-3)
#' genes <- colnames(sim$pa)[which(apply(sim$pa, 2, sd)>0.2)]
#' plot_tree_pa(sim$tree, sim$pa, genes=genes, 
#'              label_genes=FALSE, label_tips=FALSE, cols='black')
#'
#' @export
plot_tree_pa <- function(tree, pa, genes=colnames(pa), 
                         align=TRUE,
                         order=TRUE,
                         plot_titles=c('phylogeny', 'gene presence/absence'),
                         text_size=14,
                         label_genes=TRUE,
                         label_tips=TRUE,
                         cols=NULL){
  
  if (is.null(cols)){
    cols <- c('#a6cee3','#1f78b4','#b2df8a','#33a02c','#fb9a99','#e31a1c',
                 '#fdbf6f','#ff7f00','#cab2d6','#6a3d9a','#ffff99','#b15928')
  }
  cols <- rep(cols, ceiling(length(genes)/length(cols)))
  
  if('ggtree' %in% class(tree)){
    ggt <- tree
    ntips <- sum(tree$data$isTip)
  } else {
    if (!class(tree)=='phylo') stop('tree must be either a phylo or ggtree object!')
    ttree <- tree
    if(!label_tips){
      ttree$tip.label <- rep(NA, length(ttree$tip.label)) 
    }
    if (align){
      ggt <- ggtree::ggtree(ttree) + 
        ggtree::geom_tiplab(align=TRUE, col='grey')
    } else {
      ggt <- ggtree::ggtree(ttree)
    }
    ntips <- length(tree$tip.label)
  }
  ggt <- ggt + 
    ggplot2::ggtitle(plot_titles[[1]]) +
    ggplot2::theme(plot.title = ggplot2::element_text(size=text_size)) +
    ggplot2::scale_y_continuous(limits = c(0, ntips+1)) +
    ggplot2::theme(plot.margin = ggplot2::margin(t=0, b=0, r=0, l=0, unit = "pt"),
                   plot.title = ggplot2::element_text(hjust = 0.5))

  subset_pa <- pa[, colnames(pa) %in% genes, drop=FALSE]
  
  if ((ncol(subset_pa)>2) & order){
      d <- stats::as.dist(ncol(subset_pa) - tcrossprod(t(subset_pa)))
      h <- stats::hclust(d, method = 'average')
      subset_pa <- subset_pa[,h$order,drop=FALSE]
  }
  
  padf <- tibble::tibble(isolate = rep(rownames(subset_pa), ncol(subset_pa)),
                        gene = rep(colnames(subset_pa), each=nrow(subset_pa)),
                        presence = c(subset_pa))
  padf$gene <- factor(padf$gene, levels=colnames(subset_pa))
  padf <- padf[padf$presence>0, , drop=FALSE]
  
  # plot links against tree
  if(!'ggtree' %in% class(tree)){
    d <- ggtree::fortify(tree)
  } else {
    d <- ggt$data
  }
  d <- d[d$isTip, , drop=FALSE]

  tip_order <- with(d, label[order(y, decreasing=TRUE)])
  nodes <- tibble::tibble(nodes = tip_order, pos = ntips:1)

  padf$height <- nodes$pos[match(padf$isolate, nodes$nodes)]
  
  gg <- ggplot2::ggplot(padf, ggplot2::aes(x=.data$gene, y=.data$height, fill=.data$gene)) +
    ggplot2::geom_tile(col='white') +
    ggplot2::scale_fill_manual(values=cols) +
    ggplot2::scale_y_continuous(limits = c(0, ntips+1)) +
    ggplot2::theme_bw(base_size = text_size) +
    ggplot2::theme(axis.title.y=ggplot2::element_blank(),
          axis.text.y=ggplot2::element_blank(),
          axis.ticks.y=ggplot2::element_blank(),
          axis.ticks.x=ggplot2::element_blank()) +
    ggplot2::theme(panel.background = ggplot2::element_blank(),
                   panel.grid.major.y = ggplot2::element_blank(),
                   panel.grid.minor.y = ggplot2::element_blank(),
                   panel.border = ggplot2::element_blank(),
                   legend.position = "none") +
    ggplot2::ggtitle(plot_titles[[2]]) +
    ggplot2::theme(plot.title = ggplot2::element_text(size=text_size)) +
    ggplot2::theme(plot.margin = ggplot2::margin(t=0, b=0, r=0, l=0, unit = "pt"),
                   plot.title = ggplot2::element_text(hjust = 0.5))
  
  if (!label_genes){
    gg <- gg +
      ggplot2::theme(axis.text.x=ggplot2::element_blank(),
                     axis.ticks.x=ggplot2::element_blank(),
                     panel.grid.major.x = ggplot2::element_blank(),
                     panel.grid.minor.x = ggplot2::element_blank()) 
  }
  
  pp <- patchwork::wrap_plots(ggt, gg) + 
    patchwork::plot_layout(nrow = 1)
  return(pp)
  
}
