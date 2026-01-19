#' get_gsc_weights
#'
#' @param tree A phylogenetic tree object (class 'phylo' from package 'ape').
#' @param normalise Logical. If TRUE, weights are scaled to sum to the number of sequences (standard practice).
#' 
#' @return A named vector of weights for each sequence (tip) in the tree.
get_gsc_weights <- function(tree, normalise = TRUE) {
  
  # Validation and Setup
  if (!inherits(tree, "phylo")) stop("Input must be an object of class 'phylo'")
  if (!ape::is.rooted(tree)) {
    warning("Tree is unrooted. Rooting at midpoint for calculation.")
    tree <- phangorn::midpoint(tree)
  }
  
  # Calculate weights following Gerstein et al., 1994
  temp_tree <- tree
  temp_tree$edge.length <- purrr::map_dbl(1:nrow(tree$edge), ~{
    # w <- sum(ape::extract.clade(tree, tree$edge[.x,1])$edge.length)
    ntips <- length(phangorn::Descendants(tree, tree$edge[.x,2], type='tips')[[1]])
    return(tree$edge.length[.x]/ntips)
  })
  
  # Extract just the distances from the root to everything else
  weights <- ape::node.depth.edgelength(temp_tree)[1:ape::Ntip(temp_tree)]
  
  # Optional: Normalise weights
  # Standard normalization: Sum of weights = Number of Sequences
  # This keeps the "average" weight at 1.0
  if (normalise) {
    weights <- weights/sum(weights) * length(tree$tip.label)
  }
  
  return(weights)
}
