#' sfreemap
#'
#' Simulation free stochastic character mapping on a phylogenetic tree.
#' 
#' @description This function has been adapted and simplified from the code in the sfreemap R package.
#'
#' @param tree 
#' @param pa 
#' @param type 
#' @param model 
#' @param method 
#' @param prior
#' @param tol
#'
#' @return result
#'
#' @examples
#'
#' tree <- ape::rtree(10)
#' pa <- matrix(rbinom(10, 1, 0.5), ncol=1)
#' rownames(pa) <- tree$tip.label
#' sfreemap(tree, pa)
#' 
#'
#' @export
sfreemap <- function(tree, pa, type="standard", model="ER"){
  
  # Reorder the tree so the root is the first row of the matrix.
  # We save the original order to make sure we have the result
  # in same order of the tree;
  tree <- ape::reorder.phylo(tree, order='pruningwise')
  
  ngenes <- ncol(pa)
  dt <- matrix(0, nrow = nrow(tree$edge), ncol = ngenes, dimnames = list(NULL, colnames(pa)))
  
  cat('Estimating expected gene gain/loss events...\n')
  for ( i in 1:ngenes ){
    tip_states <- pa[,i]
    QP <- Q_empirical(tree, tip_states, model)
    # calculate eigen vectors
    Q_eigen <- eigen(QP$Q, TRUE, only.values = FALSE, symmetric = TRUE)
    Q_eigen$vectors_inv <- solve(Q_eigen$vectors)
    
    # build state matrix
    tip_states_matrix <- matrix(0, nrow = length(tip_states), ncol = 2)
    tip_states_matrix[cbind(1:length(tip_states), tip_states+1)] <- 1
   
    # calculate transition probabilities
    tp <- transition_probabilities(tree$edge.length, Q_eigen)
    
    # calculate fractional likelihoods
    fl <- fractional_likelihoods(tree$edge, tip_states_matrix, QP$prior, tp)
    
    #calculate the dwelling times
    dt[,i] <- rowSums(dwelling_times(tree, Q_eigen, QP$Q, fl))
    
    cat(paste0(round(i / ngenes * 100), '% completed\r'))
    if (i == ngenes) cat('Done\n')
  }
  
  new_tree <- ape::reorder.phylo(tree, order='cladewise')
  new_tree[['mapped.edge.lmt']] <- dt[match(paste(new_tree$edge[,1], new_tree$edge[,2]),
                                            paste(tree$edge[,1], tree$edge[,2])),]
  
  return(new_tree)
}


transition_probabilities <- function(edges, Q_eigen){
  a <- exp(matrix(edges) %*% matrix(Q_eigen$values, nrow=1))
  return(purrr::map(1:nrow(a), ~ Q_eigen$vectors_inv %*% (Q_eigen$vectors * a[.x,])))
}

Q_empirical <- function(tree, tip_states, model='ER'){
  
  if (model=='ER'){
    m <- matrix(c(0, 1, 1, 0), 2)
  } else if (model=='ARD') {
    m <- matrix(c(0, 1, 2, 0), 2)
  }
  
  fit <- ape::ace(tip_states, tree, model = m, type = 'discrete')
  Q <- matrix(fit$rates[fit$index.matrix], 2)
  Q[c(1,4)] <- -Q[c(2,3)]
  
  return(list(
    Q=Q,
    prior=c(0.5, 0.5),
    logL=fit$loglik
  ))
}


