#' simulate_pan
#'
#' @description Simulates a pangenome presence/absence matrix and accompanying core phylogeny
#'
#' @param ngenomes number of genomes to simulate
#' @param rate rate of gene gain/loss events
#' @param genome_length number of genes
#' @param fn_error_rate false negative error rate
#' @param fp_error_rate false positive error rate
#' @param mean_trans_size mean size of the number of genes gained or lost per event
#' @param Q the rate matrix (default=1)
#' @param bf the base frequencies (default=c(0.3,0.2))
#'
#' @return A list comprising of a binary gene presence/absence matrix and a core phylogeny 
#'
#' @examples
#'
#' pan <- simulate_pan()
#' 
#'
#' @export
simulate_pan <- function(ngenomes=50, rate=1e-5, genome_length=2000, 
                         fn_error_rate=1, fp_error_rate=3, 
                         mean_trans_size=3, Q=1, bf=c(0.3,0.2)){
  
  # simulate tree
  tree <- ape::rtree(ngenomes)
  
  #simulate genes
  data <- as.character(phangorn::simSeq(tree, l = genome_length, type="USER", 
                              levels = c(0,1), bf=bf, 
                              Q=Q, rate = rate))
  
  
  is_acc <- which(apply(data, 2, stats::sd)>0)
  is_core <- which(apply(data, 2, stats::sd)==0)
  
  #convert non-variable to core
  data[, is_core] <- 1
  
  #model number of genes involved in each event
  acc <- purrr::map2(is_acc, 1+stats::rpois(length(is_acc), lambda = mean_trans_size-1), ~{
    matrix(rep(data[,.x], .y), byrow = FALSE, ncol = .y)})
  acc <- acc[purrr::map_lgl(acc, ~ nrow(.x)>0)]
  data <- cbind(data[, is_core], 
                do.call(cbind, acc))
  
  #add errors
  # false negatives
  for (i in 1:nrow(data)){
    nerrors <- stats::rpois(1, fn_error_rate)
    if (nerrors>0){
      data[i,sample(which(data[i,] == '1'), nerrors, replace = FALSE)] <- '0'  
    }
  }
  
  #false positives
  data <- cbind(data, 
                do.call(cbind, 
                        purrr::imap(stats::rpois(nrow(data), fp_error_rate), 
                             ~ {
                               if(.y>0){
                                 tmp <- matrix('0', ncol = .x, nrow = ngenomes)
                                 tmp[.y,] <- '1'
                                 return(tmp)
                               } else{
                                 return(NULL)
                               }
                             })))
  colnames(data) <- 1:ncol(data)
  class(data) <- "numeric"
  
  return(list(pa=data, tree=tree))
  
}
