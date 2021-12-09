#' read_rtab
#'
#' @description Loads an Rtab pangenome gene presence/absence file as produced by Panaroo, Roary and other prokaryotic pangenome programs
#'
#' @param file a binary tab delimited gene presence/absence matrix (.Rtab) with gene names in the first column
#'
#' @return a binary matrix
#'
#' @examples
#'
#' file <- system.file("extdata", "gene_presence_absence.Rtab", package = "panstripe")
#' pa <- read_rtab(file)
#'
#' @export
read_rtab <- function(file){
  
  # check inputs
  if (!file.exists(file)) stop('File does not exist!')
  
  df <- readr::read_tsv(file, col_names = TRUE)
  pa <- data.matrix(df[, 2:ncol(df), drop=FALSE])
  rownames(pa) <- unlist(df[,1])
  return(t(pa))
  
}
