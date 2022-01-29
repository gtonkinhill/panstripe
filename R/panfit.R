#' new_panfit
#' 
#' @useDynLib panstripe, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' 
#' @description class for storing output of main panstripe fitting function
#'
#' @param summary data.frame or tibble
#' @param model glm
#' @param data data.frame or tibble
#' @param bootrap_replicates data.frame or tibble
#' @param tree phylo
#' @param pa matrix, data.frame or tibble
#' @return A panfit object
#'
new_panfit <- function(summary, model, data, ci_samples, 
                       tree, pa){
  
  stopifnot(any(class(summary) %in% c('NULL','tbl','data.frame')))
  stopifnot(any(class(model) %in% c('NULL','glm','cplm','cpglm','bcplm')))
  stopifnot(any(class(data) %in% c('NULL','tbl','data.frame')))
  stopifnot(any(class(ci_samples) %in% c('NULL','tbl','data.frame')))
  stopifnot(any(class(tree) %in% c('NULL','phylo')))
  stopifnot(any(class(pa) %in% c('NULL','tbl','data.frame','matrix')))
  
  return(
    structure(
      list(
        summary = summary,
        model = model,
        data = data,
        ci_samples=ci_samples,
        tree=tree,
        pa=pa),
      class = "panfit"
    )
  )
}

#' print.panfit
#' 
#' @description `print.panfit()` prints a summary of pangenome fit.
#'
#' @param x panfit
print.panfit <- function(x){
  print(x$summary)
  invisible(x)
}

#' validate_panfit
#' 
#' @description `validate_panfit()` checks a tibble for internal consistency.
#' Correct behavior can be guaranteed only if this function
#' runs without raising an error.
#'
#' @param x panfit
validate_panfit <- function(x) {
  
  # check summary
  if(!any(class(x$summary) %in% c('NULL','tbl','data.frame'))) stop("Invalid class for `summary`", call. = FALSE)
  if (!is.null(x$summary)){
    if(ncol(x$summary)!=7) {
      stop(
        "There must be 7 columns in `summary` data.frame",
        call. = FALSE
      )
    }
    if(!all(c('term','estimate','std.error','statistic',
              '2.5%',
              '97.5%') %in% colnames(x$summary))) {
      stop(
        "Invalid column names in `summary` data.frame",
        call. = FALSE
      )
    }
  }
  
  # check model
  if(!any(class(x$model) %in% c('NULL','glm','cplm','cpglm','bcplm'))) stop("Invalid class for `model`", call. = FALSE)
  
  # check data
  if(!any(class(x$data) %in% c('NULL','tbl','data.frame'))) stop("Invalid class for `data`", call. = FALSE)
  if (!is.null(x$data)){
    if(ncol(x$data)!=3) {
      stop(
        "There must be 4 columns in `data` data.frame",
        call. = FALSE
      )
    }
    if(!all(c('acc','core','istip') %in% colnames(x$data))) {
      stop(
        "Invalid column names in `data` data.frame",
        call. = FALSE
      )
    }
  }
  
  # check ci samplles
  if(!any(class(x$bootrap_replicates) %in% c('NULL','tbl','data.frame'))) stop("Invalid class for `ci_samples`", call. = FALSE)
  
}