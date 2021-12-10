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
new_panfit <- function(summary, model, data, bootrap_replicates, 
                       tree, pa){
  
  stopifnot(any(class(summary) %in% c('NULL','tbl','data.frame')))
  stopifnot(any(class(model) %in% c('NULL','glm','glmmTMB')))
  stopifnot(any(class(data) %in% c('NULL','tbl','data.frame')))
  stopifnot(any(class(bootrap_replicates) %in% c('NULL','tbl','data.frame')))
  stopifnot(any(class(tree) %in% c('NULL','phylo')))
  stopifnot(any(class(pa) %in% c('NULL','tbl','data.frame','matrix')))
  
  return(
    structure(
      list(
        summary = summary,
        model = model,
        data = data,
        bootrap_replicates=bootrap_replicates,
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
              'p.value','bootstrap CI (2.5%)',
              'bootstrap CI (97.5%)') %in% colnames(x$summary))) {
      stop(
        "Invalid column names in `summary` data.frame",
        call. = FALSE
      )
    }
  }
  
  # check model
  if(!any(class(x$model) %in% c('NULL','glm','glmmTMB'))) stop("Invalid class for `model`", call. = FALSE)
  
  # check data
  if(!any(class(x$data) %in% c('NULL','tbl','data.frame'))) stop("Invalid class for `data`", call. = FALSE)
  if (!is.null(x$data)){
    if(ncol(x$data)!=4) {
      stop(
        "There must be 4 columns in `data` data.frame",
        call. = FALSE
      )
    }
    if(!all(c('acc','core','istip','height') %in% colnames(x$data))) {
      stop(
        "Invalid column names in `data` data.frame",
        call. = FALSE
      )
    }
  }
  
  # check bootrap_replicates
  if(!any(class(x$bootrap_replicates) %in% c('NULL','tbl','data.frame'))) stop("Invalid class for `bootrap_replicates`", call. = FALSE)
  if (!is.null(x$bootrap_replicates)){
    if(ncol(x$bootrap_replicates)!=11) {
      stop(
        "There must be 12 columns in `bootrap_replicates` data.frame",
        call. = FALSE
      )
    }
    if(!all(c('rep','core','tmean','tpoisson.lambda',
              'tgamma.mean','tgamma.phi','model.xi',
              'model.dispersion.estimate','model.tip.estimate',
              'model.core.estimate','model.height.estimate') %in% colnames(x$bootrap_replicates))) {
      stop(
        "Invalid column names in `bootrap_replicates` data.frame",
        call. = FALSE
      )
    }
  }
  
}