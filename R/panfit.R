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
#' @param boot_samples data.frame or tibble
#' @param tree phylo
#' @param pa matrix, data.frame or tibble
#' @param anc_states branch level ancestral state matrix
#' @return A panfit object
#'
new_panfit <- function(summary, model, data, boot_samples, 
                       tree, pa, anc_states=NULL){
  
  stopifnot(any(class(summary) %in% c('NULL','tbl','data.frame')))
  stopifnot(any(class(model) %in% c('NULL','glm','glmmTMB')))
  stopifnot(any(class(data) %in% c('NULL','tbl','data.frame')))
  stopifnot(any(class(boot_samples) %in% c('NULL','boot')))
  stopifnot(any(class(tree) %in% c('NULL','phylo')))
  stopifnot(any(class(pa) %in% c('NULL','tbl','data.frame','matrix')))
  
  return(
    structure(
      list(
        summary = summary,
        model = model,
        data = data,
        ci_samples=boot_samples,
        tree=tree,
        pa=pa,
        anc_states=anc_states),
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
    if(!ncol(x$summary) %in% c(7,8)) {
      stop(
        "There must be 7 or 8 columns in `summary` data.frame",
        call. = FALSE
      )
    }
    if(!all(c('term','estimate','p.value','std.error',
              'bootstrap CI 2.5%',
              'bootstrap CI 97.5%') %in% colnames(x$summary))) {
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
    if(!all(c('acc','core','istip') %in% colnames(x$data))) {
      stop(
        "Invalid column names in `data` data.frame",
        call. = FALSE
      )
    }
  }
  
  # check ci samplles
  if(!class(x$boot_samples) %in% c('NULL','boot')) stop("Invalid class for `boot_samples`", call. = FALSE)
  
}

