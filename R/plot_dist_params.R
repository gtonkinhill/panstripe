#' plot_dist_params
#'
#' Plots the fitted pangenome tweedie regression model.
#' 
#' @description 
#'
#' @param res the result of running the `fit_tweedie` function. Multiple results can be passed as a named list.
#' @param plot (default=FALSE)
#'
#' @return result
#'
#' @examples
#'
#' sim <- simulate_pan(rate=1e-3)
#' res <- fit_tweedie(sim$pa, sim$tree)
#' plot_dist_params(res)
#'
#' @export
plot_dist_params <- function(res, plot=TRUE){
  
  if ((length(res)==5) & all(names(res)==c("points","model_fit","fit_data", "dist_params","boot_reps"))){
    res <- list(pangeome=res)
  }
  
  dist_dat <- purrr::imap_dfr(res, ~{
    stopifnot(all(names(.x)==c("points","model_fit","fit_data", "dist_params","boot_reps")))
    .x <- .x$dist_params
    .x <- rbind(
      tibble::tibble(
        core=.x$core,
        parameter='Poisson mean',
        value=.x$`Poisson mean`,
        lower=.x$`Poisson mean lower`,
        upper=.x$`Poisson mean upper`
      ),
      tibble::tibble(
        core=.x$core,
        parameter='Gamma mean',
        value=.x$`Gamma mean`,
        lower=.x$`Gamma mean lower`,
        upper=.x$`Gamma mean upper`
      )
    )
    .x$pangenome=.y
    return(.x)
  })
  
  if (!plot){
    return(dist_dat)
  }
  
  if (length(res)>1) {
    gg <- ggplot2::ggplot(dist_dat, ggplot2::aes(x=core, y=value, col=pangenome)) +
      ggplot2::geom_line() +
      ggplot2::geom_ribbon(ggplot2::aes(ymin=lower, ymax=upper, fill=pangenome), alpha=0.3) +
      ggplot2::facet_wrap(~parameter, ncol = 2, scales = 'free_y')
  } else {
    gg <- ggplot2::ggplot(dist_dat, ggplot2::aes(x=core, y=value)) +
      ggplot2::geom_line() +
      ggplot2::geom_ribbon(ggplot2::aes(ymin=lower, ymax=upper), alpha=0.3) +
      ggplot2::facet_wrap(~parameter, ncol = 2, scales = 'free_y')
  }
  
  gg <- gg + 
    ggplot2::xlab('core phylogentic branch distance') +
    ggplot2::ylab('estimated parameter value') +
    ggplot2::theme_bw(base_size = 14)
  
  return(gg)
  
}
