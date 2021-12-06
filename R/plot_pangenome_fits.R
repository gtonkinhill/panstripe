#' plot_pangenome_fits
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
#' res <- fit_tweedie(sim$pa, sim$tree, nboot=10)
#' plot_pangenome_fits(list(a=res,b=res), color_pallete=6)
#'
#' @export
plot_pangenome_fits <- function(res, 
                                plot=TRUE, 
                                center=FALSE,
                                legend=TRUE,
                                text_size=14,
                                color_pallete=6){
  
  if (length(res)==5){
    if (all(names(res)==c("points","model_fit","fit_data", "dist_params", "boot_reps"))){
      res <- list(pangeome=res)
    }
  }
  
  fit_dat <- purrr::imap_dfr(res, ~{
    stopifnot(all(names(.x)==c("points","model_fit","fit_data", "dist_params", "boot_reps")))
    if (center){
      offset <- .x$fit_data$val[which.min(.x$fit_data$core)]
      .x$fit_data$val <- .x$fit_data$val - offset
      .x$fit_data$lower <- .x$fit_data$lower - offset
      .x$fit_data$upper <- .x$fit_data$upper - offset
    }
    tibble::add_column(.x$fit_data, pangenome=.y, .before=1)
  })
  
  point_dat <- purrr::imap_dfr(res, ~{
    tibble::add_column(.x$points, pangenome=.y, .before=1)
  })
  
  if (!plot){
    return(list(
      fit_dat=fit_dat,
      point_dat=point_dat
    ))
  }
  
  if (length(res)>1) {
    gg <- ggplot2::ggplot(fit_dat, ggplot2::aes(x=core, y=val, colour=pangenome)) +
      ggplot2::geom_line() +
      ggplot2::geom_ribbon(ggplot2::aes(ymin=lower, ymax=upper, fill=pangenome), alpha=0.3) +
      ggplot2::geom_point(data = point_dat, ggplot2::aes(y=acc, colour=pangenome))
  } else {
    gg <- ggplot2::ggplot(fit_dat, ggplot2::aes(x=core, y=val)) +
      ggplot2::geom_line() +
      ggplot2::geom_ribbon(ggplot2::aes(ymin=lower, ymax=upper), alpha=0.3) +
      ggplot2::geom_point(data = point_dat, ggplot2::aes(y=acc))
  }
  
  gg <- gg + 
    ggplot2::scale_colour_brewer(type = 'qual', palette = color_pallete) +
    ggplot2::scale_fill_brewer(type = 'qual', palette = color_pallete) +
    ggplot2::theme_bw(base_size = text_size) +
    ggplot2::xlab('core phylogentic branch distance') +
    ggplot2::ylab('estimated accessory distance on branch')

  return(gg)
  
}
