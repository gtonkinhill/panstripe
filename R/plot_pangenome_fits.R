#' plot_pangenome_fits
#'
#' Plots the fitted pangenome tweedie regression model.
#' 
#' @description 
#'
#' @param fit the result of running the `panstripe` function. Multiple fits can be passed as a named list.
#' @param boot whether to use estimated SE for confidence intervals (default) or bootstrap estimates.
#' @param plot whether to generate the plot (default) or return a data.frame
#' @param legend toggles the display of the legend on and off
#' @param text_size the base text size of the plot (default=14)
#' @param color_pallete the pallete number passed to `scale_fill_brewer`
#'
#' @return either a ggplot2 object or a `data.frame` with the data needed to recreate the plot
#'
#' @examples
#'
#' sim <- simulate_pan(rate=0)
#' fA <- panstripe(sim$pa, sim$tree, nboot=10)
#' plot_pangenome_fits(fA, color_pallete=6)
#' sim <- simulate_pan(rate=1e-2)
#' fB <- panstripe(sim$pa, sim$tree, nboot=10)
#' plot_pangenome_fits(list(a=fA,b=fB), color_pallete=6, boot=FALSE)
#'
#' @export
plot_pangenome_fits <- function(fit,
                                boot=FALSE,
                                plot=TRUE, 
                                legend=TRUE,
                                text_size=14,
                                color_pallete=6){
  # check inputs
  if (class(fit)!='panfit'){
    purrr::map(fit, ~{
      if (class(.x) != 'panfit') stop('fit is not of class `panfit`!')
      validate_panfit(.x)
    })
  } else {
    validate_panfit(fit)
    fit <- list(fit)
  }
  
  if (boot){
    fit_dat <- purrr::imap_dfr(fit, ~{
      .x$bootrap_replicates %>%
        dplyr::filter(converged) %>%
        dplyr::group_by(core) %>%
        dplyr::summarise(
          val=tmean[which(rep==1)],
          lower=quantile(tmean, 0.025),
          upper=quantile(tmean, 0.975)) %>%
        tibble::add_column(pangenome=.y, .before = 1)
    })
  } else {
    fit_dat <- purrr::imap_dfr(fit, ~{
      p <- predict(.x$model,
                   data.frame(core=seq(0, max(.x$data$core), max(.x$data$core)/100), 
                              height = seq(0, max(.x$data$core), max(.x$data$core)/100), 
                              istip=FALSE), 
                   se.fit = TRUE)
      ilink <- family(.x$model)$linkinv
      
      tibble::tibble(
        pangenome=.y,
        core=seq(0, max(.x$data$core), max(.x$data$core)/100),
        val=ilink(p$fit),
        lower=ilink(p$fit - (2 * p$se.fit)),
        upper=ilink(p$fit + (2 * p$se.fit))
      )
    })
  }
  
  
  point_dat <- purrr::imap_dfr(fit, ~{
    tibble::add_column(.x$data, pangenome=.y, .before=1)
  })
  
  if (!plot){
    return(list(
      fit_dat=fit_dat,
      point_dat=point_dat
    ))
  }
  
  if (length(fit)>1) {
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
