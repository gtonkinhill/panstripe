#' plot_pairwise
#'
#' @importFrom rlang .data
#'
#' @description Plots the pairwise cumulative branch length versus the cumulative number of gene gain and loss events. 
#'
#' @param fit the result of running the `panstripe` function. Multiple fits can be passed as a named list.
#' @param plot whether to generate the plot (default) or return a data.frame
#' @param smooth whether to add a trend line to the plot
#' @param legend toggles the display of the legend on and off
#' @param text_size the base text size of the plot (default=14)
#' @param color_pallete the pallete number passed to `scale_fill_brewer`
#' @param facet whether or not to generate separate plots for each pangenome (default=FALSE).
#' @param smooth adds a linear line to aid the eye in seeing patters. This is NOT the model fit inferred by Panstripe. (default=FALSE)
#'
#' @return either a ggplot2 object or a `data.frame` with the data needed to recreate the plot
#'
#' @examples
#'
#' sim <- simulate_pan(rate=1e-4)
#' fA <- panstripe(sim$pa, sim$tree, nboot=0)
#' plot_pangenome_cumulative(fA, color_pallete=6)
#' sim <- simulate_pan(rate=1e-3)
#' fB <- panstripe(sim$pa, sim$tree, nboot=0)
#' fit <- fB
#' plot=TRUE
#' legend=TRUE
#' text_size=14
#' color_pallete=6
#' facet=FALSE
#' plot_pairwise(list(a=fA,b=fB), color_pallete=6, smooth=TRUE)
#' 
#' @export
plot_pairwise <- function(fit,
                          plot=TRUE, 
                          smooth=FALSE,
                          legend=TRUE,
                          text_size=14,
                          color_pallete=6,
                          facet=FALSE){
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
  
  plot_data <- purrr::imap_dfr(fit, ~{
    
    temp_tree <- .x$tree
    temp_tree$edge.length <- .x$data$acc
    acc_dist <- ape::cophenetic.phylo(temp_tree)
    
    pred_tree <- .x$tree
    pred_tree$edge.length <- predict(.x$model, type = 'response')
    
    return(tibble::tibble(
      pangenome=.y,
      core = rep(ape::cophenetic.phylo(.x$tree)[upper.tri(acc_dist, diag = FALSE)]),
      acc = c(acc_dist[upper.tri(acc_dist, diag = FALSE)]),
      type= rep(c('real'), each=sum(upper.tri(acc_dist, diag = FALSE)))
    ))
    
    # return(tibble::tibble(
    #   pangenome=.y,
    #   core = rep(ape::cophenetic.phylo(.x$tree)[upper.tri(acc_dist, diag = FALSE)], 2),
    #   acc = c(acc_dist[upper.tri(acc_dist, diag = FALSE)],
    #            ape::cophenetic.phylo(pred_tree)[upper.tri(acc_dist, diag = FALSE)]),
    #   type= rep(c('real', 'predicted'), each=sum(upper.tri(acc_dist, diag = FALSE)))
    # ))
  })
  
  
  if (!plot){
    return(plot_data)
  }
  
  if (length(fit)>1) {
    gg <- ggplot2::ggplot(plot_data, ggplot2::aes(x=.data$core, y=.data$acc, colour=.data$pangenome)) +
      ggplot2::geom_point()
    if (facet){
      gg <- gg + ggplot2::facet_wrap(~pangenome, ncol = 1)
    }
  } else {
    gg <- ggplot2::ggplot(plot_data, ggplot2::aes(x=.data$core, y=.data$acc)) +
      ggplot2::geom_point()
  }
  
  if (smooth){
    warning("Adding trend line using 'ggplot2::geom_smooth'. This is not the panstripe fit!")
    gg <- gg + ggplot2::geom_smooth(method ='lm',
                                    level = 0.95)
  }
  
  gg <- gg + 
    ggplot2::scale_colour_brewer(type = 'qual', palette = color_pallete) +
    ggplot2::scale_fill_brewer(type = 'qual', palette = color_pallete) +
    ggplot2::theme_bw(base_size = text_size) +
    ggplot2::xlab('core phylogentic branch distance') + 
    ggplot2::ylab('accessory distance')
  
  gg
  
  return(gg)
  
}
