#' plot_acc
#'
#' Plots a pangenome accumulation curve.
#' 
#' @description 
#'
#' @param pa a gene precence/absence matrix. Alternatively a list of matrices can be given.
#' @param text_size adjusts the size of text in the plot
#' @param color_pallete the colour pallete to use. A number between 1 & 9. See 'scale_colour_brewer' for more details
#'
#' @return a plot of the pangenome accumulation curve 
#'
#' @examples
#'
#' sim <- simulate_pan(rate=1e-3)
#' pa <- sim$pa
#' plot_acc(pa, color_pallete=7)
#'
#' @export
plot_acc <- function(pa,
                     nperm=100,
                     legend=TRUE,
                     text_size=14,
                     color_pallete=6){
  
  if (!is.list(pa)){
    pa <- list(pangenome=pa)
  }
  
  plotdf <- purrr::imap_dfr(pa, ~{
    .x <- t(.x)
    purrr::map_dfr(1:nperm, function(i){
      ppa <- .x[sample(nrow(.x), replace = FALSE), sample(ncol(.x), replace = FALSE)]
      cumlative <- rowSums(apply(ppa, 1, cumsum)>0)
      cumlative <- cumlative-cumlative[[1]]
      df <- tibble::tibble(N = 1:length(cumlative),
                           naccessory = cumlative,
                           permutation = i)
      return(df)
    }) %>%
      tibble::add_column(pangenome=.y)
  })
  
  plotdf <- plotdf %>%
    dplyr::group_by(N,pangenome) %>%
    dplyr::summarise(
      `accessory size` = mean(naccessory),
      std = sd(naccessory)
    )
  
  gg <- ggplot2::ggplot(plotdf, ggplot2::aes(N, `accessory size`, col=pangenome, fill=pangenome)) + 
    ggplot2::geom_ribbon(ggplot2::aes(ymin = `accessory size` - std,
                    ymax = `accessory size` + std),
                alpha=0.2, col=NA) + 
    ggplot2::scale_color_brewer(type = 'qual', palette = color_pallete) +
    ggplot2::scale_fill_brewer(type = 'qual', palette = color_pallete) +
    ggplot2::geom_line(size = 1) +
    ggplot2::theme_bw(base_size = text_size) +
    ggplot2::xlab("Number of genomes") +
    ggplot2::ylab("Accessory size")
  
  if (!legend){
    gg <- gg + ggplot2::theme(legend.position = "none")
  }
  
  return(gg)
  
}
