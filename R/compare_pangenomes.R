#' compare_pangenomes
#'
#' Compares pangenomes using anova
#' 
#' @description 
#'
#' @param fitA
#' @param fitB
#'
#' @return result
#'
#' @examples
#'
#' simA <- simulate_pan()
#' simB <- simulate_pan(rate=1e-3)
#' resA <- fit_tweedie(simA$pa, simA$tree)
#' resB <- fit_tweedie(simB$pa, simB$tree)
#' comp <- compare_pangenomes(resA, resB)
#'
#' @export
compare_pangenomes <- function(resA, resB, ...){
  
  dat <- purrr::imap_dfr(list(resA, resB, ...), ~{
    tibble::add_column(.x$points, pangenome=LETTERS[[.y]], .before=1)
  })
  
  invisible(suppressWarnings(capture.output({
    ef <- tweedie::tweedie.profile(acc ~ core + pangenome, data = dat, p.vec = seq(1.1,1.9,0.1),
                                   do.smooth = FALSE, method="series")
  })))
  
  m <- stats::glm(acc ~ core + pangenome, data = dat, family = statmod::tweedie(var.power = ef$p.max, link.power = 0))
  
  a <- stats::aov(m)
  
  tuk <- suppressWarnings(stats::TukeyHSD(a, which='pangenome'))
  
  return(list(aov=broom::tidy(a),
              tukey=tuk$pangenome))
  
}