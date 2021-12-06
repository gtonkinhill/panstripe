#' compare_pangenomes
#'
#' Compares pangenomes using anova
#' 
#' @description 
#'
#' @param fitA
#' @param fitB
#' @param intercept
#'
#' @return result
#'
#' @examples
#'
#' simA <- simulate_pan(rate=0, fp_error_rate=1, fn_error_rate=1)
#' simB <- simulate_pan(rate=0, fp_error_rate=5, fn_error_rate=5)
#' resA <- fit_tweedie(simA$pa, simA$tree, nboot=0)
#' resB <- fit_tweedie(simB$pa, simB$tree, nboot=0)
#' plot_pangenome_fits(list(a=resA, b=resB))
#' comp <- compare_pangenomes(resA, resB)
#'
#' @export
compare_pangenomes <- function(resA, resB){
  
  dat <- purrr::imap_dfr(list(resA, resB), ~{
    tibble::add_column(.x$points, pangenome=LETTERS[[.y]], .before=1)
  })
  
  model <- as.formula("acc ~ core*istip + height + height:core + istip:pangenome + core:pangenome ")
  
  if (sum(dat$acc[!dat$istip])<3){
    warning("No gene gain/loss events inferred in ancestral branches! Setting tweedie power=1")
    ef <- list(p.max=1)
  } else {
    invisible(suppressWarnings(capture.output({
      ef <- tweedie::tweedie.profile(model, 
                                     data = dat[!dat$istip, , drop=FALSE], p.vec = seq(1.1,1.9,0.1),
                                     do.smooth = TRUE, method="series")
    })))
  }
  
  m <- stats::glm(model, data = dat, 
                  family = statmod::tweedie(var.power = ef$p.max, link.power = 0))
  
  a <- stats::anova(m, test="F")

  s <- broom::tidy(a) %>%
    dplyr::filter(term %in% c('istip:pangenome', 'core:pangenome'))
  s$term <- ifelse(s$term=='core:pangenome', 'core', 'tip')
  
  return(list(
    summary=s,
    model=m
  ))
  
}
