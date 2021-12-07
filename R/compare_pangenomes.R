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
#' fA <- panstripe(simA$pa, simA$tree, nboot=0)
#' fB <- panstripe(simB$pa, simB$tree, nboot=0)
#' comp <- compare_pangenomes(fA, fB)
#'
#' @export
compare_pangenomes <- function(fitA, fitB){
  
  # input checks
  if (class(fitA) != 'panfit') stop('fitA is not of class `panfit`!')
  validate_panfit(fitA)
  if (class(fitB) != 'panfit') stop('fitB is not of class `panfit`!')
  validate_panfit(fitB)
  
  #combine data
  dat <- purrr::imap_dfr(list(fitA, fitB), ~{
    tibble::add_column(.x$data, pangenome=LETTERS[[.y]], .before=1)
  })
  
  # fit model
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

  s <- broom::tidy(a[c('Df','F','Pr(>F)')]) %>%
    dplyr::filter(term %in% c('istip:pangenome', 'core:pangenome'))
  s$term <- ifelse(s$term=='core:pangenome', 'core', 'tip')
  
  return(list(
    summary=s,
    model=m
  ))
  
}
