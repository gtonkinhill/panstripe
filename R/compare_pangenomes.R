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
#' simA <- simulate_pan(rate=1e-3)
#' simB <- simulate_pan(rate=1e-2)
#' resA <- fit_tweedie(simA$pa, simA$tree)
#' resB <- fit_tweedie(simB$pa, simB$tree)
#' plot_dist_params(list(a=resA, b=resB))
#' comp <- compare_pangenomes(resA, resB)
#'
#' @export
compare_pangenomes <- function(resA, resB){
  
  dat <- purrr::imap_dfr(list(resA, resB), ~{
    tibble::add_column(.x$points, pangenome=LETTERS[[.y]], .before=1)
  })
  
  invisible(suppressWarnings(capture.output({
    ef <- tweedie::tweedie.profile(acc ~ core + pangenome, data = dat, p.vec = seq(1.1,1.9,0.1),
                                   do.smooth = FALSE, method="series")
  })))
  
  m <- stats::glm(acc ~ core*pangenome, data = dat, family = statmod::tweedie(var.power = ef$p.max, link.power = 0))
  a <- stats::aov(m)
  
  A <- purrr::map_dfr(1:nrow(resA$boot_reps), ~{
    convert_tweedie(xi = resA$boot_reps$model.xi[[.x]],
                      mu = exp(resA$boot_reps$model.intercept.estimate[[.x]] + resA$boot_reps$model.core.estimate[[.x]]),
                      phi = resA$boot_reps$model.dispersion.estimate[[.x]])
  })
  
  B <- purrr::map_dfr(1:nrow(resB$boot_reps), ~{
    convert_tweedie(xi = resB$boot_reps$model.xi[[.x]],
                    mu = exp(resB$boot_reps$model.intercept.estimate[[.x]] + resB$boot_reps$model.core.estimate[[.x]]),
                    phi = resB$boot_reps$model.dispersion.estimate[[.x]])
  })
  
  bootdf <- tibble::tibble(
    poisA = A$poisson.lambda,
    poisB = B$poisson.lambda,
    gammaA = A$gamma.mean,
    gammaB = B$gamma.mean
  )
  
  bootstrap_pvalues <- tibble::tibble(
    core=calc_bootp(resA$boot_reps$model.core.estimate, resB$boot_reps$model.core.estimate),
    poisson.mean=calc_bootp(bootdf$poisA, bootdf$poisB),
    gamma.mean=calc_bootp(bootdf$gammaA, bootdf$gammaB)
  )

  return(list(aov=broom::tidy(a),
              bootstrap_pvalues=bootstrap_pvalues)) 
  
}

calc_bootp <- function(a, b){
  intA <- c(min(a), max(a))
  intB <- c(min(b), max(b))
  
  left <- max(intA[[1]], intB[[1]])
  right <- min(intA[[2]], intB[[2]])
  
  return((sum((a<=right) & (a>=left)) + sum((b<=right) & (b>=left)) + 1)/(length(a) + length(b) +1))
}

calc_poisson_mean <- function(x){
  exp((x$model.core.estimate*0+x$model.intercept.estimate)*(2-x$model.xi))/(x$model.dispersion.estimate*(2-x$model.xi))
}

calc_gamma_mean <- function(x){
  (x$model.xi-1)*exp((x$model.core.estimate*0+x$model.intercept.estimate)*(x$model.xi-1))
}
