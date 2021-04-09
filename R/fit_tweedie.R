#' fit_tweedie
#'
#' Fits a tweedie...
#' 
#' @description 
#'
#' @param pa a binary gene presence/absence matrix with genes as columns and genomes as rows. The rownames should match the tip.labels of the corresponding phylogeny.
#' @param tree a core gene phylogney of class \link{phylo}
#' @param nboot the number of bootstrap replicates to perform (default=100)
#'
#' @return result
#'
#' @examples
#'
#' sim <- simulate_pan()
#' pa <- sim$pa
#' tree <- sim$tree
#' nboot <- 100
#' res <- fit_tweedie(sim$pa, sim$tree)
#' 
#'
#' @export
fit_tweedie <- function(pa, tree, nboot=100){
  
  index <- which(apply(pa, 2, function(x) length(unique(x)))>1)
  sf <- sfreemap(tree, pa = pa[,index], model = "ER")
  
  dat <- data.frame(
    acc=rowSums(sf$mapped.edge.lmt),
    core=sf$edge.length
  )
  
  invisible(capture.output({ef <- tweedie::tweedie.profile(acc ~ core, data = dat, p.vec = seq(1.1,1.9,0.1),
                                                  do.smooth = FALSE, method="series")}))

  m <- glm(acc ~ core, data = dat, family = statmod::tweedie(var.power = ef$p.max, link.power = 0))
  
  cat('Running bootstrap...\n')
  suppressWarnings({suppressMessages({
    boot_reps <- purrr::map_dfr(1:nboot, ~{
      tdat <- dat[sample(1:nrow(dat), size = nrow(dat), replace = TRUE),]
      invisible(capture.output({tef <- tweedie::tweedie.profile(acc ~ core, data = tdat, p.vec = seq(1.1,1.9,0.1),
                                                       do.smooth = FALSE, method="series")}))
      # tef <- list(p.max=1.8)
      tm <- glm(acc ~ core, data = tdat, family = statmod::tweedie(var.power = tef$p.max, link.power = 0))
      stm <- summary(tm)
      tp <- predict(tm, 
                    data.frame(core=seq(0, max(dat$core), max(dat$core)/100)), 
                    type="response", se.fit=TRUE)
      tout <- tweedie::tweedie.convert(xi=tef$p.max, mu=tp$fit, phi=tp$residual.scale^2)
      
      cat(paste0(round(.x / nboot * 100), '% completed\r'))
      if (.x == nboot) cat('Done\n')
      
      data.frame(
        rep=rep(.x,length(tp$fit)),
        core=seq(0, max(dat$core), max(dat$core)/100),
        tmean=tp$fit,
        tpoisson.lambda=tout$poisson.lambda,
        tgamma.mean=tout$gamma.mean,
        tgamma.phi=tout$gamma.phi,
        model.xi=tef$p.max,
        model.dispersion.estimate=stm$dispersion,
        model.core.estimate=stm$coefficients[2,1],
        model.intercept.estimate=stm$coefficients[1,1]
      )
      

    })
  })})
  
  fit <- boot_reps %>%
    dplyr::group_by(core) %>%
    dplyr::summarise(
      mean=quantile(tmean, 0.5),
      lower=quantile(tmean, 0.025),
      upper=quantile(tmean, 0.975))
  
  dist_params <- boot_reps %>% 
    dplyr::group_by(core) %>%
    dplyr::summarise(
      "Poisson mean"=quantile(tpoisson.lambda, 0.5),
      "Poisson mean lower"=quantile(tpoisson.lambda, 0.025),
      "Poisson mean upper"=quantile(tpoisson.lambda, 0.975),
      "Gamma mean"=quantile(tgamma.mean, 0.5),
      "Gamma mean lower"=quantile(tgamma.mean, 0.025),
      "Gamma mean upper"=quantile(tgamma.mean, 0.975),
      "Gamma dispersion"=quantile(tgamma.phi, 0.5),
      "Gamma dispersion lower"=quantile(tgamma.phi, 0.025),
      "Gamma dispersion upper"=quantile(tgamma.phi, 0.975)
    )
  
  br <- boot_reps[!duplicated(boot_reps$rep), grepl('model', colnames(boot_reps))]
  
  return(list(
    points = dat,
    model_fit = broom::tidy(m),
    fit_data = fit,
    dist_params = dist_params,
    boot_reps=br)
    )
  
}
