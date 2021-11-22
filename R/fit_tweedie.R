#' fit_tweedie
#'
#' Fits a tweedie...
#' 
#' @description 
#'
#' @param pa a binary gene presence/absence matrix with genes as columns and genomes as rows. The rownames should match the tip.labels of the corresponding phylogeny.
#' @param tree a core gene phylogney of class \link{phylo}
#' @param nboot the number of bootstrap replicates to perform (default=100)
#' @param quiet whether to print progress information (default=FALSE)
#'
#' @return result
#'
#' @examples
#'
#' sim <- simulate_pan(rate=1e-3,mean_trans_size=3)
#' pa <- sim$pa
#' tree <- sim$tree
#' nboot <- 100
#' res <- fit_tweedie(sim$pa, sim$tree)
#'
#' @export
fit_tweedie <- function(pa, tree, nboot=100, quiet=FALSE, stochastic.mapping=FALSE){
  
  index <- which(apply(pa, 2, function(x) length(unique(x)))>1)
  
  sf <- NA
  if (stochastic.mapping){
    sf <- sfreemap(tree, pa = pa[,index], model = "ER", quiet=quiet)
    lmt <- sf$mapped.edge.lmt
    n <- nrow(lmt)-1
    for (i in 1:ncol(lmt)){
      thresh <- sort(lmt[,i],partial=n)[[n]]
      lmt[lmt[,i]<thresh,i] <- -Inf
    }
    
    dat <- data.frame(
      acc=exp(matrixStats::rowLogSumExps(lmt)),
      core=sf$edge.length,
      istip=tree$edge[,2]<=length(tree$tip.label)
    )
  } else {
    # use maximum parsimony
    anc_states <- do.call(cbind, map(index, ~{
      anc <- matrix(0, ncol = 2, nrow = length(tree$tip.label) + tree$Nnode)
      anc[(length(tree$tip.label)+1):nrow(anc),] <- castor::asr_max_parsimony(tree, pa[,.x]+1, Nstates = 2)$ancestral_likelihoods
      anc[cbind(1:length(tree$tip.label),pa[,.x]+1)] <- 1
      return(anc[,1])
    }))
    pmismatch <- do.call(cbind, map(1:nrow(tree$edge), ~{(1-anc_states[tree$edge[.x,1],])*anc_states[tree$edge[.x,2],] +
        anc_states[tree$edge[.x,1],]*(1-anc_states[tree$edge[.x,2],])}))
    dat <- tibble(
      acc=colSums(pmismatch),
      core=tree$edge.length,
      istip=tree$edge[,2]<=length(tree$tip.label)
    )
  }

  invisible(capture.output({ef <- tweedie::tweedie.profile(acc ~ core*istip, data = dat, p.vec = seq(1.1,2,0.1),
                                                  do.smooth = TRUE, method="series", do.ci = TRUE)}))

  m <- glm(acc ~ core*istip, data = dat, family = statmod::tweedie(var.power = ef$p.max, link.power = 0))
  
  if (!quiet) cat('Running bootstrap...\n')
  suppressWarnings({suppressMessages({
    boot_reps <- purrr::map_dfr(1:nboot, ~{
      tdat <- dat[sample(1:nrow(dat), size = nrow(dat), replace = TRUE),]
      invisible(capture.output({tef <- tweedie::tweedie.profile(acc ~ core*istip, data = tdat, p.vec = seq(1.1,1.9,0.1),
                                                       do.smooth = FALSE, method="series")}))
      tm <- glm(acc ~ core*istip, data = tdat, family = statmod::tweedie(var.power = tef$p.max, link.power = 0))
      stm <- summary(tm)
      tp <- predict(tm, 
                    data.frame(core=seq(0, max(dat$core), max(dat$core)/100), istip=FALSE), 
                    type="response")
      
      tout <- convert_tweedie(xi=tef$p.max, mu=tp, phi=tef$phi.max)
      
      if (!quiet) cat(paste0(round(.x / nboot * 100), '% completed\r'))
      if (!quiet & (.x == nboot)) cat('Done\n')
      
      data.frame(
        rep=rep(.x,length(tp)),
        core=seq(0, max(dat$core), max(dat$core)/100),
        tmean=tp,
        tpoisson.lambda=tout$poisson.lambda,
        tgamma.mean=tout$gamma.mean,
        tgamma.phi=tout$gamma.phi,
        model.xi=tef$p.max,
        model.dispersion.estimate=stm$dispersion,
        model.tip.estimate=stm$coefficients[3,1],
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
    boot_reps=br,
    anc=sf)
    )
  
}


convert_tweedie <- function(xi, mu, phi){
  return(list(
    poisson.lambda=(mu^(2-xi))/(phi*(2-xi)),
    gamma.mean=(2-xi)*phi*(mu^(xi-1)),
    gamma.phi=(2-xi)*(xi-1)*(phi^2)*(mu^(2*(xi-1)))
  ))
}
