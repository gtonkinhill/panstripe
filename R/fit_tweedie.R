#' fit_tweedie
#'
#' Fits a tweedie...
#' 
#' @description 
#'
#' @param pa a binary gene presence/absence matrix with genes as columns and genomes as rows. The rownames should match the tip.labels of the corresponding phylogeny.
#' @param tree a core gene phylogney of class \link{phylo}
#' @param nboot the number of bootstrap replicates to perform (default=100)
#' @param max_height
#' @param mrsd Most recent sampling date
#' @param quiet whether to print progress information (default=FALSE)
#'
#' @return result
#'
#' @examples
#'
#' sim <- simulate_pan(rate=1e-3,mean_trans_size=2,fn_error_rate=,fp_error_rate=10)
#' pa <- sim$pa
#' tree <- sim$tree
#' nboot <- 100
#' res <- fit_tweedie(sim$pa, sim$tree, nboot=10, max_height = NA)
#' res$model_fit
#'
#' @export
fit_tweedie <- function(pa, tree, nboot=100, stochastic.mapping=FALSE, max_height=NA, mrsd=NA, quiet=FALSE){
  
  #check inputs
  if (nrow(pa) != length(tree$tip.label)) stop('number of taxa in tree and p/a matrix need to be the same!')
  if(!all(rownames(pa) %in% tree$tip.label)) stop('taxa labels in p/a matrix and tree do not match!')
  
  pa <- pa[match(tree$tip.label, rownames(pa)),,drop=FALSE]
  index <- which(apply(pa, 2, function(x) length(unique(x)))>1)
  
  if (stochastic.mapping){
    sf <- sfreemap(tree, pa = pa[,index], model = "ER", quiet=quiet)
    lmt <- sf$mapped.edge.lmt
    lmt[lmt<log(0.5)] <- -Inf
    
    dat <- data.frame(
      acc=exp(matrixStats::rowLogSumExps(lmt)),
      core=sf$edge.length,
      istip=tree$edge[,2]<=length(tree$tip.label)
    )
  } else{
    # use maximum parsimony
    anc_states <- do.call(cbind, purrr::map(index, ~{
      return(panstripe:::asr_max_parsimony(tree, pa[,.x]+1, Nstates = 2)$change)
    }))
    dat <- tibble::tibble(
      acc=rowSums(anc_states)[tree$edge[,2]],
      core=tree$edge.length,
      istip=tree$edge[,2]<=length(tree$tip.label) 
    )
  }
  
  #filter out edges that start above height if requested
  depth <- ape::node.depth.edgelength(tree)
  if (is.na(mrsd)){
    mrsd <- max(depth)
  }
  height <- mrsd-depth
  dat$height <- height[tree$edge[,1]]
  if (!is.na(max_height)){
    dat <- dat[height[tree$edge[,1]] <= max_height, ,drop=FALSE]
  }
  
  
  if (sum(dat$acc[!dat$istip])<3){
    warning("No gene gain/loss events inferred in ancestral branches! Setting tweedie power=1")
    ef <- list(p.max=1)
  } else {
    invisible(capture.output({ef <- tweedie::tweedie.profile(acc ~ istip*core + height + height:core , data = dat[!dat$istip,,drop=FALSE], 
                                                             p.vec = seq(1.1,2,0.1),
                                                             do.smooth = TRUE, method="series", do.ci = TRUE)}))
  }
  
  m <- glm(acc ~ istip*core + height + height:core , data = dat, family = statmod::tweedie(var.power = ef$p.max, link.power = 0))
  # spline <- smooth.spline(dat$core, dat$height)
  
  if (!quiet) cat('Running bootstrap...\n')
  suppressWarnings({suppressMessages({
    boot_reps <- purrr::map_dfr(1:nboot, ~{
      if (.x==1){
        tdat <- dat
      } else {
        tdat <- rbind(dat[sample(c(1:nrow(dat))[dat$istip], size = sum(dat$istip), replace = TRUE),],
                      dat[sample(c(1:nrow(dat))[!dat$istip], size = sum(!dat$istip), replace = TRUE),])
      }
      
      if (sum(tdat$acc[!tdat$istip])<3){
        tef <- list(p.max=1)
      } else {
        invisible(capture.output({tef <- tweedie::tweedie.profile(acc ~ istip*core + height + height:core , data = tdat[!tdat$istip,,drop=FALSE],
                                                                  p.vec = seq(1.1,1.9,0.1),
                                                                  do.smooth = FALSE, method="series")}))
      }
      
      tm <- glm(acc ~ istip*core + height + height:core , data = tdat, family = statmod::tweedie(var.power = tef$p.max, link.power = 0))
      stm <- summary(tm)
      tp <- predict(tm, 
                    data.frame(core=seq(0, max(dat$core), max(dat$core)/100), 
                               height = seq(0, max(dat$core), max(dat$core)/100), 
                               istip=FALSE), 
                    type="response")
      
      tout <- convert_tweedie(xi=tef$p.max, mu=tp, phi=stm$dispersion)
      
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
        model.tip.estimate=coefficients(tm)['istipTRUE'],
        model.core.estimate=coefficients(tm)['core'],
        model.intercept.estimate=coefficients(tm)['(Intercept)'],
        converged=tm$converged
      )
      
      
    })
  })})
  
  fit <- boot_reps %>%
    dplyr::filter(converged) %>%
    dplyr::group_by(core) %>%
    dplyr::summarise(
      val=tmean[which(rep==1)],
      lower=quantile(tmean, 0.025),
      upper=quantile(tmean, 0.975))
  
  dist_params <- boot_reps %>% 
    dplyr::filter(converged) %>%
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
    model = m,
    data = dat,
    fit_data = fit,
    dist_params = dist_params,
    boot_reps=br)
  )
  
}


convert_tweedie <- function(xi, mu, phi){
  return(list(
    poisson.lambda=(mu^(2-xi))/(phi*(2-xi)),
    gamma.mean=(2-xi)*phi*(mu^(xi-1)),
    gamma.phi=(2-xi)*(xi-1)*(phi^2)*(mu^(2*(xi-1)))
  ))
}
