#' panstripe
#'
#' Fits a Tweedie distribution to the inferred gene gain and loss events along each branch of a given phylogeny. 
#' Includes covariates to control for the impact of annotation errors and the depth of ancestral branches.
#' 
#' @description 
#'
#' @param pa a binary gene presence/absence matrix with genes as columns and genomes as rows. The rownames should match the tip.labels of the corresponding phylogeny.
#' @param tree a core gene phylogeny of class \link{phylo}
#' @param nboot the number of bootstrap replicates to perform (default=100)
#' @param max_height filter out branches with nodes starting above this height (default=NA)
#' @param mrsd most recent sampling date used to calculate height (default= maximum depth of tree)
#' @param quiet whether to print progress information (default=FALSE)
#' @param stochastic.mapping use stochastic mapping in place of maximum parsimony for ancestral state reconstruction (experimental)
#'
#' @return a panfit object with the resulting parameter estimates and bootstrap replicates
#'
#' @examples
#'
#' sim <- simulate_pan(rate=1e-3, mean_trans_size=2, fn_error_rate=2, fp_error_rate=10)
#' pa <- sim$pa
#' tree <- sim$tree
#' nboot <- 10
#' res <- panstripe(sim$pa, sim$tree, nboot=10, max_height = NA)
#' res
#'
#' @export
panstripe <- function(pa, tree, nboot=100, max_height=NA, mrsd=NA, quiet=FALSE, stochastic.mapping=FALSE){
  
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
    invisible(capture.output({ef <- tweedie::tweedie.profile(acc ~ core + height + height:core , data = dat[!dat$istip,,drop=FALSE], 
                                                             p.vec = seq(1.1,2,0.1),
                                                             do.smooth = TRUE, method="series", do.ci = TRUE)}))
  }
  
  if((ef$p.max<1) || (ef$p.max>=2)) stop(paste0('Invalid p.max: ', ef$p.max))
  m <- glm(acc ~ istip*core + height + height:core , data = dat, family = statmod::tweedie(var.power = ef$p.max, link.power = 0))
  
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
        invisible(capture.output({tef <- tweedie::tweedie.profile(acc ~ core + height, data = tdat[!tdat$istip,,drop=FALSE],
                                                                  p.vec = seq(1.1,1.9,0.1),
                                                                  do.smooth = FALSE, method="series")}))
      }

      if((tef$p.max<1) || (tef$p.max>=2)) stop(paste0('Invalid p.max: ', tef$p.max))
      tm <- glm(acc ~ istip*core + height + height:core , data = tdat, family = statmod::tweedie(var.power = tef$p.max, link.power = 0))
      stm <- summary(tm)
      tp <- predict(tm, 
                    data.frame(core=seq(0, max(dat$core), max(dat$core)/100), 
                               height = seq(0, max(dat$core), max(dat$core)/100), 
                               istip=FALSE), 
                    type="response")
      
      tout <- panstripe:::convert_tweedie(xi=tef$p.max, mu=tp, phi=stm$dispersion)
      
      if (!quiet) cat(paste0(round(.x / nboot * 100), '% completed\r'))
      if (!quiet & (.x == nboot)) cat('Done\n')
      variances <- diag(vcov(tm))
      
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
        model.height.estimate=coefficients(tm)['height'],
        converged=tm$converged
      )
      
      
    })
  })})
  
  # summarise results and calculate bootstrap confidence intervals
  s <- broom::tidy(m) %>%
    dplyr::filter(term %in% c('istipTRUE', 'core', 'height'))
  s$term[grepl('istip', s$term)] <- 'tip'
  
  keep <- !duplicated(boot_reps$rep)
  index <- which(boot_reps$rep[keep]==1)
  boot.ci <-  rbind(norm_boot(boot_reps$model.core.estimate[keep], index = index),
                    norm_boot(boot_reps$model.tip.estimate[keep], index = index),
                    norm_boot(boot_reps$model.height.estimate[keep], index = index))
  colnames(boot.ci) <- c('bootstrap CI (2.5%)', 'bootstrap CI (97.5%)')
  
  return(
    new_panfit(
      summary = cbind(s,boot.ci),
      model = m,
      data = dat,
      bootrap_replicates=boot_reps,
      tree=tree,
      pa=pa
    )
  )
  
}


convert_tweedie <- function(xi, mu, phi){
  return(list(
    poisson.lambda=(mu^(2-xi))/(phi*(2-xi)),
    gamma.mean=(2-xi)*phi*(mu^(xi-1)),
    gamma.phi=(2-xi)*(xi-1)*(phi^2)*(mu^(2*(xi-1)))
  ))
}

norm_boot <- function(t, index, conf=0.95){
  t0 <- t[[index]]
  var.t0 <- var(t)
  bias <- mean(t) - t0
  merr <- sqrt(var.t0) * qnorm((1 + conf)/2)
  return(c(t0 - bias - merr, 
           t0 - bias + merr))
}
