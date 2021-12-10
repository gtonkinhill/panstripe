#' panstripe
#'
#' @description Fits a Tweedie distribution to the inferred gene gain and loss events along each branch of a given phylogeny. 
#' Includes covariates to control for the impact of annotation errors and the depth of ancestral branches.
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
#' sim <- simulate_pan(rate=0, mean_trans_size=2, fn_error_rate=2, fp_error_rate=2)
#' pa <- sim$pa
#' tree <- sim$tree
#' nboot <- 10
#' res <- panstripe(sim$pa, sim$tree, nboot=10, max_height = NA)
#' res$summary
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
  
  # fit model
  if (sum(dat$acc[!dat$istip])<3){
    warning("Very few gene gain/loss events inferred in ancestral branches! Ignoring core term and fitting Poisson!")
    m <- glmmTMB::glmmTMB(acc ~ istip, data = dat, family = poisson)
  } else {
    m <- glmmTMB::glmmTMB(acc ~ istip*core + height + height:core, data = dat, family = glmmTMB::tweedie)
  }
  sm <- summary(m)$coefficients$cond %>% 
    tibble::as_tibble(rownames = 'term')
  
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
        warning("Very few gene gain/loss events inferred in ancestral branches! Ignoring core term and fitting Poisson!")
        tm <- glmmTMB::glmmTMB(acc ~ istip, data = tdat, family = poisson)
        tpower <- NA
        dispersion <- NA
        stm <- summary(tm)$coefficients$cond
      } else {
        tm <- glmmTMB::glmmTMB(acc ~ istip*core + height + height:core, data = tdat, family = glmmTMB::tweedie)
        tpower <- unname(plogis(tm$fit$par["thetaf"]) + 1)
        dispersion <- unname(exp(glmmTMB::fixef(tm)$disp))
        stm <- summary(tm)$coefficients$cond
      }
      
      

      tp <- predict(tm, 
                    data.frame(core=seq(0, max(dat$core), max(dat$core)/100), 
                               height = seq(0, max(dat$core), max(dat$core)/100), 
                               istip=FALSE), 
                    type="response")
      
      tout <- convert_tweedie(xi=tpower,
                              mu=tp,
                              phi=dispersion)
      
      if (!quiet) cat(paste0(round(.x / nboot * 100), '% completed\r'))
      if (!quiet & (.x == nboot)) cat('Done\n')

      data.frame(
        rep=rep(.x,length(tp)),
        core=seq(0, max(dat$core), max(dat$core)/100),
        tmean=tp,
        tpoisson.lambda=tout$poisson.lambda,
        tgamma.mean=tout$gamma.mean,
        tgamma.phi=tout$gamma.phi,
        model.xi=tpower,
        model.dispersion.estimate=dispersion,
        model.tip.estimate=stm[which(rownames(stm)=='istipTRUE'),1],
        model.core.estimate=ifelse('core' %in% rownames(stm), stm[which(rownames(stm)=='core'),1], NA),
        model.height.estimate=ifelse('height' %in% rownames(stm), stm[which(rownames(stm)=='height'),1], NA)
      )

    })
  })})
  
  # summarise results and calculate bootstrap confidence intervals
  sm <- sm[sm$term %in% c('istipTRUE', 'core', 'height'), ,drop=FALSE]
  sm$term[grepl('istip', sm$term)] <- 'tip'
  sm <- sm[match(c('tip', 'core', 'height'), sm$term),]
  colnames(sm) <- c('term','estimate','std.error','statistic','p.value')
  
  keep <- !duplicated(boot_reps$rep)
  index <- which(boot_reps$rep[keep]==1)
  boot.ci <-  rbind(norm_boot(boot_reps$model.tip.estimate[keep], index = index),
                    norm_boot(boot_reps$model.core.estimate[keep], index = index),
                    norm_boot(boot_reps$model.height.estimate[keep], index = index))
  colnames(boot.ci) <- c('bootstrap CI (2.5%)', 'bootstrap CI (97.5%)')
  
  return(
    new_panfit(
      summary = cbind(sm,boot.ci),
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
  var.t0 <- stats::var(t)
  bias <- mean(t) - t0
  merr <- sqrt(var.t0) * stats::qnorm((1 + conf)/2)
  return(c(t0 - bias - merr, 
           t0 - bias + merr))
}
