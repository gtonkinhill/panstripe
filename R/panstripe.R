#' panstripe
#'
#' @import cplm
#'
#' @description Fits a Tweedie distribution to the inferred gene gain and loss events along each branch of a given phylogeny. 
#' Includes covariates to control for the impact of annotation errors and the depth of ancestral branches.
#'
#' @param pa a binary gene presence/absence matrix with genes as columns and genomes as rows. The rownames should match the tip.labels of the corresponding phylogeny.
#' @param tree a core gene phylogeny of class \link{phylo}
#' @param nboot the number of bootstrap replicates to perform (default=100)
#' @param quiet whether to print progress information (default=FALSE)
#' @param family the family used by glm. One of 'Tweedie', 'Poisson', 'Gamma' or 'Gaussian'. (default='Tweedie')
#' @param stochastic.mapping use stochastic mapping in place of maximum parsimony for ancestral state reconstruction (experimental)
#' @param method approach used to fit GLM. Can be either 'ML' for Maximum Likelihood or 'Bayesian'.
#'
#' @return a panfit object with the resulting parameter estimates and bootstrap replicates
#'
#' @examples
#'
#' sim <- simulate_pan(rate=0, mean_trans_size=2, fn_error_rate=2, fp_error_rate=2)
#' pa <- sim$pa
#' tree <- sim$tree
#' nboot <- 10
#' family <- 'Tweedie'
#' res <- panstripe(sim$pa, sim$tree, nboot=10)
#' res$summary
#'
#' @export
panstripe <- function(pa, tree, nboot=100,
                      quiet=FALSE, stochastic.mapping=FALSE, 
                      family='Tweedie', method='ML'){
  
  #check inputs
  if (nrow(pa) != length(tree$tip.label)) stop('number of taxa in tree and p/a matrix need to be the same!')
  if(!all(rownames(pa) %in% tree$tip.label)) stop('taxa labels in p/a matrix and tree do not match!')
  if(!method %in% c('Bayesian','ML')) stop("'method' must be either 'ML' or 'Bayesian'")
  if((method=='Bayesian') && (family!="Tweedie")) stop("Bayesian implementation only available for 'Tweedie' family")
  
  if (!(is.character(family) && (family=="Tweedie"))){
    if (is.character(family)) 
      family <- get(family, mode = "function", envir = parent.frame())
    if (is.function(family)) 
      family <- family()
    if (is.null(family$family)) {
      print(family)
      stop("'family' not recognized")
    }
  }
  
  pa <- pa[match(tree$tip.label, rownames(pa)),,drop=FALSE]
  index <- which(apply(pa, 2, function(x) length(unique(x)))>1)
  
  if (stochastic.mapping){
    sf <- panstripe:::sfreemap(tree, pa = pa[,index], model = "ER", quiet=quiet)
    lmt <- sf$mapped.edge.lmt
    lmt[lmt<log(0.5)] <- -Inf
    
    temp_tree <- tree
    temp_tree$edge.length <- exp(matrixStats::rowLogSumExps(lmt))
  } else{
    # use maximum parsimony
    temp_tree <- tree
    anc_states <- do.call(cbind, purrr::map(index, ~{
      return(panstripe:::asr_max_parsimony(tree, pa[,.x]+1, Nstates = 2)$change)
    }))
    temp_tree$edge.length <- rowSums(anc_states)[tree$edge[,2]]
  }
  
  dat <- tibble::tibble(
    acc=ape::node.depth.edgelength(temp_tree),
    core=ape::node.depth.edgelength(tree),
    istip=rep(c(TRUE, FALSE), c(tree$Nnode+1, tree$Nnode))
  )
  
  # check for all 0's
  if (method=="ML"){
    if (sum(dat$acc[!dat$istip])==0) {
      warning("No gene gains/losses identified on internal branches! Separation may be a problem and setting method='Bayesian' method is recommended.")
    } else if (sum(dat$acc[dat$istip])==0) {
      warning("No gene gains/losses identified at phylogeny tips! Separation may be a problem and setting method='Bayesian' method is recommended.")
    }
  }

  # fit model
  formula <- stats::as.formula("acc ~ istip + core")
  
  # fit model
  if (method=="Bayesian"){
    m <- bcplm(formula, data = dat)

    sm <- cbind(m$summary$statistics[,c('Mean','SD')],
                m$summary$quantiles[,c('2.5%','97.5%')])%>%
      tibble::as_tibble(rownames = 'term') %>%
      tibble::add_column(t.value=NA, .before = '2.5%') %>%
      tibble::add_column(p.value=NA, .before = '2.5%')
    
    samples <- purrr::imap_dfr(m$sims.list, ~{
      tibble::as_tibble(.x) %>%
        tibble::add_column(sample=1:nrow(.x), .before=1) %>%
        tibble::add_column(chain=.y, .before=1)
    })
    
  } else {
    if (is.character(family) && (family=="Tweedie")){
      m <- cpglm(formula, data = dat) 
    } else {
      m <- stats::glm(formula, data = dat, family=family)
    }
    
    invisible(capture.output({sm <- summary(m)$coefficients %>% 
      tibble::as_tibble(rownames = 'term')}))
    
    if (!quiet) cat('Running bootstrap...\n')
    boot_reps <- purrr::map_dfr(1:nboot, ~{
      if (.x==1){
        tdat <- dat
      } else {
        tdat <- dat[sample(nrow(dat), size = nrow(dat), replace = TRUE),]
      }
      
      if (is.character(family) && (family=="Tweedie")){
        tm <- cpglm(formula, data = tdat)
        dispersion <- tm$phi
        tpower <- tm$p
      } else {
        tm <- stats::glm(formula, data = tdat, family=family)
        tpower <- NA
        dispersion <- NA
      }
      
      invisible(capture.output({stm <- summary(tm)$coefficients %>% 
        tibble::as_tibble(rownames = 'term')}))

      tp <- predict(tm, 
                    data.frame(core=seq(0, max(dat$core), max(dat$core)/100), 
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
        model.intercept.estimate=stm$Estimate[which(stm$term=='(Intercept)')],
        model.tip.estimate=stm$Estimate[which(stm$term=='istipTRUE')],
        model.core.estimate=stm$Estimate[which(stm$term=='core')]
      )
    })
    
    keep <- !duplicated(boot_reps$rep)
    index <- which(boot_reps$rep[keep]==1)
    samples <-  tibble::as_tibble(rbind(norm_boot(boot_reps$model.intercept.estimate[keep], index = index),
                      norm_boot(boot_reps$model.tip.estimate[keep], index = index),
                      norm_boot(boot_reps$model.core.estimate[keep], index = index)), .name_repair = 'minimal')
    sm <- cbind(sm,samples)
  }
  
  # summarise results and calculate bootstrap confidence intervals
  sm <- sm[sm$term %in% c('istipTRUE', 'core'), ,drop=FALSE]
  sm$term[grepl('istip', sm$term)] <- 'tip'
  sm <- sm[match(c('tip', 'core'), sm$term),]
  colnames(sm) <- c('term','estimate','std.error','statistic','p.value','2.5%','97.5%')
  
  return(
    new_panfit(
      summary = sm,
      model = m,
      data = dat,
      ci_samples=samples,
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
