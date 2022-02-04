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
#' @param ci_type the method used to calculate the bootstrap CI (default='bca'). See \link[boot]{boot.ci} for more details.
#' @param conf A scalar indicating the confidence level of the required intervals (default=0.95).
#'
#' @return a panfit object with the resulting parameter estimates and bootstrap replicates
#'
#' @examples
#'
#' sim <- simulate_pan(rate=1e-4, mean_trans_size=5, fn_error_rate=1, fp_error_rate=2)
#' pa <- sim$pa
#' tree <- sim$tree
#' nboot <- 100
#' family <- 'Tweedie'
#' ci_type='perc'
#' boot_type='branch'
#' conf=0.95
#' res <- panstripe(sim$pa, sim$tree, ci_type='norm', boot_type='branch', nboot=100)
#' res$summary
#'
#' @export
panstripe <- function(pa, tree, nboot=1000,
                      quiet=FALSE, stochastic.mapping=FALSE, 
                      min_depth=NULL,
                      family='Tweedie', ci_type='bca', boot_type='branch', conf=0.95){
  
  #check inputs
  if (nrow(pa) != length(tree$tip.label)) stop('number of taxa in tree and p/a matrix need to be the same!')
  if(!all(rownames(pa) %in% tree$tip.label)) stop('taxa labels in p/a matrix and tree do not match!')
  
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
  
  # add depth and filter if requested
  dat$depth <- ape::node.depth.edgelength(tree)[tree$edge[,1]]
  
  if (!is.null(min_depth)){
    dat <- dat[dat$depth>=min_depth,]
  }
  
  # check for all 0's
  if (sum(dat$acc[!dat$istip])==0) {
    warning("No gene gains/losses identified on internal branches! Separation may be a problem.")
  } else if (sum(dat$acc[dat$istip])==0) {
    warning("No gene gains/losses identified at phylogeny tips! Separation may be a problem.")
  }
  
  # fit model
  model <- stats::as.formula("acc ~ istip + core + depth + istip:core")
  
  # fit model
  if (is.character(family) && (family=="Tweedie")){
    m <- fit_tweedie(model, data = dat)
    coef <- c(m$coefficients, m$p, m$phi)
  } else {
    m <- stats::glm(model, data = dat, family=family)
    coef <- c(m$coefficients, NA, NA)
  }
  
  sm <- summary(m)
  coef_names <- strsplit(as.character(model), ' \\+ ')[[3]]
  s <- tibble::tibble(
    term = c('Intercept', coef_names, 'p', 'phi'),
    estimate = coef,
    std.error=c(sm$coefficients[,2], NA, NA),
    statistic=c(sm$coefficients[,3], NA, NA),
    p.value = c(sm$coefficients[,4], NA, NA)
  )
  
  # run bootstrap
  if (nboot>1){
    if (boot_type=='gene'){
      boot_reps <- boot::boot(t(anc_states), fit_model,
                              R = nboot,
                              stype='i',
                              tree=tree, min_depth=min_depth, model=model, family=family, boot_type='gene')
    } else {
      boot_reps <- boot::boot(dat, fit_model,
                              R = nboot,
                              stype='i',
                              tree=tree, min_depth=min_depth, model=model, family=family, boot_type='branch')
    }
    
    # calculate CIs and p-values
    mindex <- length(coef_names) + 1
    if (is.character(family) && (family=="Tweedie")) mindex <- mindex + 2
    ci <- do.call(rbind, purrr::map(1:mindex, ~{
      transformation <- 'identity'
      if (s$term[[.x]]=='p') transformation <- 'logit'
      if (s$term[[.x]]=='phi') transformation <- 'inverse'
      
      boot_ci_pval(boot_reps, index=.x, type=ci_type,
                                       theta_null=0, ci_conf=conf,
                                       transformation=transformation)
      
    }))
    
    if (is.character(family) && (family=="Tweedie")){
      s$`bootstrap CI 2.5%` <- signif(as.numeric(ci[,1]), 5)
      s$`bootstrap CI 97.5%` <- signif(as.numeric(ci[,2]), 5)
    } else {
      s$`bootstrap CI 2.5%` <- c(signif(as.numeric(ci[,1]), 5), NA, NA)
      s$`bootstrap CI 97.5%` <- c(signif(as.numeric(ci[,2]), 5), NA, NA)
    }
  } else {
    s$`bootstrap CI 2.5%` <- NA
    s$`bootstrap CI 97.5%` <- NA
    boot_reps <- NULL
  }
  
  return(
    panstripe:::new_panfit(
      summary = s,
      model = m,
      data = dat,
      boot_samples=boot_reps,
      tree=tree,
      pa=pa
    )
  )
}

twd_llk <- function(p, model, data) {
  tm <- stats::glm(model, data, family = statmod::tweedie(var.power = p, link.power = 0))
  tsm <- summary(tm)
  -sum(log(tweedie::dtweedie(y = tm$y, mu = tm$fitted.values, phi = tsm$dispersion, power = p)))
}

fit_tweedie <- function(model, data){
  op <- stats::optimise(twd_llk, lower = 1, upper = 2, model=model, data=data)
  tm <- stats::glm(model, data, family = statmod::tweedie(var.power = op$minimum, link.power = 0))
  stm <- summary(tm)
  tm$p <- op$minimum
  tm$phi <- stm$dispersion
  return(tm)
}

fit_model <- function(d, indices=1:nrow(d), tree=NULL, min_depth=NULL, model, family, boot_type){
  stopifnot(length(indices)==nrow(d))
  stopifnot(boot_type %in% c('branch', 'gene'))

  if (boot_type == 'branch') {
    tdat <- d[indices,]
  } else {
    tdat <- tibble::tibble(
      acc=colSums(d[indices,tree$edge[,2]]),
      core=tree$edge.length,
      istip=tree$edge[,2]<=length(tree$tip.label) 
    )
    # add depth
    tdat$depth <- ape::node.depth.edgelength(tree)[tree$edge[,2]]
    if (!is.null(min_depth)){
      tdat <- tdat[tdat$min_depth>=min_depth,] 
    }
  }
  
  if (is.character(family) && (family=="Tweedie")){
    tm <- fit_tweedie(model, data = tdat)
    coef <- c(tm$coefficients, tm$p, tm$phi)
  } else {
    tm <- stats::glm(model, data = tdat, family=family)
    coef <- c(tm$coefficients, NA, NA)
  }

  return(coef)
}

boot_ci_pval <- function(boot_res, index, type, 
                         theta_null=0, 
                         precision=NULL, 
                         ci_conf=0.95,
                         transformation='identity',
                         calc_pval=FALSE) {
  
  if (!transformation %in% c('identity','logit','inverse')) stop('Invalid transformation!')
  
  if (transformation=='logit'){
    ll <- stats::make.link('logit')
    h <- function(x) ll$linkfun(x-1)
    hinv <- function(x) ll$linkinv(x)+1
    hdot <- function(x) ll$mu.eta(x)
  } else if (transformation=='inverse'){
    h <- function(x) 1/x
    hinv <- function(x) 1/x
    hdot <- function(x) -1/x^2
  } else {
    h <- function(x) x
    hinv <- function(x) x
    hdot <- function(x) 1
  }
  
  if (is.null(precision)){
    precision = 1/boot_res$R  
  }
  
  if (calc_pval){
    #calculate p-value by inverting the corresponding confidence intervals, as described in Section 3.12 in Hall (1992)
    alpha_seq <- seq(precision, 1 - 1e-16, precision)
    ci <- suppressWarnings({boot::boot.ci(boot_res, conf = 1 - alpha_seq, 
                                          type = type, index=index,
                                          h = h, hinv = hinv, hdot = hdot)})
    bounds <- switch(type, 
                     norm = ci$normal[, 2:3],
                     basic = ci$basic[,4:5],
                     stud = ci$student[, 4:5],
                     perc = ci$percent[, 4:5], 
                     bca = ci$bca[, 4:5])
    alpha <- alpha_seq[which.min(theta_null >= bounds[, 1] & 
                                   theta_null <= bounds[, 2])]
  }
  
  #calculate CI
  ci <- suppressWarnings({boot::boot.ci(boot_res, conf = ci_conf, 
                                        type = type, index=index,
                                        h = h, hinv = hinv, hdot = hdot)})
  bounds <- switch(type, 
                   norm = ci$normal[, 2:3],
                   basic = ci$basic[,4:5],
                   stud = ci$student[, 4:5],
                   perc = ci$percent[, 4:5], 
                   bca = ci$bca[, 4:5])
  
  if (calc_pval){
    return(c(alpha, sort(bounds)))
  } else {
    return(sort(bounds))
  }
}

