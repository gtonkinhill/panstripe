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
#' @param boot_type the method used to calculate the bootstrap CI (default='bca'). See \link[boot]{boot.ci} for more details.
#' @param conf A scalar indicating the confidence level of the required intervals (default=0.95).
#'
#' @return a panfit object with the resulting parameter estimates and bootstrap replicates
#'
#' @examples
#'
#' sim <- simulate_pan(rate=1e-3, mean_trans_size=5, fn_error_rate=1, fp_error_rate=1)
#' pa <- sim$pa
#' tree <- sim$tree
#' nboot <- 100
#' family <- 'Tweedie'
#' boot_type='bca'
#' conf=0.95
#' res <- panstripe(sim$pa, sim$tree, boot_type='perc')
#' res$summary
#'
#' @export
panstripe <- function(pa, tree, nboot=1000,
                      quiet=FALSE, stochastic.mapping=FALSE, 
                      family='Tweedie', boot_type='bca', conf=0.95){
  
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
    sf <- panstripe:::sfreemap(tree, pa = pa[,index], model = "ER", quiet=quiet)
    lmt <- sf$mapped.edge.lmt
    lmt[lmt<log(0.5)] <- -Inf
    
    temp_tree <- tree
    temp_tree$edge.length <- exp(matrixStats::rowLogSumExps(lmt))
  } else{
    # use maximum parsimony
    temp_tree <- tree
    temp_tree2 <- tree
    anc_states <- do.call(cbind, purrr::map(index, ~{
      return(panstripe:::asr_max_parsimony(tree, pa[,.x]+1, Nstates = 2)$change)
    }))
    temp_tree$edge.length <- rowSums(anc_states)[tree$edge[,2]]
    temp_tree2$edge.length <- rowSums(anc_states)[tree$edge[,1]]
  }
  
  dat <- tibble::tibble(
    acc=ape::node.depth.edgelength(temp_tree),
    # prev_acc=ape::node.depth.edgelength(temp_tree2),
    core=ape::node.depth.edgelength(tree),
    istip=rep(c(TRUE, FALSE), c(tree$Nnode+1, tree$Nnode))
  )
  
  # check for all 0's
  if (sum(dat$acc[!dat$istip])==0) {
    warning("No gene gains/losses identified on internal branches! Separation may be a problem.")
  } else if (sum(dat$acc[dat$istip])==0) {
    warning("No gene gains/losses identified at phylogeny tips! Separation may be a problem.")
  }
  
  
  # fit model
  model <- stats::as.formula("acc ~ istip + core")
  
  # fit model
  if (is.character(family) && (family=="Tweedie")){
    m <- fit_tweedie(model, data = dat)
    # m <- cplm::cpglm(model, data = dat)
    coef <- c(m$coefficients, m$p, m$phi)
  } else {
    m <- stats::glm(model, data = dat, family=family)
    coef <- c(m$coefficients, NA, NA)
  }
  names(coef) <- c('Intercept', 'tip', 'core', 'p', 'phi')

  # run bootstrap
  boot_reps <- boot::boot(t(anc_states), fit_model,
             R = nboot,
             stype='i',
             tree=tree, model=model, family=family)
  
  # calculate CIs and p-values
  summary <- purrr::map_dfr(1:5, ~{
    transformation <- 'identity'
    if (.x==4) transformation <- 'logit'
    if (.x==5) transformation <- 'inverse'
    
    setNames(c(names(coef)[[.x]], coef[[.x]],
      boot_ci_pval(boot_reps, index=.x, type=boot_type,
                 theta_null=0, ci_conf=conf,
                 transformation=transformation)),
      c('term', 'estimate', 'p-value', '2.5%', '97.5%'))
  
  })
  summary$estimate <- signif(as.numeric(summary$estimate), 4)
  summary$`2.5%` <- signif(as.numeric(summary$`2.5%`), 5)
  summary$`97.5%` <- signif(as.numeric(summary$`97.5%`), 5)
  summary$`p-value`[4:5] <- NA
  
  return(
    panstripe:::new_panfit(
      summary = summary,
      model = m,
      data = dat,
      boot_samples=boot_reps,
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

fit_model <- function(anc, indices=1:ncol(anc), tree, model, family){
  stopifnot(length(indices)==nrow(anc))
  temp_tree <- tree
  temp_tree2 <- tree
  temp_tree$edge.length <- colSums(anc[indices,tree$edge[,2]])
  temp_tree2$edge.length <- colSums(anc[indices,tree$edge[,1]])
  tdat <- tibble::tibble(
    acc=ape::node.depth.edgelength(temp_tree),
    # prev_acc=ape::node.depth.edgelength(temp_tree2),
    core=ape::node.depth.edgelength(tree),
    istip=rep(c(TRUE, FALSE), c(tree$Nnode+1, tree$Nnode))
  )
  
  if (is.character(family) && (family=="Tweedie")){
    tm <- fit_tweedie(model, data = tdat)
    # tm <- cpglm(model, data = tdat)
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
                         transformation='identity') {
  
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
  
  return(c(alpha, sort(bounds)))
}

