#' dglm_mod
#'
#' @description Function taken from the dglm package by Peter Dunn and Gordon Smyth. Adapted to be callable by other functions.
#'
dglm_mod <- function (formula = NULL, dformula = ~1, family = stats::gaussian, 
                  dlink = "log", data = NULL, subset = NULL, weights = NULL, 
                  contrasts = NULL, method = "ml", mustart = NULL, betastart = NULL, 
                  etastart = NULL, phistart = NULL,
                  ykeep = TRUE, xkeep = FALSE, zkeep = FALSE, tweedie.var.power=NA, ...) 
{
  call <- match.call()
  if (is.character(family)) {
    family <- get(family, mode = "function", envir = parent.frame())
  }
  if (is.function(family)) {
    family <- family()
  }
  if (is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }
  mean.mframe <- model.frame(formula = formula, data = data)
  y <- stats::model.response(mean.mframe, "numeric")
  if (is.null(dim(y))) {
    N <- length(y)
  }
  else {
    N <- dim(y)[1]
  }
  nobs <- N
  mterms <- attr(mean.mframe, "terms")
  X <- stats::model.matrix(mterms, mean.mframe, contrasts)
  weights <- stats::model.weights(mean.mframe)
  if (is.null(weights)) 
    weights <- rep(1, N)
  if (!is.null(weights) && any(weights < 0)) {
    stop("negative weights not allowed")
  }
  offset <- stats::model.offset(mean.mframe)
  if (is.null(offset)) 
    offset <- rep(0, N)
  if (!is.null(offset) && length(offset) != NROW(y)) {
    stop(gettextf("number of offsets is %d should equal %d (number of observations)", 
                  length(offset), NROW(y)), domain = NA)
  }
  var.mframe <- model.frame(formula = dformula, data = data)
  dterms <- attr(var.mframe, "terms")
  Z <- stats::model.matrix(dterms, var.mframe, contrasts)
  doffset <- stats::model.extract(var.mframe, offset)
  if (is.null(doffset)) 
    doffset <- rep(0, N)
  name.dlink <- substitute(dlink)
  if (is.name(name.dlink)) {
    if (is.character(dlink)) {
      name.dlink <- dlink
    }
    else {
      dlink <- name.dlink <- as.character(name.dlink)
    }
  }
  else {
    if (is.call(name.dlink)) 
      name.dlink <- deparse(name.dlink)
  }
  if (!is.null(name.dlink)) 
    name.dlink <- name.dlink
  if (family$family == "Tweedie") {
    tweedie.p <- tweedie.var.power
  }
  Digamma <- family$family == "Gamma" || (family$family == 
                                            "Tweedie" && tweedie.p == 2)
  if (Digamma) {
    linkinv <- stats::make.link(name.dlink)$linkinv
    linkfun <- stats::make.link(name.dlink)$linkfun
    mu.eta <- stats::make.link(name.dlink)$mu.eta
    valid.eta <- stats::make.link(name.dlink)$valid.eta
    init <- expression({
      if (any(y <= 0)) {
        print(y)
        print(any(y <= 0))
        stop("non-positive values not allowed for the DM gamma family")
      }
      n <- rep.int(1, nobs)
      mustart <- y
    })
    dfamily <- structure(list(family = "Digamma", variance = varfun.digamma, 
                              dev.resids = function(y, mu, wt) {
                                wt * unitdeviance.digamma(y, mu)
                              }, aic = function(y, n, mu, wt, dev) NA, link = name.dlink, 
                              linkfun = linkfun, linkinv = linkinv, mu.eta = mu.eta, 
                              initialize = init, validmu = function(mu) {
                                all(mu > 0)
                              }, valideta = valid.eta))
  }
  else {
    eval(substitute(dfamily <- Gamma(link = lk), list(lk = name.dlink)))
  }
  dlink <- as.character(dfamily$link)
  logdlink <- dlink == "log"
  if (!is.null(method)) {
    name.method <- substitute(method)
    if (!is.character(name.method)) 
      name.method <- deparse(name.method)
    list.methods <- c("ml", "reml", "ML", "REML", "Ml", "Reml")
    i.method <- pmatch(method, list.methods, nomatch = 0)
    if (!i.method) 
      stop("Method must be ml or reml")
    method <- switch(i.method, "ml", "reml", "ml", "reml", 
                     "ml", "reml")
  }
  reml <- method == "reml"
  if (is.null(mustart)) {
    etastart <- NULL
    eval(family$initialize)
    mu <- mustart
    mustart <- NULL
  }
  if (!is.null(betastart)) {
    eta <- X %*% betastart
    mu <- family$linkinv(eta + offset)
  }
  else {
    if (!is.null(mustart)) {
      mu <- mustart
      eta <- family$linkfun(mu) - offset
    }
    else {
      eta <- stats::lm.fit(X, family$linkfun(mu) - offset, 
                           singular.ok = TRUE)$fitted.values
      mu <- family$linkinv(eta + offset)
    }
  }
  d <- family$dev.resids(y, mu, weights)
  if (!is.null(phistart)) {
    phi <- phistart
    deta <- dfamily$linkfun(phi) - doffset
  }
  else {
    deta <- stats::lm.fit(Z, dfamily$linkfun(d + (d == 0)/6) - 
                            doffset, singular.ok = TRUE)$fitted.values
    if (logdlink) 
      deta <- deta + 1.27036
    phi <- dfamily$linkinv(deta + offset)
  }
  if (any(phi <= 0)) {
    warning("Some values for  phi  are non-positive, suggesting an inappropriate model", 
            "Try a different link function.\n")
  }
  zm <- as.vector(eta + (y - mu)/family$mu.eta(eta))
  wm <- as.vector(eval(family$variance(mu)) * weights/phi)
  mfit <- stats::lm.wfit(X, zm, wm, method = "qr", singular.ok = TRUE)
  eta <- mfit$fitted.values
  mu <- family$linkinv(eta + offset)
  if (family$family == "Tweedie") {
    # message("p:", tweedie.p, "\n")
    if ((tweedie.p > 0) & (any(mu < 0))) {
      warning("Some values for  mu  are negative, suggesting an inappropriate model.", 
              "Try a different link function.\n")
    }
  } else {
    tweedie.p=NULL
  }
  d <- family$dev.resids(y, mu, weights)
  const <- dglm.constant(y, family, weights, tweedie.p = tweedie.p)
  if (Digamma) {
    h <- 2 * (lgamma(weights/phi) + (1 + log(phi/weights)) * 
                weights/phi)
  }
  else {
    h <- log(phi/weights)
  }
  m2loglik <- const + sum(h + d/phi)
  if (reml) 
    m2loglik <- m2loglik + 2 * log(abs(prod(diag(mfit$R))))
  m2loglikold <- m2loglik + 1
  epsilon <- 1e-07
  maxit <- 50
  trace <- FALSE
  iter <- 0
  while (abs(m2loglikold - m2loglik)/(abs(m2loglikold) + 1) > 
         epsilon && iter < maxit) {
    hdot <- 1/dfamily$mu.eta(deta)
    if (Digamma) {
      delta <- 2 * weights * (log(weights/phi) - digamma(weights/phi))
      u <- 2 * (weights^2) * (trigamma(weights/phi) - phi/weights)
      fdot <- (phi^2)/u * hdot
    }
    else {
      delta <- phi
      u <- (phi^2)
      fdot <- hdot
    }
    wd <- 1/((fdot^2) * u)
    if (reml) {
      h <- stats::hat(mfit$qr)
      delta <- delta - phi * h
      wd <- wd - 2 * (h/hdot^2/phi^2) + h^2
    }
    if (any(wd < 0)) {
      warning(" Some weights are negative; temporarily fixing.  This may be a sign of an inappropriate model.\n")
      wd[wd < 0] <- 0
    }
    if (any(is.infinite(wd))) {
      warning(" Some weights are negative; temporarily fixing.  This may be a sign of an inappropriate model.\n")
      wd[is.infinite(wd)] <- 100
    }
    zd <- deta + (d - delta) * fdot
    dfit <- stats::lm.wfit(Z, zd, wd, method = "qr", singular.ok = TRUE)
    deta <- dfit$fitted.values
    phi <- dfamily$linkinv(deta + doffset)
    if (any(is.infinite(phi))) {
      warning("*** Some values for  phi  are infinite, suggesting an inappropriate model", 
              "Try a different link function.  Making an attempt to continue...\n")
      phi[is.infinite(phi)] <- 10
    }
    zm <- eta + (y - mu)/family$mu.eta(eta)
    fam.wt <- weights * family$variance(mu)
    wm <- fam.wt/phi
    mfit <- stats::lm.wfit(X, zm, wm, method = "qr", singular.ok = TRUE)
    eta <- mfit$fitted.values
    mu <- family$linkinv(eta + offset)
    if (family$family == "Tweedie") {
      if ((tweedie.p > 0) & (any(mu < 0))) {
        warning("*** Some values for  mu  are negative, suggesting an inappropriate model.", 
                "Try a different link function.  Making an attempt to continue...\n")
        mu[mu <= 0] <- 1
      }
    }
    d <- family$dev.resids(y, mu, weights)
    m2loglikold <- m2loglik
    if (Digamma) {
      h <- 2 * (lgamma(weights/phi) + (1 + log(phi/weights)) * 
                  weights/phi)
    }
    else {
      h <- log(phi/weights)
    }
    m2loglik <- const + sum(h + d/phi)
    if (reml) {
      m2loglik <- m2loglik + 2 * log(abs(prod(diag(mfit$R))))
    }
    iter <- iter + 1
    if (trace) 
      message("DGLM iteration ", iter, ": -2*log-likelihood = ", 
              format(round(m2loglik, 4)), " \n", sep = "")
  }
  mfit$call <- call
  mfit$formula <- formula
  mfit$terms <- mterms
  mfit$model <- mean.mframe
  mfit$family <- family
  mfit$linear.predictors <- mfit$fitted.values + offset
  mfit$fitted.values <- mu
  mfit$prior.weights <- weights
  mfit$contrasts <- attr(X, "contrasts")
  intercept <- attr(mterms, "intercept")
  mfit$df.null <- N - sum(weights == 0) - as.integer(intercept)
  mfit$deviance <- sum(d/phi)
  mfit$aic <- NA
  # mfit$null.deviance <- stats::glm.fit(x = X, y = y, weights = weights/phi, 
  #                                      offset = offset, family = family)
  # if (length(mfit$null.deviance) > 1) 
  #   mfit$null.deviance <- mfit$null.deviance$null.deviance
  if (ykeep) 
    mfit$y <- y
  if (xkeep) 
    mfit$x <- X
  class(mfit) <- c("glm", "lm")
  call$formula <- dformula
  dfit$terms <- dterms
  dfit$model <- var.mframe
  dfit$family <- dfamily
  dfit$prior.weights <- rep(1, N)
  dfit$linear.predictors <- dfit$fitted.values + doffset
  dfit$fitted.values <- phi
  dfit$aic <- NA
  call$dformula <- NULL
  call$family <- call(dfamily$family, link = name.dlink)
  dfit$call <- call
  dfit$residuals <- dfamily$dev.resid(d, phi, wt = rep(1/2, 
                                                       N))
  dfit$deviance <- sum(dfit$residuals)
  # dfit$null.deviance <- stats::glm.fit(x = Z, y = d, weights = rep(1/2, 
  #                                                                  N), offset = doffset, family = dfamily)
  # if (length(dfit$null.deviance) > 1) 
  #   dfit$null.deviance <- dfit$null.deviance$null.deviance
  if (ykeep) 
    dfit$y <- d
  if (zkeep) 
    dfit$z <- Z
  dfit$formula <- as.vector(attr(dterms, "formula"))
  dfit$iter <- iter
  class(dfit) <- c("glm", "lm")
  out <- c(mfit, list(dispersion.fit = dfit, iter = iter, method = method, 
                      m2loglik = m2loglik))
  class(out) <- c("dglm", "glm", "lm")
  out
}


dglm.constant <- function (y, family, weights = 1, tweedie.p=NULL) 
{
  const <- switch(family$family[1], 
                  Gaussian = length(y) * log(2 * pi), 
                  Poisson = 2 * sum(y - y * ifelse(y > 0, log(y), 0) + lgamma(y + 1)), Gamma = 2 * sum(log(y)), 
                  `Inverse Gaussian` = sum(log(2 * pi * y^3)), 
                  Tweedie = switch(match(tweedie.p, c(0, 1, 2, 3), nomatch = 0), 
                                   length(y) * log(2 * pi), 
                                   2 * sum(y - y * ifelse(y > 0, log(y), 0) + lgamma(y + 1)), 
                                   2 * sum(log(y)), sum(log(2 * pi * y^3))), 
                  binomial = -2 *  sum(lgamma(weights + 1) - lgamma(weights * y + 1) - 
                                         lgamma(weights * (1 - y) + 1) + weights * 
                                         (y * ifelse(y > 0, log(y), 0) + (1 - y) * ifelse(1 - y > 0, log(1 - y), 0))
                                       ) + sum(log(weights)))
  if (is.null(const)) {
    V <- family$variance(y)
    if (any(V == 0)) 
      V[V == 0] <- family$variance(y[V == 0] + 1/6)
    const <- sum(log(2 * pi * V))
    if (length(V) == 1 && length(y) > 1) 
      const <- length(y) * const
  }
  const
}

anova.dglm.basic <- function(object, tweedie.power){
  reduced <- dglm_mod(formula = object$formula, 
                        dformula = ~1, 
                        family = object$family, 
                        data = object$model, 
                        tweedie.var.power = tweedie.power)
  Chisq <- reduced$m2loglik - object$m2loglik
  return(tibble::tibble(
    Chisq = Chisq,
    p.value = 1 - stats::pchisq(Chisq, 1)))
}
