#' plot_pangenome_fits
#'
#' @importFrom rlang .data
#'
#' @description Plots the fitted pangenome tweedie regression model.
#'
#' @param fit the result of running the `panstripe` function. Multiple fits can be passed as a named list.
#' @param type either 'cumulative' (default) or 'individual'. Indicates whether or not to plot the cumulative sum of core and accessory branch lengths.
#' @param ci whether or not to include the confidence interval ribbon (default=TRUE)
#' @param plot whether to generate the plot (default) or return a data.frame
#' @param legend toggles the display of the legend on and off
#' @param text_size the base text size of the plot (default=14)
#' @param color_pallete the pallete number passed to `scale_fill_brewer`
#'
#' @return either a ggplot2 object or a `data.frame` with the data needed to recreate the plot
#'
#' @examples
#'
#' sim <- simulate_pan(rate=0)
#' fA <- panstripe(sim$pa, sim$tree, nboot=10, ci_type='perc')
#' plot_pangenome_fits(fA, color_pallete=6, type='individual')
#' plot_pangenome_fits(fA, color_pallete=6)
#' sim <- simulate_pan(rate=1e-2)
#' fB <- panstripe(sim$pa, sim$tree, nboot=10, ci_type='perc')
#' plot_pangenome_fits(list(a=fA,b=fB), color_pallete=6, ci=TRUE, type='individual')
#' plot_pangenome_fits(list(a=fA,b=fB), color_pallete=6, ci=TRUE, type='cumulative')
#'
#' @export
plot_pangenome_fits <- function(fit,
                                type='cumulative',
                                ci=TRUE,
                                plot=TRUE, 
                                legend=TRUE,
                                text_size=14,
                                color_pallete=6){
  # check inputs
  if (class(fit)!='panfit'){
    purrr::map(fit, ~{
      if (class(.x) != 'panfit') stop('fit is not of class `panfit`!')
      panstripe:::validate_panfit(.x)
    })
  } else {
    panstripe:::validate_panfit(fit)
    fit <- list(fit)
  }
  
  if (!type %in% c('cumulative','individual')) stop("type must be one of either 'cumulative' or 'individual'!")
  
  nsim <- 1000
  if (type=='cumulative') {
    
    fit_dat <- purrr::imap_dfr(fit, ~{
      
      # ilink <- family(.x$model)$linkinv
      # p <- predict(.x$model, type = 'link', se.fit = TRUE)
      # df <- tibble::tibble(
      #   pangenome=.y,
      #   core=ape::node.depth.edgelength(.x$tree),
      #   acc=node_depth_edge_weight(.x$tree, ilink(p$fit)),
      #   lower=node_depth_edge_weight(.x$tree, ilink(p$fit - 2*p$se.fit)),
      #   upper=node_depth_edge_weight(.x$tree, ilink(p$fit + 2*p$se.fit))
      # )
      
      d <- .x$model$data
      d$istip <- FALSE
      p <- predict(.x$model, newdata = d, type = 'response')
      df <- purrr::map_dfr(1:nsim, function(i){
        sim <- tweedie::rtweedie(length(.x$model$fitted.values),
                          mu=p,
                          phi = .x$model$phi,
                          power = .x$model$p)
        tibble::tibble(
          edge_index=1:(2*.x$tree$Nnode+1),
          core=ape::node.depth.edgelength(.x$tree),
          simacc=node_depth_edge_weight(.x$tree, sim)
        )
      }) %>%
        dplyr::group_by(edge_index) %>%
        dplyr::summarise(
          core=unique(core),
          acc=quantile(simacc, 0.5),
          lower=quantile(simacc, 0.025),
          upper=quantile(simacc, 0.975)
        ) %>%
        tibble::add_column(pangenome=.y, .before=1) %>%
        dplyr::arrange(core)
      
      df$edge_index <- NULL
      df$acc <- pmax(0, stats::smooth.spline(x = df$core, y = df$acc, cv=TRUE, all.knots=TRUE)$y)
      df$lower <- pmax(0, stats::smooth.spline(x = df$core, y = df$lower, cv=TRUE, all.knots=TRUE)$y)
      df$upper <- pmax(0, stats::smooth.spline(x = df$core, y = df$upper, cv=TRUE, all.knots=TRUE)$y)
      
      return(df)
    })
      
    point_dat <- purrr::imap_dfr(fit, ~{
      temp_tree <- .x$tree
      temp_tree$edge.length <- .x$data$acc
      
      tibble::tibble(
        pangenome=.y,
        acc = ape::node.depth.edgelength(temp_tree),
        core = ape::node.depth.edgelength(.x$tree),
        istip = rep(c(TRUE, FALSE), c(.x$tree$Nnode+1, .x$tree$Nnode))
      )
    })
  } else {
    
    fit_dat <- purrr::imap_dfr(fit, ~{
      
      d <- .x$model$data
      d$istip <- FALSE
      p <- predict(.x$model, newdata = d, type = 'response')
      df <- purrr::map_dfr(1:nsim, function(i){
        sim <- tweedie::rtweedie(length(.x$model$fitted.values), 
                                 mu=p, 
                                 phi = .x$model$phi, 
                                 power = .x$model$p)
        tibble::tibble(
          edge_index=1:length(sim),
          core=.x$model$data$core,
          simacc=sim
        )
      }) %>% 
        dplyr::group_by(edge_index) %>%
        dplyr::summarise(
          core=unique(core),
          acc=quantile(simacc, 0.5),
          lower=quantile(simacc, 0.025),
          upper=quantile(simacc, 0.975)
        ) %>% 
        tibble::add_column(pangenome=.y, .before=1) %>% 
        dplyr::arrange(core)
      
      df$edge_index <- NULL
      df$acc <- pmax(0, stats::smooth.spline(x = df$core, y = df$acc, cv=TRUE, all.knots=TRUE)$y)
      df$lower <- pmax(0, stats::smooth.spline(x = df$core, y = df$lower, cv=TRUE, all.knots=TRUE)$y)
      df$upper <- pmax(0, stats::smooth.spline(x = df$core, y = df$upper, cv=TRUE, all.knots=TRUE)$y)
      
      return(df)
      
      
      
      sims <- do.call(rbind, purrr::map(1:nsim, function(i){
        tweedie::rtweedie(length(.x$model$fitted.values), 
                          mu=.x$model$fitted.values, 
                          phi = .x$model$phi, 
                          power = .x$model$p)
      }))
      
      ci <- apply(sims, 2, function(sim) quantile(sim, c(0.025, 0.975)))
      temp_tree <- .x$tree
      
      df <- tibble::tibble(
        pangenome=.y,
        core=.x$model$data$core,
        acc=colMeans(sims),
        lower=ci[1,],
        upper=ci[2,]) %>% 
        dplyr::arrange(core)
      
      df$acc <- pmax(0, stats::smooth.spline(x = df$core, y = df$acc)$y)
      df$lower <- pmax(0, stats::smooth.spline(x = df$core, y = df$lower)$y)
      df$upper <- pmax(0, stats::smooth.spline(x = df$core, y = df$upper)$y)
      
      return(df)
    })
    
    point_dat <- purrr::imap_dfr(fit, ~{
      .x$data %>% 
        tibble::add_column(pangenome=.y, .before = 1)
    })
    
  }
  
  if (!plot){
    return(list(
      fit_dat=fit_dat,
      point_dat=point_dat
    ))
  }
  
  if (length(fit)>1) {
    gg <- ggplot2::ggplot(fit_dat, ggplot2::aes(x=.data$core, y=.data$acc, colour=.data$pangenome)) +
      ggplot2::geom_line() +
      ggplot2::geom_point(data = point_dat, ggplot2::aes(y=.data$acc, colour=.data$pangenome, shape=istip))
    if (ci){
      gg <- gg + ggplot2::geom_ribbon(ggplot2::aes(ymin=.data$lower, ymax=.data$upper, fill=.data$pangenome), alpha=0.3)
    }
  } else {
    gg <- ggplot2::ggplot(fit_dat, ggplot2::aes(x=.data$core, y=.data$acc)) +
      ggplot2::geom_line() +
      ggplot2::geom_point(data = point_dat, ggplot2::aes(y=.data$acc, shape=istip))
    if (ci){
      gg <- gg + ggplot2::geom_ribbon(ggplot2::aes(ymin=.data$lower, ymax=.data$upper), alpha=0.3)
    }
  }
  
  gg <- gg + 
    ggplot2::scale_colour_brewer(type = 'qual', palette = color_pallete) +
    ggplot2::scale_fill_brewer(type = 'qual', palette = color_pallete) +
    ggplot2::theme_bw(base_size = text_size) +
    ggplot2::xlab('core phylogentic branch distance')
gg
  if (type=='cumulative'){
    gg + ggplot2::ylab('cumulative accessory distance')
  } else {
    gg + ggplot2::ylab('accessory distance')
  }
  
  return(gg)
  
}

node_depth_edge_weight <- function(tree, edge_weight){
  tree$edge.length <- edge_weight
  ape::node.depth.edgelength(tree)
}
