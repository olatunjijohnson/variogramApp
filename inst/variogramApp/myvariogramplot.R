myvariogramplot <-  function (x, max.dist, vario.col = "all", scaled = FALSE, var.lines = FALSE, 
          envelope.obj = NULL, pts.range.cex, bin.cloud = FALSE, ...) 
{
  fc.call <- match.call()
  fc.call$max.dist <- fc.call$vario.col <- fc.call$scaled <- fc.call$pts.range.cex <- fc.call$bin.cloud <- fc.call$envelope.obj <- NULL
  Ldots <- list(...)
  argsnames <- c(names(formals(plot.default)), names(par()))
  names(Ldots) <- argsnames[charmatch(names(Ldots), argsnames)]
  if (missing(max.dist)) 
    max.dist <- max(x$u)
  if (is.null(Ldots$xlim)) 
    fc.call$xlim <- c(0, max.dist)
  if (bin.cloud == TRUE && all(is.na(x$bin.cloud))) 
    stop("plot.variogram: object must be a binned variogram with option bin.cloud=TRUE")
  if (bin.cloud == TRUE && any(!is.na(x$bin.cloud))) 
    boxplot(x$bin.cloud, varwidth = TRUE, xlab = "distance", 
            ylab = paste(x$estimator.type, "variogram"))
  else {
    if (!missing(pts.range.cex)) {
      cex.min <- min(pts.range.cex)
      cex.max <- max(pts.range.cex)
      if (cex.min != cex.max) {
        pts.prop <- TRUE
        sqn <- sqrt(x$n[x$u <= max.dist])
        pts.cex <- cex.min + ((sqn - min(sqn)) * (cex.max - 
                                                    cex.min)/(max(sqn) - min(sqn)))
      }
      else pts.prop <- FALSE
    }
    else pts.prop <- FALSE
    u <- x$u[x$u <= max.dist]
    v <- x$v
    if (is.vector(v) | length(v) == length(x$u)) 
      v <- matrix(v, ncol = 1)
    v <- v[x$u <= max.dist, , drop = FALSE]
    if (vario.col == "all") 
      vario.col <- 1:dim(v)[2]
    else if (mode(vario.col) != "numeric" | any(vario.col > 
                                                ncol(v))) 
      stop("argument vario.col must be equals to \"all\" or a vector indicating the column numbers to be plotted")
    v <- v[, vario.col, drop = FALSE]
    if (scaled) 
      v <- t(t(v)/x$var.mark[vario.col])
    if (is.null(list(...)$ylim)) {
      ymax <- max(v)
      if (!is.null(envelope.obj)) 
        ymax <- max(c(envelope.obj$v.upper, ymax))
      fc.call$ylim <- c(0, ymax)
    }
    if (ncol(v) == 1) {
      fc.call[[1]] <- as.name("plot")
      fc.call$x <- data.frame(distance = u, semivariance = as.vector(v))
      if (pts.prop) 
        fc.call$cex <- pts.cex
      # eval(fc.call, sys.frame(sys.parent()))
      pl <- ggplot(fc.call$x) + geom_point(aes(x = distance, y = semivariance))   # , size=np
    }
    else {
      fc.call[[1]] <- as.name("matplot")
      fc.call$x <- u
      fc.call$y <- v
      if (is.null(Ldots$xlab)) 
        fc.call$xlab <- "distance"
      if (is.null(Ldots$ylab)) 
        fc.call$ylab <- "semivariance"
      if (is.null(Ldots$ty)) {
        if (x$output.type == "bin") 
          fc.call$type <- "p"
        if (x$output.type == "smooth") 
          fc.call$type <- "l"
        if (x$output.type == "cloud") 
          fc.call$type <- "p"
      }
      eval(fc.call, sys.frame(sys.parent()))
    }
    if (var.lines) {
      if (scaled) 
        abline(h = 1, lty = 3)
      else abline(h = x$var.mark, lty = 3)
    }
    if (!is.null(envelope.obj)) {
      lines(u, envelope.obj$v.lower, lty = 4)
      lines(u, envelope.obj$v.upper, lty = 4)
    }
  }
  return(pl)
}

variog_envelope <- function (geodata, coords = geodata$coords, data = geodata$data, 
                             obj.variog, nsim = 99, save.sim = FALSE, messages) 
{
  call.fc <- match.call()
  if (missing(geodata)) 
    geodata <- list(coords = coords, data = data)
  if (missing(messages)) 
    messages.screen <- as.logical(ifelse(is.null(getOption("geoR.messages")), 
                                         TRUE, getOption("geoR.messages")))
  else messages.screen <- messages
  obj.variog$v <- NULL
  if ((is.matrix(data) | is.data.frame(data))) 
    if (ncol(data) > 1) 
      stop("envelops can be computed for only one data set at once")
  if (!is.null(obj.variog$estimator.type)) 
    estimator.type <- obj.variog$estimator.type
  else estimator.type <- "classical"
  if (abs(obj.variog$lambda - 1) > 1e-04) {
    if (abs(obj.variog$lambda) < 1e-04) 
      data <- log(data)
    else data <- ((data^obj.variog$lambda) - 1)/obj.variog$lambda
  }
  xmat <- unclass(trend.spatial(trend = obj.variog$trend, geodata = geodata))
  if (obj.variog$trend != "cte") {
    if (is.vector(data)) {
      data <- lm(data ~ xmat + 0)$residuals
      names(data) <- NULL
    }
    else {
      only.res <- function(y, x) {
        lm(y ~ xmat + 0)$residuals
      }
      data <- apply(data, 2, only.res, x = xmat)
    }
  }
  if (messages.screen) 
    cat(paste("variog.env: generating", nsim, "simulations by permutating data values\n"))
  simula <- list(coords = coords)
  n.data <- length(data)
  perm.f <- function(i, data, n.data) {
    return(data[sample(1:n.data)])
  }
  simula$data <- apply(as.matrix(1:nsim), 1, perm.f, data = data, 
                       n.data = n.data)
  if (messages.screen) 
    cat(paste("variog.env: computing the empirical variogram for the", 
              nsim, "simulations\n"))
  nbins <- length(obj.variog$bins.lim) - 1
  if (obj.variog$direction == "omnidirectional") {
    bin.f <- function(sim) {
      cbin <- vbin <- sdbin <- rep(0, nbins)
      temp <- .C("binit", as.integer(obj.variog$n.data), 
                 as.double(as.vector(coords[, 1])), as.double(as.vector(coords[, 
                                                                               2])), as.double(as.vector(sim)), as.integer(nbins), 
                 as.double(as.vector(obj.variog$bins.lim)), as.integer(estimator.type == 
                                                                         "modulus"), as.double(max(obj.variog$u)), as.double(cbin), 
                 vbin = as.double(vbin), as.integer(FALSE), as.double(sdbin), 
                 PACKAGE = "geoR")$vbin
      return(temp)
    }
    simula.bins <- apply(simula$data, 2, bin.f)
  }
  else {
    variog.vbin <- function(x, ...) {
      variog(geodata = geodata, 
             data = x, uvec = obj.variog$uvec, estimator.type = obj.variog$estimator.type, 
             nugget.tolerance = obj.variog$nugget.tolerance, max.dist = obj.variog$max.dist, 
             pairs.min = obj.variog$pairs.min, direction = obj.variog$direction, 
             tolerance = obj.variog$tolerance, messages.screen = FALSE,...)$v
    }
    simula.bins <- apply(simula$data, 2, variog.vbin)
  }
  simula.bins <- simula.bins[obj.variog$ind.bin, ]
  if (save.sim == FALSE) 
    simula$data <- NULL
  if (messages.screen) 
    cat("variog.env: computing the envelops\n")
  limits <- apply(simula.bins, 1, quantile, prob = c(0.025, 0.975))
  res.env <- list(u = obj.variog$u, v.lower = limits[1, ], 
                  v.upper = limits[2, ])
  if (save.sim) 
    res.env$simulations <- simula$data
  res.env$call <- call.fc
  oldClass(res.env) <- "variogram.envelope"
  return(res.env)
}

# Calculate and plot the variogram
ggvario <- function(coords, 
                    data, 
                    bins = 15, 
                    maxdist = max(dist(coords))/3, 
                    uvec = NULL, 
                    nsim = 999,
                    color = "royalblue1", 
                    xlab = "distance", 
                    show_nbins = F, envelop=F, theo=F, my.l=1) {
  require(geoR)
  coords <- as.matrix(coords)
  min_dist <- min(dist(coords))
  if(is.null(uvec)) uvec <- seq(min_dist, maxdist, l = bins)
  empvario <- variog(coords = coords, data = data, uvec = uvec, messages = F)
  if(envelop==FALSE){
    dfvario <- data.frame(distance = empvario$u, empirical = empvario$v,
                          nbins = empvario$n)
    p1 <- ggplot(dfvario, aes(y = empirical, x = distance, label = nbins)) +
      geom_point(col = "black", fill = color, shape = 21, size = 3) +
      scale_x_continuous(name = xlab, limits = c(0, uvec[length(uvec)]),
                         breaks = round(seq(0, uvec[length(uvec)], l = 6))) +
      # scale_y_continuous(name = "semivariance",
      #                    limits = c(0,  max(dfvario$empirical))) +
      ggtitle("Empirical semivariogram") 
      # theme_classic()
    p2 <- p1 + geom_text(vjust = 1, nudge_y = - diff(range(dfvario$empirical)) / 22)
    
    if(theo == T){
      fit.dat <- data.frame(xx= seq(0, maxdist, length.out = 101), 
                            yy = theo.var.fit(u= seq(0, maxdist, length.out = 101), my.l = my.l))
      maxofy <- max(max(fit.dat$yy), max(dfvario$empirical))
      if(show_nbins){
        p3 <- p2 + geom_line(data = fit.dat , aes(x = xx, y = yy), inherit.aes = FALSE) +
          scale_y_continuous(name = "semivariance", limits = c(0,  maxofy)) 
        p3
      } else{
        p3 <- p1 + geom_line(data = fit.dat, aes(x = xx, y = yy), inherit.aes = FALSE) +
          scale_y_continuous(name = "semivariance", limits = c(0,  maxofy))
        p3
      } 
    }else{
      if(show_nbins){
        p3 <- p2 +  scale_y_continuous(name = "semivariance", limits = c(0,  max(dfvario$empirical))) 
        p3
      } else{
        p3 <- p1 +  scale_y_continuous(name = "semivariance", limits = c(0,  max(dfvario$empirical))) 
        p3
      }
    }
  }else{
    envmc <- variog_envelope(coords = coords, data = data, 
                             obj.variog = empvario, nsim = nsim, messages = F)
    dfvario <- data.frame(distance = empvario$u, empirical = empvario$v,
                          lowemp = envmc$v.lower, upemp = envmc$v.upper, 
                          nbins = empvario$n)
    p1 <- ggplot(dfvario, aes(y = empirical, x = distance, label = nbins)) +
      geom_ribbon(aes(ymin = lowemp, ymax = upemp), fill = color, alpha = .3) +
      geom_point(aes(y = empirical), col = "black", fill = color, shape = 21, size = 3) +
      scale_x_continuous(name = xlab, limits = c(0, uvec[length(uvec)]),
                         breaks = round(seq(0, uvec[length(uvec)], l = 6))) +
      scale_y_continuous(name = "semivariance", 
                         #breaks = round(seq(0, max(dfvario$upemp, dfvario$empirical), l = 6), 1), 
                         limits = c(0, max(dfvario$upemp, dfvario$empirical))) +
      ggtitle("Empirical semivariogram") 
      # theme_classic()
    p2 <- p1 + geom_text(vjust = 1, nudge_y = - diff(range(dfvario$empirical)) / 22)
    if(show_nbins) p2 else p1
  }
}



thr.var <- function (x, max.dist, scaled = FALSE, ...) 
{
  my.l <- list()
  if (missing(max.dist)) {
    my.l$max.dist <- x$max.dist
    if (is.null(my.l$max.dist)) 
      stop("argument max.dist needed for this object")
  }
  else my.l$max.dist <- max.dist
  if (any(x$cov.model == c("matern", "powered.exponential", 
                           "cauchy", "gencauchy", "gneiting.matern"))) 
    my.l$kappa <- x$kappa
  else kappa <- NULL
  if (is.vector(x$cov.pars)) 
    my.l$sill.total <- x$nugget + x$cov.pars[1]
  else my.l$sill.total <- x$nugget + sum(x$cov.pars[, 1])
  my.l$nugget <- x$nugget
  my.l$cov.pars <- x$cov.pars
  my.l$cov.model <- x$cov.model
  if (scaled) {
    if (is.vector(x$cov.model)) 
      my.l$cov.pars[1] <- my.l$cov.pars[1]/my.l$sill.total
    else my.l$cov.pars[, 1] <- my.l$cov.cov.pars[, 1]/my.l$sill.total
    my.l$sill.total <- 1
  }
  gamma.f <- function(x, my.l) {
    if (any(my.l$cov.model == c("linear", "power"))) 
      return(my.l$nugget + my.l$cov.pars[1] * (x^my.l$cov.pars[2]))
    else return(my.l$sill.total - cov.spatial(x, cov.model = my.l$cov.model, 
                                              kappa = my.l$kappa, cov.pars = my.l$cov.pars))
  }
  dd <- gamma.f(x= seq(0, my.l$max.dist, length.out = 101), my.l = my.l)
  # curve(gamma.f(x, my.l = my.l), from = 0, to = my.l$max.dist, 
  #       add = TRUE, ...)
  return(dd)
}

theo.var.fit <- function(u, my.l=NULL, cov.model=NULL,  sigma=NULL, phi=NULL, tau=NULL, kappa=NULL){
  if(is.null(my.l)){
    my.l <- list()
    my.l$cov.model <- cov.model
    my.l$cov.pars <- c(sigma, phi)
    my.l$nugget <- tau
    if (any(my.l$cov.model == c("matern", "powered.exponential", 
                                "cauchy", "gencauchy", "gneiting.matern"))) 
      my.l$kappa <- kappa
    else my.l$kappa <- NULL
  }
  my.l$sill.total <- my.l$nugget + my.l$cov.pars[1]
  gamma.f <- function(x, my.l) {
    if (any(my.l$cov.model == c("linear", "power"))) 
      return(my.l$nugget + my.l$cov.pars[1] * (x^my.l$cov.pars[2]))
    else return(my.l$sill.total - cov.spatial(x, cov.model = my.l$cov.model, 
                                              kappa = my.l$kappa, cov.pars = my.l$cov.pars))
  }
  
  dd <- gamma.f(x= u, my.l = my.l)
  return(dd)
}
