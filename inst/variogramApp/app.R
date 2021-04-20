#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(geoR)
library(ggplot2)
library(dplyr)
require(raster)
require(tmap)
require("tmaptools")
# source('~/Documents/Lancaster/Lancaster Job/projects/LOALOA/rcode/loaloa/myvariogramplot.R')

################ functions ###########################
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
  if(dim(coords)[1] < 2){
    # if(is.null(uvec))  uvec <- seq(min_dist, maxdist, l = bins)
    
    fit.dat <- data.frame(xx= seq(0, maxdist, length.out = 101), 
                          yy = theo.var.fit(u= seq(0, maxdist, length.out = 101), my.l = my.l))
    maxofy <- max(fit.dat$yy)
    pp <- ggplot() + geom_line(data = fit.dat, aes(x = xx, y = yy), inherit.aes = FALSE) +
      scale_y_continuous(name = "semivariance", limits = c(0,  maxofy))
    return(pp)
    } 
  else{
    if(min_dist < maxdist){
      if(is.null(uvec)) uvec <- seq(min_dist, maxdist, l = bins)
      }else{
        if(is.null(uvec))  uvec <- seq(maxdist, min_dist, l = bins)
        }
  }
  
  empvario <- variog(coords = coords, data = data, uvec = uvec, messages = F)
  # plot(empvario)
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


###### function ##################
geo_sample <- function(xy.all,n,delta,k=0) {
  deltasq<-delta*delta
  N<-dim(xy.all)[1]
  index<-1:N
  index.sample<-sample(index,n,replace=FALSE)
  xy.sample<-xy.all[index.sample,]
  
  for (i in 2:n) {
    dmin<-0
    while (dmin<deltasq) {
      take<-sample(index,1)
      dvec<-(xy.all[take,1]-xy.sample[,1])^2+(xy.all[take,2]-xy.sample[,2])^2
      dmin<-min(dvec)
    }
    
    xy.sample[i,]<-xy.all[take,]
  }
  if(k>0) {
    take<-matrix(sample(1:n,2*k,replace=FALSE),k,2)
    for (j in 1:k) {
      take1<-take[j,1]; take2<-take[j,2]
      xy1<-c(xy.sample[take1,])
      dvec<-(xy.all[,1]-xy1[1])^2+(xy.all[,2]-xy1[2])^2
      neighbour<-order(dvec)[2]
      xy.sample[take2,]<-xy.all[neighbour,]
    }
  }
  
  all.locs <- paste(xy.all[,1], xy.all[,2], sep="_")
  selected.locs <-  paste(xy.sample[,1], xy.sample[,2], sep="_")
  sapply(1:(length(all.locs)), function(i) all.locs[i] %in% selected.locs)
}


##### set the names of the correlation function ##########
# cor.names <- c("matern", "exponential", "gaussian", "spherical", "circular", 
#                "cubic", "wave", "power", "powered.exponential", "cauchy", "gencauchy", 
#                "gneiting", "gneiting.matern", "pure.nugget")
cor.names <- c("matern", "exponential", "gaussian", "spherical", "circular", 
               "cubic", "wave", "power", "powered.exponential", "cauchy", 
               "gneiting",  "pure.nugget")
choices = data.frame(
  var = cor.names,
  num = 1:length(cor.names)
)
# List of choices for selectInput
mylist <- as.list(choices$num)
# Name it
names(mylist) <- choices$var
##################################################### function to output the expresiion for the cor functions
cor.expr <- function(x){
  if(x==1){
    return(withMathJax("Matern covariance function:  
    $$C(u; \\sigma^2 \\phi, \\kappa) = \\sigma^2 (1/(2^{\\kappa-1} \\Gamma(\\kappa))) (u/\\phi)^\\kappa) K_{\\kappa}(u/\\phi)$$
                       $$\\text{where u is the distance, } \\sigma^2 \\text{ is a variance parameter, } \\phi \\text{ is a scale parameter and }$$ $$\\kappa \\text{ is the smoothness parameter}$$"))
  }else if(x==2){
    return(withMathJax("Exponential covariance function:  $$C(u; \\sigma^2, \\phi) = \\sigma^2 \\exp(-u/\\phi)$$  $$\\text{where u is the distance, }
                       \\sigma^2 \\text{ is a variance parameter and } \\phi \\text{ is a scale parameter}$$"))
  }else if(x==3){
    return(withMathJax("Gaussian covariance function:  $$C(u; \\sigma^2, \\phi) = \\sigma^2 \\exp(-u/\\phi)^2$$  $$\\text{where u is the distance, }
                       \\sigma^2 \\text{ is a variance parameter and } \\phi \\text{ is a scale parameter}$$"))
  }else if(x==4){
    return(withMathJax("Spherical covariance function:  $$C(u; \\sigma^2, \\phi) = \\sigma^2  (1-1.5(u/\\phi) + 0.5(u/\\phi)^3) ~~~~ \\text{if } u \\leq \\phi, ~ 0 \\text{ otherwise}$$  
    $$\\text{where u is the distance, } \\sigma^2 \\text{ is a variance parameter and } \\phi \\text{ is a scale parameter}$$"))
  }else if(x==5){
    return(withMathJax("Circular covariance function:  $$C(u; \\sigma^2, \\phi) = \\sigma^2  (1- \\gamma(u)) ~~~~ \\text{if } u \\leq \\phi, ~ 0 \\text{ otherwise}$$ 
    $$\\gamma(\\theta)= 2 ((\\theta * \\sqrt{1-\\theta^2} + \\sin^{-1} \\sqrt{\\theta}))/\\pi \\text{ and } \\theta = min(u/\\phi,1) $$
    $$\\text{where u is the distance, } \\sigma^2 \\text{ is a variance parameter and } \\phi \\text{ is a scale parameter}$$"))
  }else if(x==6){
    return(withMathJax("Cubic covariance function:  $$C(u; \\sigma^2, \\phi) = \\sigma^2  (1- (7((u/\\phi)^2) - 8.75((u/\\phi)^3) + 3.5 ((u/\\phi)^5) - 0.75((u/\\phi)^7)))$$ 
    $$ \\text{if } u \\leq \\phi, ~ 0 \\text{ otherwise}$$ 
    $$\\text{where u is the distance, } \\sigma^2 \\text{ is a variance parameter and } \\phi \\text{ is a scale parameter}$$"))
  }else if(x==7){
    return(withMathJax("Wave covariance function:  $$C(u; \\sigma^2, \\phi) = \\sigma^2  (\\phi/u) \\sin(u/\\phi)$$ 
    $$\\text{where u is the distance, } \\sigma^2 \\text{ is a variance parameter and } \\phi \\text{ is a scale parameter}$$"))
  }else if(x==8){
    return(withMathJax("Power covariance function:  $$C(u; \\sigma^2, \\phi) = \\sigma^2  (A - u^\\phi)$$ 
    $$\\text{where A is a constant, u is the distance, } \\sigma^2 \\text{ is a variance parameter and }$$ 
                       $$\\phi \\text{ is a scale parameter}$$"))
  }else if(x==9){
    return(withMathJax("Powered exponential covariance function:  $$C(u; \\sigma^2, \\phi, \\kappa) = \\sigma^2 \\exp(-u/\\phi)^\\kappa ~~~~  \\text{if } 0 < \\kappa \\leq 2 $$  
    $$\\text{where u is the distance, } \\sigma^2 \\text{ is a variance parameter, } \\phi \\text{ is a scale parameter and }$$ $$\\kappa \\text{ is the smoothness parameter}$$"))
  }else if(x==10){
    return(withMathJax("Cauchy covariance function:  $$C(u; \\sigma^2, \\phi, \\kappa) = \\sigma^2 [1+(u/\\phi)^2]^(-\\kappa)$$  
    $$\\text{where u is the distance, } \\sigma^2 \\text{ is a variance parameter, } \\phi \\text{ is a scale parameter and }$$ $$\\kappa \\text{ is the smoothness parameter}$$"))
  }else if(x==11){
    return(withMathJax("Gneiting covariance function:  $$C(u; \\sigma^2) = \\sigma^2 (1 + 8 s u + 25 s^2 u^2 + 32 s^3 u^3)(1-s u)^8  ~~~~  \\text{if } 0 \\leq s, ~ u \\leq 1, ~ 0 \\text{ otherwise}$$  
    $$\\text{where u is the distance, } \\sigma^2 \\text{ is a variance parameter, }$$"))
  }else if(x==12){
    return(withMathJax("Pure nugget covariance function:  $$C(u; \\sigma^2) = \\sigma^2 k$$  
    $$\\text{where k is a constant value. This model corresponds to no spatial correlation.}$$"))
  }
}




# Define UI for application that draws a histogram
ui <- fluidPage(
  
  # Application title
  titlePanel("Variogram"),
  
  # Sidebar with a slider input for number of bins 
  
  sidebarLayout(
    sidebarPanel(
      

      selectInput(
        inputId = "covmodel",
        label = "Correlation functions",
        choices = mylist,
        selected = 2
      ),
      selectInput("dim", "Dimension",
                  choices = c("One-dimension" = 1,
                    "Two-dimension" = 2),
                  selected =2),
      
      conditionalPanel(condition = "input.dim == 2", 
                       selectInput("domain2", "Domain",
                                   choices = c("100 by 100" = 100,
                                               "150 by 150" = 150,
                                               "150 by 100" = 150100,
                                               "200 by 200" = 200),
                                   selected =150100)),
      conditionalPanel(condition = "input.dim == 1", 
                       selectInput("domain1", "Domain",
                                   choices = c("500" = 500,
                                               "1000" = 1000,
                                               "1500" = 1500),
                                   selected =1000)),
      
      numericInput("sigma", withMathJax("$$\\text{Variance parameter }(\\sigma^2)$$"), 1.8),
      numericInput("phi", withMathJax("$$\\text{Scale parameter }(\\phi)$$"), 10),
      numericInput("tau", "Variance of the nugget effect", 0.02),
      conditionalPanel(
        condition = "input.covmodel == 1 || input.covmodel == 9 || input.covmodel == 10",
        numericInput("kappa", withMathJax("$$\\text{Smoothness parameter }(\\kappa)$$"), 0.5),
      ),
      selectInput("sampling", "Sampling method",
                  c("Random sampling" = "random",
                    "Regular sampling" = "regular")),
      conditionalPanel(
        condition = "input.sampling == 'regular'",
        numericInput("delta", "Minimum distance", 3),
      ),
      numericInput("sample", "Number of samples", 200),

    
      
      sliderInput(inputId = "maxdist",
                  label = "Distance:",
                  min = 0,
                  max = 142,
                  value = 50, step=1),
      
      sliderInput(inputId = "nbins",
                  label = "Number of bins:",
                  min = 0,
                  max = 100,
                  value = 15, step=1),
      
      numericInput("nsim", "number of simulation", 1),
      actionButton("showmap2", "Perform simulation"),
      conditionalPanel(condition = "input.nsim > 1", 
                       sliderInput('myslider', 
                                   'Show animation', 
                                   min=1, 
                                   max=3,  
                                   value=1, 
                                   step=1,
                                   animate=animationOptions(interval = 1000,loop=F)
                                   )
                       )
      

      
      # conditionalPanel(condition = "input.tabselected==1",
      #                  actionButton("showmap", "Perform simulation")),
      # 
      # conditionalPanel(condition = "input.tabselected==2",
      #                  numericInput("nsim", "number of simulation", 5),
      #                  actionButton("showmap2", "Perform simulation"),
      #                  sliderInput('myslider', 
      #                              'Show animation', 
      #                              min=1, 
      #                              max=3, # count all frames(pictures) 
      #                              value=1, 
      #                              step=1,
      #                              animate=animationOptions(interval = 1000,loop=F)
      #                  ))
      
      
    ),
    
    
    
    # Show a plot of the generated distribution
    mainPanel(
      uiOutput("formula"),
      plotOutput("map2")  , 
      plotOutput("plot2")

      # tabsetPanel(type="pills", 
      #             tabPanel("Explore", value = 1,       
      #                      plotOutput("map")  , 
      #                      plotOutput("plot")),
      #             tabPanel("Animation", value = 2,       
      #                      plotOutput("map2")  , 
      #                      plotOutput("plot2")),
      #             id="tabselected"
      # )
      
    )
  )
)

# Define server logic required to draw a histogram
server <- function(input, output, session) {
  # options(shiny.maxRequestSize=200*1024^2)    # to increase the size of the input
  observe({
    updateSliderInput(session, "myslider", min= 1, step= 1, value= 1, max=input$nsim)
  })
  
  observe({
    if(input$dim == 1){
      if(input$domain1==500){
        variog_extent <- 500
      } else if(input$domain1==1000){
        variog_extent <- 1000
      }else if(input$domain1==1500){
        variog_extent <- 1500
      }
    }else{
      if(input$domain2==100){
        variog_extent <- 141
      } else if(input$domain2==150){
        variog_extent <- 212
      }else if(input$domain2==200){
        variog_extent <- 282
      }else if(input$domain2==150100){
        variog_extent <- 180
      }
    }
    
    updateSliderInput(session, "maxdist",
                      label = "Distance:",
                      min = 0,
                      max = variog_extent,
                      value = round(variog_extent/3), step=1)
  })
  

  
  # sim.cont <- eventReactive(input$showmap, {
  #   coords <- expand.grid(x= seq(0, 100, length.out = 50), y= seq(0, 100, length.out = 50))
  #   phi <- input$phi
  #   sigma <- input$sigma
  #   tau <- input$tau
  #   kappa = input$kappa
  #   cov.model = cor.names[as.numeric(input$covmodel)]
  #   U <- as.matrix(dist(coords))
  #   Sigma <- cov.spatial(obj = U, cov.model = cov.model, cov.pars = c(sigma, phi), kappa = kappa)
  #   Y <- t(chol(Sigma))%*%rnorm(nrow(U)) + rnorm(n = nrow(U), sd = sqrt(tau))
  #   simu.data <- data.frame(coords, Y)
  #   # ras.df <- raster::rasterFromXYZ(simu.data)
  #   
  #   ##### Take the subset of the simulated data
  #   if(input$sampling=="random"){
  #     sam.data <- simu.data %>% sample_n(input$sample)
  #   }else{
  #     sam.data <- simu.data[geo_sample(xy.all = simu.data[, 1:2], n = input$sample, delta = input$delta), ]
  #   }
  # 
  #   
  #   ##### map the process and the subset ###############
  #   
  #   col <- rev(terrain.colors(255))
  #   m1 <- ggplot(simu.data) + geom_raster(aes(x = x, y = y, fill=Y)) +  
  #     scale_fill_gradientn(colours = col) + theme(legend.title = element_blank()) + 
  #     geom_point(data = sam.data, aes(x = x, y = y, col=Y))
  # 
  #   res <- list()
  #   # res[["simudata"]] <- simu.data
  #   res[["m1"]] <- m1
  #   
  #   
  #   
  #   ######
  #   my.l <- list()
  #   my.l$cov.model <- cov.model
  #   my.l$cov.pars <- c(sigma, phi)
  #   my.l$nugget <- tau
  #   if (any(my.l$cov.model == c("matern", "powered.exponential", 
  #                               "cauchy", "gencauchy", "gneiting.matern"))) {
  #     my.l$kappa <- kappa
  #   } else my.l$kappa <- NULL
  #   
  #   p1 <- ggvario(coords = sam.data[, 1:2], data = sam.data[,3], theo=T, my.l=my.l, 
  #                 envelop = F, bins = input$nbins, maxdist = input$maxdist)
  #   
  #   res[["p1"]] <- p1
  #   res
  # })
  
  
  # output$map = renderPlot({
  #   sim.cont()$m1
  # })
  # 
  # 
  # output$plot <- renderPlot({
  #   sim.cont()$p1
  # })
  
  
  sim.cont2 <- eventReactive(input$showmap2, {
    if(input$dim == 1){
      coords <- seq(0, as.numeric(input$domain1), by=1)
    }else{
      if(input$domain2 == 150100){
        coords <- expand.grid(x= seq(0, 150, length.out = 50), 
                              y= seq(0, 100, length.out = 50))
      }else{
        coords <- expand.grid(x= seq(0, as.numeric(input$domain2), length.out = as.numeric(input$domain2)/3), 
                              y= seq(0, as.numeric(input$domain2), length.out = as.numeric(input$domain2)/3))
      }
    }

    phi <- input$phi
    sigma <- input$sigma
    tau <- input$tau
    kappa = input$kappa
    cov.model = cor.names[as.numeric(input$covmodel)]
    U <- as.matrix(dist(coords))
    Sigma <- cov.spatial(obj = U, cov.model = cov.model, cov.pars = c(sigma, phi), kappa = kappa)
    col <- rev(terrain.colors(255))
    my.l <- list()
    my.l$cov.model <- cov.model
    my.l$cov.pars <- c(sigma, phi)
    my.l$nugget <- tau
    if (any(my.l$cov.model == c("matern", "powered.exponential", 
                                "cauchy", "gencauchy", "gneiting.matern"))) {
      my.l$kappa <- kappa
    } else my.l$kappa <- NULL
    
    
    m1 <- list()
    p1 <- list()
    for(i in 1:as.numeric(input$nsim)){
      Y <- t(chol(Sigma))%*%rnorm(nrow(U)) + rnorm(n = nrow(U), sd = sqrt(tau))

      if(input$dim == 1){
        simu.data <- data.frame(x=coords, y=coords, Y)
        ##### Take the subset of the simulated data
        if(input$sampling=="random"){
          sam.data <- simu.data %>% sample_n(input$sample)
        }else{
          sam.data <- simu.data[geo_sample(xy.all = simu.data[, 1:2], n = input$sample, delta = input$delta), ]
        }
        ##### map the process and the subset ###############
        m1[[i]] <- ggplot(simu.data) + geom_line(aes(x = x, y = Y)) +  
           theme(legend.title = element_blank()) + 
          geom_point(data = sam.data, aes(x = x, y = Y), size=2, col="red") +
          labs(x="time") 
        p1[[i]] <- ggvario(coords = sam.data[, 1:2], data = sam.data[,3], theo=T, my.l=my.l,
                           envelop = F, bins = input$nbins, maxdist = as.numeric(input$maxdist), xlab = "time")
      }else{
        simu.data <- data.frame(coords, Y)
        ##### Take the subset of the simulated data
        if(input$sampling=="random"){
          sam.data <- simu.data %>% sample_n(input$sample)
        }else{
          sam.data <- simu.data[geo_sample(xy.all = simu.data[, 1:2], n = input$sample, delta = input$delta), ]
        }
        ##### map the process and the subset ###############
        m1[[i]] <- ggplot(simu.data) + geom_raster(aes(x = x, y = y, fill=Y)) +  
          scale_fill_gradientn(colours = col) + theme(legend.title = element_blank()) + 
          geom_point(data = sam.data, aes(x = x, y = y, col=Y), size=2) + coord_fixed()
        p1[[i]] <- ggvario(coords = sam.data[, 1:2], data = sam.data[,3], theo=T, my.l=my.l,
                           envelop = F, bins = input$nbins, maxdist = as.numeric(input$maxdist))

      }
      

    }

    

    

    
    res <- list()
    # res[["simudata"]] <- simu.data
    res[["m1"]] <- m1
    res[["p1"]] <- p1
    res
  })
  
  
  output$map2 = renderPlot({
    sim.cont2()$m1[[as.numeric(input$myslider)]]
  })
  
  
  output$plot2 <- renderPlot({
    sim.cont2()$p1[[as.numeric(input$myslider)]]
  })
  
  output$formula <- renderUI({
    cor.expr(input$covmodel)
  })
  
  
}

# Run the application 
shinyApp(ui = ui, server = server)
