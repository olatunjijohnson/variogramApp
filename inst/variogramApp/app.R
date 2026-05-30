library(shiny)
library(geoR)
library(ggplot2)
library(dplyr)

# --- Helper functions --------------------------------------------------------

# Permutation-based Monte Carlo envelope for the empirical variogram.
# Significance bands are 2.5th/97.5th percentiles across nsim permutations.
variog_envelope <- function(geodata, coords = geodata$coords, data = geodata$data,
                            obj.variog, nsim = 99, save.sim = FALSE, messages) {
  call.fc <- match.call()
  if (missing(geodata))
    geodata <- list(coords = coords, data = data)
  messages.screen <- if (missing(messages))
    as.logical(ifelse(is.null(getOption("geoR.messages")), TRUE, getOption("geoR.messages")))
  else messages
  obj.variog$v <- NULL
  if (is.matrix(data) || is.data.frame(data))
    if (ncol(data) > 1)
      stop("Envelopes can be computed for only one data set at once.")
  estimator.type <- if (!is.null(obj.variog$estimator.type)) obj.variog$estimator.type else "classical"
  if (abs(obj.variog$lambda - 1) > 1e-04) {
    data <- if (abs(obj.variog$lambda) < 1e-04) log(data)
            else ((data^obj.variog$lambda) - 1) / obj.variog$lambda
  }
  xmat <- unclass(trend.spatial(trend = obj.variog$trend, geodata = geodata))
  if (obj.variog$trend != "cte") {
    if (is.vector(data)) {
      data <- lm(data ~ xmat + 0)$residuals
      names(data) <- NULL
    } else {
      data <- apply(data, 2, function(y, x) lm(y ~ xmat + 0)$residuals, x = xmat)
    }
  }
  if (messages.screen)
    cat(paste("variog.env: generating", nsim, "simulations by permuting data values\n"))
  n.data <- length(data)
  simula <- list(
    coords = coords,
    data   = apply(as.matrix(1:nsim), 1, function(i, data, n) data[sample(1:n)],
                   data = data, n = n.data)
  )
  if (messages.screen)
    cat(paste("variog.env: computing the empirical variogram for the", nsim, "simulations\n"))
  nbins <- length(obj.variog$bins.lim) - 1
  if (obj.variog$direction == "omnidirectional") {
    bin.f <- function(sim) {
      cbin <- vbin <- sdbin <- rep(0, nbins)
      .C("binit",
         as.integer(obj.variog$n.data),
         as.double(coords[, 1]), as.double(coords[, 2]),
         as.double(sim), as.integer(nbins),
         as.double(obj.variog$bins.lim),
         as.integer(estimator.type == "modulus"),
         as.double(max(obj.variog$u)), as.double(cbin),
         vbin = as.double(vbin), as.integer(FALSE), as.double(sdbin),
         PACKAGE = "geoR")$vbin
    }
    simula.bins <- apply(simula$data, 2, bin.f)
  } else {
    variog.vbin <- function(x, ...) {
      variog(geodata = geodata, data = x,
             uvec = obj.variog$uvec, estimator.type = obj.variog$estimator.type,
             nugget.tolerance = obj.variog$nugget.tolerance,
             max.dist = obj.variog$max.dist, pairs.min = obj.variog$pairs.min,
             direction = obj.variog$direction, tolerance = obj.variog$tolerance,
             messages.screen = FALSE, ...)$v
    }
    simula.bins <- apply(simula$data, 2, variog.vbin)
  }
  simula.bins <- simula.bins[obj.variog$ind.bin, ]
  if (!save.sim) simula$data <- NULL
  if (messages.screen) cat("variog.env: computing the envelopes\n")
  limits <- apply(simula.bins, 1, quantile, prob = c(0.025, 0.975))
  res.env <- list(u = obj.variog$u, v.lower = limits[1, ], v.upper = limits[2, ])
  if (save.sim) res.env$simulations <- simula$data
  res.env$call <- call.fc
  oldClass(res.env) <- "variogram.envelope"
  res.env
}

# Theoretical variogram values at distances u for a given covariance parameter list.
theo.var.fit <- function(u, my.l = NULL, cov.model = NULL, sigma = NULL,
                         phi = NULL, tau = NULL, kappa = NULL) {
  if (is.null(my.l)) {
    my.l <- list(cov.model = cov.model, cov.pars = c(sigma, phi), nugget = tau,
                 kappa = if (cov.model %in% c("matern", "powered.exponential",
                                              "cauchy", "gencauchy", "gneiting.matern"))
                           kappa else NULL)
  }
  my.l$sill.total <- my.l$nugget + my.l$cov.pars[1]
  gamma.f <- function(x, my.l) {
    if (my.l$cov.model %in% c("linear", "power"))
      my.l$nugget + my.l$cov.pars[1] * (x^my.l$cov.pars[2])
    else
      my.l$sill.total - cov.spatial(x, cov.model = my.l$cov.model,
                                    kappa = my.l$kappa, cov.pars = my.l$cov.pars)
  }
  gamma.f(x = u, my.l = my.l)
}

# Compute and plot the empirical variogram with optional MC envelope and
# theoretical overlay. Returns a ggplot object.
ggvario <- function(coords, data, bins = 15, maxdist = max(dist(coords)) / 3,
                    uvec = NULL, nsim_env = 99, color = "royalblue1",
                    xlab = "distance", show_nbins = FALSE,
                    envelop = FALSE, theo = FALSE, my.l = 1) {
  coords <- as.matrix(coords)
  min_dist <- min(dist(coords))

  # With fewer than 2 unique locations only a theoretical curve can be shown
  if (nrow(coords) < 2) {
    fit.dat <- data.frame(xx = seq(0, maxdist, length.out = 101),
                          yy = theo.var.fit(u = seq(0, maxdist, length.out = 101), my.l = my.l))
    return(
      ggplot() +
        geom_line(data = fit.dat, aes(x = xx, y = yy)) +
        scale_y_continuous(name = "semivariance", limits = c(0, max(fit.dat$yy))) +
        theme_bw()
    )
  }

  if (is.null(uvec))
    uvec <- seq(min(min_dist, maxdist), max(min_dist, maxdist), length.out = bins)

  empvario <- variog(coords = coords, data = data, uvec = uvec, messages = FALSE)

  if (!envelop) {
    dfvario <- data.frame(distance = empvario$u, empirical = empvario$v, nbins = empvario$n)
    p_base <- ggplot(dfvario, aes(y = empirical, x = distance, label = nbins)) +
      geom_point(col = "black", fill = color, shape = 21, size = 3) +
      scale_x_continuous(name = xlab, limits = c(0, uvec[length(uvec)]),
                         breaks = round(seq(0, uvec[length(uvec)], length.out = 6))) +
      ggtitle("Empirical semivariogram") +
      theme_bw()
    if (show_nbins)
      p_base <- p_base + geom_text(vjust = 1, nudge_y = -diff(range(dfvario$empirical)) / 22)
    if (theo) {
      fit.dat <- data.frame(xx = seq(0, maxdist, length.out = 101),
                            yy = theo.var.fit(u = seq(0, maxdist, length.out = 101), my.l = my.l))
      p_base +
        geom_line(data = fit.dat, aes(x = xx, y = yy), inherit.aes = FALSE) +
        scale_y_continuous(name = "semivariance",
                           limits = c(0, max(max(fit.dat$yy), max(dfvario$empirical))))
    } else {
      p_base + scale_y_continuous(name = "semivariance",
                                  limits = c(0, max(dfvario$empirical)))
    }
  } else {
    envmc <- variog_envelope(coords = coords, data = data,
                             obj.variog = empvario, nsim = nsim_env, messages = FALSE)
    dfvario <- data.frame(distance = empvario$u, empirical = empvario$v,
                          lowemp = envmc$v.lower, upemp = envmc$v.upper, nbins = empvario$n)
    p_base <- ggplot(dfvario, aes(y = empirical, x = distance, label = nbins)) +
      geom_ribbon(aes(ymin = lowemp, ymax = upemp), fill = color, alpha = 0.3) +
      geom_point(aes(y = empirical), col = "black", fill = color, shape = 21, size = 3) +
      scale_x_continuous(name = xlab, limits = c(0, uvec[length(uvec)]),
                         breaks = round(seq(0, uvec[length(uvec)], length.out = 6))) +
      scale_y_continuous(name = "semivariance",
                         limits = c(0, max(dfvario$upemp, dfvario$empirical))) +
      ggtitle("Empirical semivariogram") +
      theme_bw()
    if (show_nbins)
      p_base + geom_text(vjust = 1, nudge_y = -diff(range(dfvario$empirical)) / 22)
    else
      p_base
  }
}

# Inhibitory-plus-close-pairs spatial sampling.
# Returns a logical vector of length nrow(xy.all) — TRUE for selected points.
geo_sample <- function(xy.all, n, delta, k = 0) {
  deltasq <- delta^2
  index <- 1:nrow(xy.all)
  xy.sample <- xy.all[sample(index, n, replace = FALSE), ]
  for (i in 2:n) {
    dmin <- 0
    while (dmin < deltasq) {
      take <- sample(index, 1)
      dvec <- (xy.all[take, 1] - xy.sample[, 1])^2 + (xy.all[take, 2] - xy.sample[, 2])^2
      dmin <- min(dvec)
    }
    xy.sample[i, ] <- xy.all[take, ]
  }
  if (k > 0) {
    pairs <- matrix(sample(1:n, 2 * k, replace = FALSE), k, 2)
    for (j in 1:k) {
      xy1  <- c(xy.sample[pairs[j, 1], ])
      dvec <- (xy.all[, 1] - xy1[1])^2 + (xy.all[, 2] - xy1[2])^2
      xy.sample[pairs[j, 2], ] <- xy.all[order(dvec)[2], ]
    }
  }
  all.locs      <- paste(xy.all[, 1],     xy.all[, 2],     sep = "_")
  selected.locs <- paste(xy.sample[, 1],  xy.sample[, 2],  sep = "_")
  all.locs %in% selected.locs
}

# --- Covariance model catalogue ----------------------------------------------

cor.names <- c("matern", "exponential", "gaussian", "spherical", "circular",
               "cubic", "wave", "power", "powered.exponential", "cauchy",
               "gneiting", "pure.nugget")

mylist <- setNames(as.list(seq_along(cor.names)), cor.names)

# Returns a withMathJax UI element with the formula for covariance model x.
cor.expr <- function(x) {
  exprs <- list(
    withMathJax("Matern:
      $$C(u;\\sigma^2,\\phi,\\kappa) = \\frac{\\sigma^2}{2^{\\kappa-1}\\Gamma(\\kappa)}
        \\left(\\frac{u}{\\phi}\\right)^\\kappa K_\\kappa\\!\\left(\\frac{u}{\\phi}\\right)$$
      $$u=\\text{distance},\\; \\sigma^2=\\text{variance},\\;
        \\phi=\\text{scale},\\; \\kappa=\\text{smoothness}$$"),

    withMathJax("Exponential:
      $$C(u;\\sigma^2,\\phi) = \\sigma^2 \\exp(-u/\\phi)$$
      $$u=\\text{distance},\\; \\sigma^2=\\text{variance},\\; \\phi=\\text{scale}$$"),

    withMathJax("Gaussian:
      $$C(u;\\sigma^2,\\phi) = \\sigma^2 \\exp\\!\\left(-(u/\\phi)^2\\right)$$
      $$u=\\text{distance},\\; \\sigma^2=\\text{variance},\\; \\phi=\\text{scale}$$"),

    withMathJax("Spherical:
      $$C(u;\\sigma^2,\\phi) = \\sigma^2\\!\\left(1-\\tfrac{3}{2}\\tfrac{u}{\\phi}
        +\\tfrac{1}{2}\\left(\\tfrac{u}{\\phi}\\right)^3\\right)
        \\text{ if }u\\leq\\phi,\\; 0\\text{ otherwise}$$"),

    withMathJax("Circular:
      $$C(u;\\sigma^2,\\phi) = \\sigma^2(1-\\gamma(u))
        \\text{ if }u\\leq\\phi,\\; 0\\text{ otherwise}$$
      $$\\gamma(\\theta)=\\frac{2}{\\pi}\\!\\left(\\theta\\sqrt{1-\\theta^2}
        +\\sin^{-1}\\theta\\right),\\quad\\theta=\\min(u/\\phi,1)$$"),

    withMathJax("Cubic:
      $$C(u;\\sigma^2,\\phi) = \\sigma^2\\!\\left(1 - 7\\left(\\tfrac{u}{\\phi}\\right)^2
        + 8.75\\left(\\tfrac{u}{\\phi}\\right)^3 - 3.5\\left(\\tfrac{u}{\\phi}\\right)^5
        + 0.75\\left(\\tfrac{u}{\\phi}\\right)^7\\right)
        \\text{ if }u\\leq\\phi,\\; 0\\text{ otherwise}$$"),

    withMathJax("Wave:
      $$C(u;\\sigma^2,\\phi) = \\sigma^2\\frac{\\phi}{u}\\sin\\!\\left(\\frac{u}{\\phi}\\right)$$
      $$u=\\text{distance},\\; \\sigma^2=\\text{variance},\\; \\phi=\\text{scale}$$"),

    withMathJax("Power:
      $$C(u;\\sigma^2,\\phi) = \\sigma^2(A - u^\\phi)$$
      $$A=\\text{constant},\\; \\sigma^2=\\text{variance},\\; \\phi=\\text{scale}$$"),

    withMathJax("Powered exponential:
      $$C(u;\\sigma^2,\\phi,\\kappa) = \\sigma^2\\exp\\!\\left(-(u/\\phi)^\\kappa\\right),
        \\quad 0<\\kappa\\leq 2$$
      $$u=\\text{distance},\\; \\sigma^2=\\text{variance},\\;
        \\phi=\\text{scale},\\; \\kappa=\\text{smoothness}$$"),

    withMathJax("Cauchy:
      $$C(u;\\sigma^2,\\phi,\\kappa) = \\sigma^2\\left[1+\\left(\\frac{u}{\\phi}\\right)^2\\right]^{-\\kappa}$$
      $$u=\\text{distance},\\; \\sigma^2=\\text{variance},\\;
        \\phi=\\text{scale},\\; \\kappa=\\text{smoothness}$$"),

    withMathJax("Gneiting:
      $$C(u;\\sigma^2) = \\sigma^2(1+8su+25s^2u^2+32s^3u^3)(1-su)^8
        \\text{ if }su\\leq 1,\\; 0\\text{ otherwise}$$"),

    withMathJax("Pure nugget:
      $$C(u;\\sigma^2) = \\sigma^2 k$$
      $$\\text{No spatial correlation — entirely independent noise.}$$")
  )
  exprs[[as.integer(x)]]
}

# --- UI ----------------------------------------------------------------------

ui <- fluidPage(
  titlePanel("Gaussian Random Field & Variogram Explorer"),

  sidebarLayout(
    sidebarPanel(
      selectInput("covmodel", "Covariance model", choices = mylist, selected = 2),

      selectInput("dim", "Dimension",
                  choices = c("One-dimension" = 1, "Two-dimension" = 2), selected = 2),

      conditionalPanel(
        condition = "input.dim == 2",
        selectInput("domain2", "Domain",
                    choices = c("100 × 100" = 100, "150 × 150" = 150,
                                "150 × 100" = 150100, "200 × 200" = 200),
                    selected = 150100)
      ),
      conditionalPanel(
        condition = "input.dim == 1",
        selectInput("domain1", "Domain length",
                    choices = c("500" = 500, "1000" = 1000, "1500" = 1500),
                    selected = 1000)
      ),

      numericInput("sigma", withMathJax("$$\\text{Variance } (\\sigma^2)$$"),     1.8,  min = 1e-4),
      numericInput("phi",   withMathJax("$$\\text{Scale } (\\phi)$$"),            10,   min = 1e-4),
      numericInput("tau",   withMathJax("$$\\text{Nugget variance } (\\tau^2)$$"), 0.02, min = 0),
      conditionalPanel(
        condition = "input.covmodel == 1 || input.covmodel == 9 || input.covmodel == 10",
        numericInput("kappa", withMathJax("$$\\text{Smoothness } (\\kappa)$$"), 0.5, min = 1e-4)
      ),

      selectInput("sampling", "Sampling design",
                  choices = c("Random" = "random", "Inhibitory (min distance)" = "regular")),
      conditionalPanel(
        condition = "input.sampling == 'regular'",
        numericInput("delta", "Minimum distance between samples", 3, min = 0)
      ),
      numericInput("sample", "Number of sample locations", 100, min = 2),

      sliderInput("maxdist", "Max variogram distance:", min = 0, max = 142, value = 50, step = 1),
      sliderInput("nbins",   "Number of bins:",         min = 2, max = 100, value = 15, step = 1),

      checkboxInput("show_envelop", "Show Monte Carlo envelope (slower)", value = FALSE),
      checkboxInput("show_nbins",   "Show bin counts on plot",            value = FALSE),

      numericInput("nsim", "Number of simulations", 1, min = 1),
      actionButton("showmap2", "Run simulation", class = "btn-primary"),

      conditionalPanel(
        condition = "input.nsim > 1",
        sliderInput("myslider", "Browse simulations:",
                    min = 1, max = 3, value = 1, step = 1,
                    animate = animationOptions(interval = 1000, loop = FALSE))
      ),

      br(),
      downloadButton("downloadPlot", "Download variogram plot")
    ),

    mainPanel(
      uiOutput("formula"),
      plotOutput("map2"),
      plotOutput("plot2")
    )
  )
)

# --- Server ------------------------------------------------------------------

server <- function(input, output, session) {

  # Keep the animation slider range in sync with nsim
  observe({
    updateSliderInput(session, "myslider", min = 1, step = 1, value = 1, max = input$nsim)
  })

  # Update max-distance slider when domain or dimension changes
  observe({
    variog_extent <- if (input$dim == 1) {
      as.numeric(input$domain1)
    } else {
      # Approximate diagonal of each 2-D domain
      switch(as.character(input$domain2),
             "100"    = 141,
             "150"    = 212,
             "200"    = 282,
             "150100" = 180)
    }
    updateSliderInput(session, "maxdist",
                      label = "Max variogram distance:",
                      min = 0, max = variog_extent,
                      value = round(variog_extent / 3), step = 1)
  })

  sim.results <- eventReactive(input$showmap2, {
    validate(
      need(isTRUE(input$sigma > 0),  "Variance σ² must be positive."),
      need(isTRUE(input$phi   > 0),  "Scale ϕ must be positive."),
      need(isTRUE(input$tau   >= 0), "Nugget τ² must be non-negative."),
      need(isTRUE(input$sample >= 2), "Number of samples must be at least 2.")
    )

    phi       <- input$phi
    sigma     <- input$sigma
    tau       <- input$tau
    kappa     <- input$kappa
    cov.model <- cor.names[as.numeric(input$covmodel)]
    nsim      <- as.numeric(input$nsim)
    col       <- rev(terrain.colors(255))

    my.l <- list(
      cov.model = cov.model,
      cov.pars  = c(sigma, phi),
      nugget    = tau,
      kappa     = if (cov.model %in% c("matern", "powered.exponential",
                                       "cauchy", "gencauchy", "gneiting.matern"))
                    kappa else NULL
    )

    withProgress(message = "Simulating random field...", value = 0, {

      # Simulate all realisations at once via grf
      if (input$dim == 1) {
        n1d <- as.numeric(input$domain1)
        grid1d <- cbind(x = seq(0, n1d, by = 1), y = seq(0, n1d, by = 1))
        Y <- as.matrix(grf(grid = grid1d, cov.model = cov.model,
                           cov.pars = c(sigma, phi), nugget = tau,
                           nsim = nsim, kappa = kappa)$data)
        coords <- seq(0, n1d, by = 1)
      } else if (input$domain2 == 150100) {
        grid2d <- expand.grid(x = seq(0, 150, length.out = 50),
                              y = seq(0, 100, length.out = 50))
        Y      <- as.matrix(grf(grid = grid2d, cov.model = cov.model,
                                cov.pars = c(sigma, phi), nugget = tau,
                                nsim = nsim, kappa = kappa)$data)
        coords <- grid2d
      } else {
        d      <- as.numeric(input$domain2)
        grid2d <- expand.grid(x = seq(0, d, length.out = d / 3),
                              y = seq(0, d, length.out = d / 3))
        Y      <- as.matrix(grf(grid = grid2d, cov.model = cov.model,
                                cov.pars = c(sigma, phi), nugget = tau,
                                nsim = nsim, kappa = kappa)$data)
        coords <- grid2d
      }

      m1 <- vector("list", nsim)
      p1 <- vector("list", nsim)

      for (i in seq_len(nsim)) {
        incProgress(1 / nsim, detail = paste("Simulation", i, "of", nsim))

        if (input$dim == 1) {
          simu.data <- data.frame(x = coords, y = coords, Y = Y[, i])
          sam.data <- if (input$sampling == "random")
            simu.data %>% sample_n(input$sample)
          else
            simu.data[geo_sample(xy.all = as.matrix(simu.data[, 1:2]),
                                 n = input$sample, delta = input$delta), ]

          m1[[i]] <- ggplot(simu.data, aes(x = x, y = Y)) +
            geom_line() +
            geom_point(data = sam.data, aes(x = x, y = Y), size = 2, col = "red") +
            labs(x = "time", y = "Y") +
            theme_bw()

          p1[[i]] <- ggvario(coords = sam.data[, 1:2], data = sam.data[, 3],
                             theo = TRUE, my.l = my.l,
                             envelop    = input$show_envelop,
                             show_nbins = input$show_nbins,
                             bins = input$nbins,
                             maxdist = as.numeric(input$maxdist),
                             xlab = "time")
        } else {
          simu.data <- data.frame(coords, Y = Y[, i])
          sam.data <- if (input$sampling == "random")
            simu.data %>% sample_n(input$sample)
          else
            simu.data[geo_sample(xy.all = as.matrix(simu.data[, 1:2]),
                                 n = input$sample, delta = input$delta), ]

          m1[[i]] <- ggplot(simu.data, aes(x = x, y = y, fill = Y)) +
            geom_raster() +
            scale_fill_gradientn(colours = col) +
            geom_point(data = sam.data, aes(x = x, y = y, col = Y), size = 2) +
            coord_fixed() +
            theme_bw() +
            theme(legend.title = element_blank())

          p1[[i]] <- ggvario(coords = sam.data[, 1:2], data = sam.data[, 3],
                             theo = TRUE, my.l = my.l,
                             envelop    = input$show_envelop,
                             show_nbins = input$show_nbins,
                             bins = input$nbins,
                             maxdist = as.numeric(input$maxdist))
        }
      }
    })

    list(m1 = m1, p1 = p1)
  })

  output$map2 <- renderPlot({
    sim.results()$m1[[as.numeric(input$myslider)]]
  })

  current_vario <- reactive({
    sim.results()$p1[[as.numeric(input$myslider)]]
  })

  output$plot2 <- renderPlot({
    current_vario()
  })

  output$formula <- renderUI({
    cor.expr(input$covmodel)
  })

  output$downloadPlot <- downloadHandler(
    filename = function() paste0("variogram_", Sys.Date(), ".png"),
    content  = function(file) {
      ggsave(file, plot = current_vario(), width = 8, height = 5, dpi = 150)
    }
  )
}

shinyApp(ui = ui, server = server)
