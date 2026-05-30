# run_app.R
#
#' @title Run the variogramApp Shiny application
#' @description Launches the variogramApp interactive Shiny web application
#'   for exploring variograms and Gaussian random fields.
#' @details Opens the app in the default browser. Use Ctrl+C (or the Stop
#'   button in RStudio) to stop the app when finished.
#' @export
run_variog_app <- function() {
  appDir <- system.file('variogramApp', package='variogramApp')
  if (appDir == "") {
    stop("Could not find variogramApp. Try re-installing `variogramApp`.", call. = FALSE)
  }

  shiny::runApp(appDir, display.mode = "normal")
}
