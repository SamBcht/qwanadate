#' Launch qwanadate Shiny App
#'
#' @return
#' @export
#'
#' @examples
run_qwanadate <- function() {

  app_file <- system.file(
    "app",
    "qwanadate_app.R",
    package = "qwanadate"
  )

  if (app_file == "") {
    stop("Cannot find app script in installed package.")
  }

  e <- new.env(parent = globalenv())

  sys.source(app_file, envir = e)

  # bring ui/server into scope
  ui <- e$ui
  server <- e$server

  shiny::shinyApp(ui, server)
}
