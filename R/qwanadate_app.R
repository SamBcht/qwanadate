
#library(tidyverse)
library(dplyr)
library(tidyr)
library(readr)
library(purrr)
library(tibble)
library(stringr)

library(shiny)
library(bslib)
library(shinyFiles)
library(DT)
library(plotly)


# --------------------------
# Crossdating function
# --------------------------
crossdate_profile <- function(anat_series, dendro_data, min_overlap = 5, expected_start = NULL, lag_window = 20) {
  x <- anat_series$ringwidth_mm
  y <- dendro_data$TRW
  dendro_years <- dendro_data$year
  first_dendro_year <- dendro_years[1]
  last_dendro_year <- last(dendro_years)
  n_x <- length(x)
  n_y <- length(y)

  # Determine lag range
  if (!is.null(expected_start)) {
    lags <- (expected_start - first_dendro_year - lag_window):(expected_start - first_dendro_year + lag_window)
  } else {
    lags <- seq(-min_overlap, n_y - min_overlap)
  }

  # ---- SHORT SERIES HANDLING ----
  if (n_x < 3 && n_x >= 1) {

    results <- purrr::map_dfr(lags, function(lag) {
      idx_x <- seq_len(n_x)
      idx_y <- idx_x + lag
      valid <- idx_y >= 1 & idx_y <= n_y
      overlap <- sum(valid)

      if (overlap == 0) {
        return(tibble(
          lag = lag,
          r_val = NA_real_,
          t_val = NA_real_,
          p_val = NA_real_,
          overlap = 0,
          first_year = NA_real_,
          last_year = NA_real_,
          mean_diff = NA_real_,
          rmse = NA_real_
        ))
      }

      shifted_x <- x[valid]
      shifted_y <- y[idx_y[valid]]

      tibble(
        lag = lag,
        r_val = NA_real_,   # no correlation for short series
        t_val = NA_real_,
        p_val = NA_real_,
        overlap = overlap,
        first_year = dendro_years[min(idx_y[valid])],
        last_year = dendro_years[max(idx_y[valid])],
        mean_diff = mean(abs(shifted_x - shifted_y), na.rm = TRUE),
        rmse = sqrt(mean((shifted_x - shifted_y)^2, na.rm = TRUE)),
        adj_rmse = rmse / overlap,
        log_rmse = log1p(rmse),
        rxt = NA_real_,
        w_rmse = rmse * mean_diff
      )
    })

    # Keep valid temporal window
    results <- results %>%
      filter(first_year > first_dendro_year - 1,
             first_year + n_x < last_dendro_year + 3) %>%
      filter(lag >= 0) %>%
      mutate(
        rank_rmse = rank(rmse, na.last = "keep"),
        score = rank_rmse   # key change: RMSE drives selection
      )

    best <- results %>%
      filter(!is.na(score)) %>%
      slice_min(score, n = 1, with_ties = FALSE)

    return(list(
      profile = results,
      best_lag = best$lag,
      best_r = NA_real_,
      rank_r = NA_real_,
      best_t = NA_real_,
      rank_t = NA_real_,
      best_first_year = best$first_year,
      best_last_year = best$last_year,
      best_rmse = best$rmse,
      rank_rmse = best$rank_rmse
    ))
  }

  results <- purrr::map_dfr(lags, function(lag) {
    idx_x <- seq_len(n_x)
    idx_y <- idx_x + lag
    valid <- idx_y >= 1 & idx_y <= n_y
    overlap <- sum(valid)

    if (overlap < min_overlap) {
      tibble(
        lag = lag,
        r_val = NA_real_,
        t_val = NA_real_,
        p_val = NA_real_,
        overlap = overlap,
        first_year = NA_real_,
        last_year = NA_real_,
        mean_diff = NA_real_,
        rmse = NA_real_
      )
    } else {
      shifted_x <- x[valid]
      shifted_y <- y[idx_y[valid]]
      cor_test <- suppressWarnings(cor.test(shifted_x, shifted_y))

      tibble(
        lag = lag,
        r_val = unname(cor_test$estimate),
        t_val = unname(cor_test$statistic),
        p_val = cor_test$p.value,
        overlap = overlap,
        first_year = dendro_years[min(idx_y[valid])],
        last_year = dendro_years[max(idx_y[valid])],
        mean_diff = mean(abs(shifted_x - shifted_y), na.rm = TRUE),
        rmse = sqrt(mean((shifted_x - shifted_y)^2, na.rm = TRUE)),
        adj_rmse = rmse / overlap,
        log_rmse = log1p(rmse),
        rxt = abs(r_val)*t_val,
        w_rmse = rmse*mean_diff
      )
    }
  })

  results <- results %>%
    filter(first_year > first_dendro_year-1, first_year + n_x < last_dendro_year+3) %>%
    filter(lag >= 0) %>%
    mutate(
      rank_r = rank(-r_val),
      rank_t = rank(-t_val),
      rank_rmse = rank(rmse),
      score = rank(-r_val))  # simple ranking for best correlation

  best <- results %>%
    filter(!is.na(score)) %>%
    slice_min(score, n = 1, with_ties = FALSE)

  list(
    profile = results,
    best_lag = best$lag,
    best_r = best$r_val,
    rank_r = best$rank_r,
    best_t = best$t_val,
    rank_t = best$rank_t,
    best_first_year = best$first_year,
    best_last_year = best$last_year,
    best_rmse = best$rmse,
    rank_rmse = best$rank_rmse
  )
}

read_csv_auto <- function(file) {
  # Read first line to detect separator
  first_line <- readLines(file, n = 1)

  if (str_detect(first_line, ";")) {
    read_delim(file, delim = ";", escape_double = FALSE, trim_ws = TRUE)
  } else {
    read_csv(file, show_col_types = FALSE)
  }
}

assign_years <- function(df_img, dendro_data, lag) {
  df_img %>%
    arrange(row_number()) %>%
    mutate(
      year = dendro_data$year[row_number() + lag]
    ) %>%
    filter(!is.na(year))
}

# =========================
# UI
# =========================
ui <- page_sidebar(

  theme = bs_theme(bootswatch = "solar"),

  title = span(
    #img(src = "Logo1_QWAnatools.png", height = 50),
    "QWAnatools – Crossdating Prototype",
    windowTitle = "QWAnatools"
  ),

  sidebar = sidebar(
    accordion(
      accordion_panel(
        "Inputs",

        textInput(
          "treeID",
          "Tree ID:",
          placeholder = "e.g. L20_F25"
        ),

        shinyDirButton(
          "base_dir",
          "Select base directory",
          "Choose a folder"
        ),
        helpText("Select the folder containing qwanamiz '_outputs' directories"),

        fileInput(
          "dendro_file",
          "Select reference ringwidth",
          accept = c(".csv", ".rwl")
        ),
        helpText("CSV must contain columns 'year' and 'TRW'."),

        selectInput(
          inputId = "trw_scaling",
          label = "Dendro TRW unit",
          choices = c(
            "mm / 100" = 1,
            "mm / 10"  = 10,
            "mm"       = 100,
            "µm"       = 0.1
          ),
          selected = 1
        ),


        actionButton(
          "search_files",
          "Search files",
          class = "btn-primary",
          icon = icon("magnifying-glass")
        )
      ),

      accordion_panel(
        "Crossdating",

        numericInput(
          "lag_window",
          "Lag window (± years):",
          value = 5, min = 1, max = 30
        ),

        numericInput(
          "min_overlap",
          "Minimum overlap:",
          value = 3, min = 1, max = 20
        ),

        actionButton(
          "run_crossdating",
          "Run crossdating",
          class = "btn-danger",
          icon = icon("play")
        ),

        hr(),

        uiOutput("lag_ui")
      ),

      accordion_panel(
        "Dating validation",

        actionButton(
          "finalize_chron",
          "Create anatomical chronology",
          class = "btn-success",
          icon = icon("check")
        ),

        hr(),

        uiOutput("saving_ui")
      )

    )
  ),

  navset_underline(
    #title = "Results",
    #fillable = FALSE,
    #well = TRUE,

    #nav_panel("Status", verbatimTextOutput("status")),
    nav_panel("Input files",
              uiOutput("status"),
              hr(),
              DTOutput("file_table")),
    nav_panel("Reference TRW",
              uiOutput("treeID_override_ui"),
              DTOutput("dendro_table")),
    nav_panel(
      "Anatomical rings table",

      actionButton(
        "skip_rings",
        "Skip selected rings",
        icon = icon("minus-circle"),
        class = "btn-warning"
      ),

      actionButton(
        "reset_rings",
        "Reset skipped rings",
        icon = icon("rotate-left"),
        class = "btn-secondary"
      ),

      br(), br(),

      uiOutput("skip_message"),

      DTOutput("anat_table")
    ),

    # ---- Manual crossdating workspace ----
    nav_panel(
      "Manual crossdating",

      # ── TOP: aligned plot + summary table ──
      layout_columns(
        col_widths = c(9, 3),

        # ---- Aligned plot ----
        card(
          card_header("Aligned ringwidth series"),
          plotlyOutput("aligned_plot", height = "360px")
        ),

        # ---- Crossdating summary ----
        card(
          card_header("Crossdating summary"),
          DTOutput("crossdating_table")
        )
      ),

      hr(),

      # ── BOTTOM: profile inspection ──
      layout_columns(
        col_widths = c(6,6),
        card(
          card_header("Crossdating profile"),
          plotlyOutput("profile_plot", height = "260px")
        ),
        card(
          card_header("Crossdating table"),
          DTOutput("profile_table")
        )
      )
    ),

    nav_panel(
      "Validation",

      layout_columns(
        col_widths = c(6,6),
        card(
          card_header("Anatomical vs. dendro series"),
          plotlyOutput("final_chron_plot", height = "360px")
        ),
        card(
          card_header("Tables"),
          DTOutput("segplot_table"),
          hr(),
          DTOutput("missing_years_table")
        )
      ),
      hr(),

      # Show the dating decisions dataframe
      DTOutput("dating_decisions_table")


    )
  )
)

# =========================
# Server
# =========================
server <- function(input, output, session) {

  # ---- Drives for shinyFiles ----
  volumes <- shinyFiles::getVolumes()()
  shinyDirChoose(input, "base_dir", roots = volumes, session = session)

  base_dir <- reactive({
    req(input$base_dir)
    parseDirPath(volumes, input$base_dir)
  })

  # ---- Find "_outputs" directories and *_rings.csv files ----
  found_files <- eventReactive(input$search_files, {
    req(input$treeID, base_dir())

    output_dirs <- list.dirs(path = base_dir(), full.names = TRUE, recursive = TRUE) %>%
      keep(~ str_detect(basename(.x), paste0("^", input$treeID, ".*_outputs$")))

    if (length(output_dirs) == 0) {
      showNotification("No '_outputs' directories found for this Tree ID", type = "warning")
      return(NULL)
    }

    all_files <- map_df(output_dirs, function(dir) {
      files <- list.files(path = dir, pattern = "_rings\\.csv$", full.names = TRUE)
      if (length(files) == 0) return(NULL)
      tibble(
        output_dir = dir,
        imageID    = str_remove(basename(files), "_rings\\.csv$"),
        file_path  = files
      )
    })

    if (nrow(all_files) == 0) {
      showNotification("No *_rings.csv files found inside '_outputs' directories", type = "warning")
      return(NULL)
    }

    all_files
  })

  # ---- Status ----
  output$status <- renderUI({
    req(input$search_files)

    files <- found_files()

    # Check base directory
    if (is.null(base_dir())) {
      return(
        tags$div(
          class = "alert alert-warning",
          icon("exclamation-triangle"),
          "No base directory selected."
        )
      )
    }

    # Check if files were found
    if (is.null(files) || nrow(files) == 0) {
      return(
        tags$div(
          class = "alert alert-danger",
          icon("times-circle"),
          "No anatomical files found. Check Tree ID input"
        )
      )
    }

    # Check dendro data
    dendro_df <- dendro_raw()
    df_tree <- dendro_df %>% filter(Tree.ID == input$treeID)
    tree_warning <- NULL
    if (nrow(df_tree) > 0) {
      tree_warning <- tags$div(
        class = "alert alert-info",
        icon("check2-circle"),
        HTML(paste("Tree '", input$treeID, "' found in reference data.'"))
      )
    }
    if (nrow(df_tree) == 0) {
      tree_warning <- tags$div(
        class = "alert alert-warning",
        icon("exclamation-triangle"),
        HTML(paste("Tree ID '", input$treeID, "' not found in reference data. Please check Reference TRW panel and select the correct Tree ID.")
        ))
    }

    # Normal summary if all OK
    summary_box <- tags$div(
      class = "alert alert-success",
      icon("check-circle"),
      HTML(paste(
        "Tree ID: ", input$treeID, "<br>",
        "Base directory: ", base_dir(), "<br>",
        "'_outputs' folders found: ", length(unique(files$output_dir)), "<br>",
        "*_rings.csv files found: ", nrow(files)
      ))
    )

    # Return combined UI
    tagList(
      summary_box,
      tree_warning
    )
  })


  # ---- Table of found files ----
  output$file_table <- renderDT({
    req(found_files())

    datatable(
      found_files(),
      options = list(
        pageLength = 10,
        scrollX = TRUE
      )
    )
  })


  # ---- Reference dendro data ----
  dendro_raw <- reactive({
    req(input$dendro_file, input$treeID)

    ext <- tools::file_ext(input$dendro_file$name)

    if (tolower(ext) == "rwl") {

      # ---- RWL input ----
      rwl <- dplR::read.rwl(input$dendro_file$datapath)

      df <- rwl %>%
        as.data.frame() %>%
        rownames_to_column(var = "year") %>%
        mutate(year = as.integer(year)) %>%
        pivot_longer(
          cols = -year,
          names_to = "Tree.ID",
          values_to = "TRW"
        )

    } else {

      # ---- CSV input ----
      df <- read_delim(
        input$dendro_file$datapath,
        delim = ";",
        escape_double = FALSE,
        trim_ws = TRUE,
        col_types = cols(...1 = col_skip())
      )

      names(df) <- str_trim(names(df))

      validate(
        need("Tree.ID" %in% names(df), "Dendro CSV must contain 'Tree.ID' column"),
        need(all(c("year", "TRW") %in% names(df)),
             "Dendro CSV must contain 'year' and 'TRW'")
      )
    }

    df
  })

  dendro_try <- reactive({
    req(dendro_raw(), input$treeID)

    dendro_raw() %>%
      filter(Tree.ID == input$treeID) %>%
      arrange(year)
  })

  output$treeID_override_ui <- renderUI({
    req(input$search_files, dendro_raw(), input$treeID)

    if (nrow(dendro_try()) > 0) return(NULL)

    tagList(
      tags$div(
        class = "alert alert-warning",
        icon("exclamation-triangle"),
        " Tree ID not found in reference data"
      ),

      selectInput(
        "treeID_override",
        "Select matching Tree ID from file",
        choices = sort(unique(dendro_raw()$Tree.ID)),
        selected = NULL
      )
    )
  })

  final_tree_id <- reactive({
    if (nrow(dendro_try()) > 0) {
      input$treeID
    } else {
      req(input$treeID_override)
      input$treeID_override
    }
  })

  dendro_data <- reactive({
    req(final_tree_id(), input$trw_scaling)

    df_tree <- dendro_raw() %>%
      filter(Tree.ID == final_tree_id()) %>%
      mutate(TRW = TRW * as.numeric(input$trw_scaling)) %>%
      arrange(year)

    validate(
      need(nrow(df_tree) > 0,
           paste("No dendro data found for Tree ID:", final_tree_id()))
    )

    df_tree
  })

  # ---- Show reference TRW ----
  output$dendro_table <- renderDT({
    req(dendro_data())

    datatable(
      dendro_data(),
      options = list(
        pageLength = 10,
        scrollX = TRUE
      )
    )
  })

  # ---- Show aggregated anatomical files ----
  aggregated_anat <- reactive({
    req(found_files(), input$treeID)

    files <- found_files()

    df <- map_df(files$file_path, function(f) {

      df <- read_csv_auto(f)

      df %>%
        mutate(
          imageID = str_remove(basename(f), "_rings\\.csv$")
        ) %>%
        mutate(
          rw_diff = ringwidth - rw_from_cells,
          ringwidth_mm = rw_from_cells / 10,
          rwbis_mm = ifelse(
            abs(rw_diff) > 15,
            (earlywood_width + latewood_width) / 10,
            ringwidth / 10
          )
        ) %>%
        separate(
          imageID,
          into = c("treeID", "woodpieceID", "scanID"),
          sep = "-",
          remove = FALSE,
          fill = "right"
        ) %>%
        relocate(
          imageID,
          treeID,
          woodpieceID,
          scanID,
          year_x,
          ringwidth_mm,
          ringwidth,
          rw_from_cells,
          earlywood_width,
          latewood_width,
          rw_diff,
          .before = everything()
        )
    })

    # APPLY SKIPPING RINGS
    if (length(skipped_rows()) > 0) {
      df$ringwidth_mm[skipped_rows()] <- NA
    }

    df
  })

  qwa_rings <- reactive({
    req(aggregated_anat())

    aggregated_anat() %>%
      filter(!is.na(ringwidth_mm)) %>%
      arrange(imageID)
  })


  # ---- Show aggregated table ----
  output$anat_table <- renderDT({
    req(aggregated_anat())

    datatable(
      aggregated_anat(),
      selection = "multiple",
      options = list(
        pageLength = 25,
        scrollX = TRUE,
        scrollY = "60vh",
        dom = "tip",
        columnDefs = list(
          list(targets = which(names(aggregated_anat()) == "filePath"),
               visible = FALSE)
        )
      )
    )
  })

  skipped_rows <- reactiveVal(integer(0))

  observeEvent(input$skip_rings, {
    req(input$anat_table_rows_selected)

    current <- skipped_rows()
    new <- input$anat_table_rows_selected

    skipped_rows(unique(c(current, new)))
  })

  observeEvent(input$reset_rings, {
    skipped_rows(integer(0))
  })

  output$skip_message <- renderUI({
    req(aggregated_anat())

    skipped <- skipped_rows()

    if (length(skipped) == 0) {
      return(NULL)
    }

    df <- aggregated_anat()
    skipped_df <- df[skipped, ]

    # Build readable message
    msg <- skipped_df %>%
      group_by(imageID) %>%
      summarise(
        rings = paste(year_x, collapse = ", "),
        .groups = "drop"
      ) %>%
      mutate(text = paste0(imageID, ": ", rings)) %>%
      pull(text)

    tags$div(
      class = "alert alert-light",
      icon("exclamation-circle"),
      HTML(paste(
        "<b>Excluded from crossdating:</b><br>",
        paste(msg, collapse = "<br>")
      ))
    )
  })

  # CROSSDATING
  cross_results <- eventReactive(input$run_crossdating, {
    req(qwa_rings(), dendro_data())

    qwa_rings() %>%
      group_split(imageID) %>%
      map_df(function(df_img) {

        img_id <- unique(df_img$imageID)

        cross <- crossdate_profile(
          anat_series  = df_img,
          dendro_data  = dendro_data(),
          lag_window   = input$lag_window,
          min_overlap = input$min_overlap
        )

        tibble(
          imageID         = img_id,
          profile         = list(cross$profile),
          best_lag        = cross$best_lag,
          best_first_year = cross$best_first_year,
          best_last_year  = cross$best_last_year,
          best_r          = cross$best_r,
          rank_r          = cross$rank_r,
          best_t          = cross$best_t,
          rank_t          = cross$rank_t,
          best_rmse       = cross$best_rmse,
          rank_rmse       = cross$rank_rmse
        )
      })
  })

  selected_lags <- reactiveVal(NULL)

  observeEvent(cross_results(), {
    lags <- cross_results() %>%
      select(imageID, best_lag) %>%
      deframe()

    selected_lags(lags)
  })

  output$lag_ui <- renderUI({
    req(selected_lags())

    tagList(
      lapply(names(selected_lags()), function(img) {
        numericInput(
          inputId = paste0("lag_", img),
          label = paste("Lag for", img),
          value = selected_lags()[[img]],
          min = -50, max = 200, step = 1
        )
      })
    )
  })

  observe({
    req(selected_lags())

    for (img in names(selected_lags())) {
      local({
        img_local <- img
        observeEvent(input[[paste0("lag_", img_local)]], {
          lags <- selected_lags()
          lags[[img_local]] <- input[[paste0("lag_", img_local)]]
          selected_lags(lags)
        }, ignoreInit = TRUE)
      })
    }
  })

  dated_qwa <- reactive({
    req(qwa_rings(), dendro_data(), selected_lags())

    qwa_rings() %>%
      left_join(
        tibble(
          imageID = names(selected_lags()),
          lag     = unname(selected_lags())
        ),
        by = "imageID"
      ) %>%
      group_split(imageID) %>%
      map_df(~ assign_years(.x, dendro_data(), unique(.x$lag)))
  })

  output$aligned_plot <- renderPlotly({
    req(dated_qwa(), dendro_data())
    validate(need(nrow(dated_qwa()) > 0, "No aligned data to display"))

    #sel_img <- selected_image()

    plot_data <- bind_rows(
      dated_qwa() %>%
        select(year, ringwidth_mm, imageID) %>%
        rename(value = ringwidth_mm)%>%
        mutate(type = "anat"),

      dendro_data() %>%
        mutate(imageID = "Dendro",
               type = "dendro") %>%
        rename(value = TRW) %>%
        select(year, value, imageID, type)
    )

    plot_ly(
      data = plot_data,
      x = ~year,
      y = ~value,
      color = ~imageID,
      type = "scatter",
      mode = "lines",
      line = list(width = 1.5),
      hovertemplate = paste(
        "<b>%{legendgroup}</b><br>",
        "Year: %{x}<br>",
        "Ring width: %{y:.2f} mm<extra></extra>"
      )
    ) %>%
      layout(
        #title = paste("Aligned ringwidth series –", input$treeID),
        xaxis = list(title = "Year"),
        yaxis = list(title = "Ring width (mm)"),
        legend = list(
          orientation = "v",
          x = 1.02,
          y = 1,
          xanchor = "left",
          yanchor = "top"
        )
      )
  })





  output$crossdating_table <- renderDT({
    req(cross_results())

    table_data <- cross_results() %>%
      select(
        imageID,
        best_lag,
        best_first_year,
        best_last_year,
        best_r,
        rank_r,
        best_t,
        rank_t,
        best_rmse,
        rank_rmse
      )

    table_data %>%
      datatable(
        selection = "single",
        options = list(
          pageLength = 10,
          scrollX = TRUE,
          scrollY =TRUE,
          dom = "tip"
        ),
        rownames = FALSE
      )
  })

  selected_image <- reactive({
    req(input$crossdating_table_rows_selected)

    cross_results()$imageID[
      input$crossdating_table_rows_selected
    ]
  })

  selected_profile <- reactive({
    req(selected_image(), cross_results())

    cross_results() %>%
      filter(imageID == selected_image()) %>%
      pull(profile) %>%
      .[[1]]
  })

  output$profile_plot <- renderPlotly({
    req(selected_profile(), selected_image(), selected_lags())

    df <- selected_profile()
    img <- selected_image()
    current_lag <- selected_lags()[[img]]

    plot_ly(
      data = df,
      x = ~lag,
      y = ~r_val,
      type = "bar",
      source = "profile_click",
      marker = list(
        color = ~rmse,
        colorscale = "YlGnBu",
        reversescale = FALSE,
        colorbar = list(title = "RMSE")
      ),
      hovertemplate = paste(
        "Lag: %{x}<br>",
        "r: %{y:.3f}<br>",
        "RMSE: %{customdata[0]:.2f}<extra></extra>"
      ),
      customdata = ~rmse
    ) %>%
      add_trace(
        x = current_lag,
        y = df$r_val[df$lag == current_lag],
        type = "scatter",
        mode = "markers",
        marker = list(size = 12, color = "black", symbol = "diamond"),
        inherit = FALSE,
        showlegend = FALSE
      ) %>%
      layout(
        #title = paste("Crossdating profile –", selected_image()),
        xaxis = list(title = "Lag"),
        yaxis = list(
          title = "Correlation (r)",
          range = c(0, 1))
      )

    #event_register(p, "plotly_click")

  })

  output$profile_table <- renderDT({
    req(selected_profile(), selected_lags(), selected_image())

    current_lag <- selected_lags()[[selected_image()]]
    req(current_lag)

    selected_profile() %>%
      select(lag, r_val, t_val, rmse, rank_r, rank_rmse, overlap, first_year, last_year)  %>%
      datatable(
        selection = "single",
        options = list(pageLength = 10,
                       order = list(list(2, "desc")),
                       scrollX = TRUE,
                       scrollY =TRUE,
                       dom = "tip",
                       stripe = FALSE),
        rownames = FALSE
      ) %>%
      formatRound(columns=c('r_val', 't_val', 'rmse'), digits=3) %>%
      formatStyle(
        "lag",
        target = "row",
        border = styleEqual(
          current_lag,
          "4px solid #ffe6b3"
        ),
        fontWeight = styleEqual(
          current_lag,
          "bold"
        )
      ) %>%
      formatStyle(
        "lag",
        color = styleEqual(
          current_lag,
          "#ffe6b3"
        )
      )
  })

  observeEvent(
    plotly::event_data("plotly_click", source = "profile_click"),
    {
      ed <- plotly::event_data("plotly_click", source = "profile_click")
      req(ed, selected_image())

      lag_clicked <- ed$x
      img <- selected_image()

      lags <- selected_lags()
      lags[[img]] <- lag_clicked
      selected_lags(lags)
    }
  )

  final_chronology <- eventReactive(input$finalize_chron, {
    req(dated_qwa(), dendro_data(), selected_lags())

    anat_shifted <- dated_qwa() %>%
      select(imageID, year, ringwidth_mm)

    dendro_ref <- dendro_data() %>% select(year, TRW)

    # Choose closest ringwidth to dendro if multiple per year
    anat_final <- anat_shifted %>%
      left_join(dendro_ref, by = "year") %>%
      group_by(year) %>%
      slice_min(abs(ringwidth_mm - TRW), n = 1, with_ties = FALSE) %>%
      ungroup() %>%
      select(year, ringwidth_mm) %>%
      arrange(year)

    anat_final
  })

  missing_years <- reactive({
    req(final_chronology())

    yrs <- sort(final_chronology()$year)
    full_seq <- seq(min(yrs), max(yrs))

    tibble(
      year = setdiff(full_seq, yrs)
    )
  })

  output$missing_years_table <- renderDT({
    req(missing_years())

    datatable(
      missing_years(),
      rownames = FALSE,
      options = list(
        pageLength = 10,
        dom = "tip"
      )
    )
  })


  segplot_data <- reactive({
    req(dated_qwa())
    validate(need(nrow(dated_qwa()) > 0, "No dated QWA data"))

    dated_qwa() %>%
      group_by(imageID) %>%
      summarise(
        start_year = min(year, na.rm = TRUE),
        end_year   = max(year, na.rm = TRUE),
        n_rings    = n(),
        .groups = "drop"
      ) %>%
      arrange(start_year) %>%
      mutate(seg_row = row_number())
  })


  output$segplot_table <- renderDT({
    req(segplot_data())

    datatable(
      segplot_data(),
      rownames = FALSE,
      options = list(
        pageLength = 10,
        scrollX = TRUE,
        dom = "tip"
      )
    )
  })


  output$final_chron_plot <- renderPlotly({
    req(final_chronology(), dendro_data(), segplot_data())

    anat_final <- final_chronology()
    dendro_ref <- dendro_data()
    seg_df <- segplot_data() %>%
      mutate(y_seg = rep(c(1, 2, 3), length.out = n()))

    missing_years_vec <- missing_years()$year

    y_max <- max(
      c(
        anat_final$ringwidth_mm,
        dendro_ref$TRW
      ),
      na.rm = TRUE
    )

    # ---- Base chronology plot ----
    p <- plot_ly() %>%
      add_lines(
        data = anat_final,
        x = ~year,
        y = ~ringwidth_mm,
        name = "Anatomical",
        line = list(color = "darkseagreen"),
        hovertemplate =
          "Anatomical<br>Year: %{x}<br>RW: %{y:.2f} mm<extra></extra>"
      ) %>%
      add_lines(
        data = dendro_ref,
        x = ~year,
        y = ~TRW,
        name = "Dendro reference",
        line = list(color = "darkslategray"),
        hovertemplate =
          "Dendro<br>Year: %{x}<br>RW: %{y:.2f} mm<extra></extra>"
      )


    # ---- Add image segments on y2 ----
    for (i in seq_len(nrow(seg_df))) {
      p <- p %>%
        add_segments(
          x = seg_df$start_year[i],
          xend = seg_df$end_year[i],
          y = seg_df$y_seg[i],
          yend = seg_df$y_seg[i],
          yaxis = "y2",
          line = list(width = 6,
                      color = "lightseagreen"),
          hoverinfo = "text",
          text = paste0(
            "Image: ", seg_df$imageID[i], "<br>",
            "From: ", seg_df$start_year[i], "<br>",
            "To: ", seg_df$end_year[i], "<br>",
            "Rings: ", seg_df$n_rings[i]
          ),
          showlegend = FALSE
        )
    }

    # ---- Layout ----
    p %>%
      layout(
        xaxis = list(title = "Calendar year"),

        yaxis = list(
          title = "Ring width (mm)",
          rangemode = "tozero"
        ),

        shapes = lapply(missing_years_vec, function(x) {
          list(
            type = "line",
            x0 = x, x1 = x,
            y0 = 0, y1 = y_max,  # ymax is your top of y-axis
            line = list(dash = "dash", color = "#CB4B16")
          )
        }),

        yaxis2 = list(
          title = "",
          overlaying = "y",
          side = "right",
          range = c(0, y_max),
          showgrid = FALSE,
          showticklabels = FALSE
        ),

        legend = list(
          orientation = "h",
          x = 0.5,
          xanchor = "center",
          y = 1.05,
          yanchor = "bottom"
        )
      )
  })

  output$saving_ui <- renderUI({
    req(input$finalize_chron)

    tagList(
      h5("Save dating results"),

      shinyDirButton(
        "save_dir",
        "Choose output directory",
        "Select folder"
      ),

      textInput(
        "save_filename",
        "File name",
        value = paste0(input$treeID, "_dated.csv")
      ),

      actionButton(
        "save_chron",
        "Save dating file",
        class = "btn-primary",
        icon = icon("floppy-disk")
      )
    )
  })

  shinyDirChoose(
    input,
    "save_dir",
    roots = volumes,
    session = session
  )

  save_dir_path <- reactive({
    req(input$save_dir)
    parseDirPath(volumes, input$save_dir)
  })

  ring_structure <- reactive({
    req(aggregated_anat())

    aggregated_anat() %>%
      group_by(imageID) %>%
      summarise(
        n_rings_total = n(),
        n_rings_valid = sum(!is.na(ringwidth_mm)),

        first_valid_ring_x = min(year_x[!is.na(ringwidth_mm)]),
        last_valid_ring_x  = max(year_x[!is.na(ringwidth_mm)]),

        n_missing_start = first_valid_ring_x - min(year_x),
        n_missing_end   = max(year_x) - last_valid_ring_x,

        .groups = "drop"
      )
  })

  selected_cross_stats <- reactive({
    req(cross_results(), selected_lags())

    cross_results() %>%
      select(imageID, profile) %>%
      mutate(
        selected_lag = selected_lags()[imageID]
      ) %>%
      rowwise() %>%
      mutate(
        stats = list(
          profile %>% filter(lag == selected_lag)
        )
      ) %>%
      unnest(stats) %>%
      select(
        imageID,
        selected_lag,
        r_val,
        t_val,
        rmse,
        overlap,
        first_year,
        last_year
      ) %>%
      ungroup()
  })

  dating_decisions <- reactive({
    req(
      selected_cross_stats(),
      ring_structure(),
      dendro_data()
    )

    dendro_start <- min(dendro_data()$year)

    selected_cross_stats() %>%
      left_join(ring_structure(), by = "imageID") %>%
      mutate(
        first_calendar_year =
          dendro_start + selected_lag,

        last_calendar_year =
          first_calendar_year +
          (last_valid_ring_x - first_valid_ring_x),

        dating_method = "manual_crossdating"
      ) %>%
      relocate(
        imageID,
        dating_method,
        selected_lag,
        first_calendar_year,
        last_calendar_year
      )
  })

  output$dating_decisions_table <- renderDT({
    req(dating_decisions())

    datatable(
      dating_decisions(),
      rownames = FALSE,
      options = list(
        pageLength = 15,
        scrollX = TRUE,
        scrollY = "300px",
        dom = "tip",
        stripe = TRUE
      )
    )
  })


  observeEvent(input$save_chron, {
    req(dating_decisions(), save_dir_path(), input$save_filename)

    out_path <- file.path(save_dir_path(), input$save_filename)

    if (file.exists(out_path)) {
      showNotification(
        "File already exists — choose another name",
        type = "warning"
      )
      return()
    }

    write_csv(dating_decisions(), out_path)

    showModal(
      modalDialog(
        title = div(style = "color:#2c7a2c;", "Dating successfully saved"),

        div(
          style = "font-size:18px; text-align:center;",
          "All crossdated series have been saved.",
          br(), br(),
          "You can safely continue or export your results."
        ),

        easyClose = TRUE,
        footer = modalButton("OK"),
        size = "m"
      )
    )
  })



}

# =========================
# Run the app
# =========================
shinyApp(ui, server)
