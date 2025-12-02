# modules/igv_module.R
library(shiny)
library(igvShiny)
library(shinycssloaders)

igvModuleUI <- function(id) {
  ns <- NS(id)
  sidebarLayout(
    sidebarPanel(
      tagList(
        fluidRow(
          column(8, selectInput(ns("species"), NULL, choices = "Mouse", selected = "Mouse")),
          # column(4, actionButton(ns("refreshGenome"), "Refresh", class = "btn-primary"))
        ),
        fluidRow(
          column(8, textInput(ns("roi"), NULL, value = "chr5:142,902,115-142,907,754", placeholder = "e.g chr5:142,902,115-142,907,754")),
          column(4, actionButton(ns("searchButton"), "Search"))
        ),
        fluidRow(
          column(12, sliderInput(ns("subsample_percentage"), "Subsample Percentage:", min = 10, max = 100, value = 100, step = 10, width = '80%'))
        ),
        uiOutput(ns("sample_ui")),
        tags$div(style = "font-weight: bold;", "Params:"),
        checkboxInput(ns("filter_mapq"), "Filter out MAPQ<60", value = FALSE),
      ),
      width = 3
    ),
    mainPanel(
      shinycssloaders::withSpinner(
        igvShinyOutput(ns("igvShiny"), height = "800px", width = '1300px'),
        type = 4, color = "#0dc5c1"
      ),
    )
  )
}

igvModuleServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    output$igvShiny <- renderIgvShiny({
      genomeOptions <- igvShiny::parseAndValidateGenomeSpec(
        genomeName = "mm10",
        initialLocus = "chr5:142,902,115-142,907,754"
      )
      igvShiny(genomeOptions)
    })
    
    output$sample_ui <- renderUI({
      tagList(
        checkboxGroupInput(ns("samples"), "Add Samples:", 
                          choices = c(
                            "Brain",
                            "Thyroid",
                            "Thymus",
                            "Heart",
                            "Lung",
                            "Liver",
                            "Spleen",
                            "Pancreas",
                            "Stomach",
                            "Small Intestine",
                            "Large Intestine",
                            "Adrenal Gland",
                            "Kidney",
                            "Muscle",
                            "Adipose Tissue",
                            "Bone Marrow",
                            "Testis",
                            "Sperm"),
                          ),
        # helpText("Note: Large files will be subsampled (20%)")
      )
    })
    
    loaded_bams <- reactiveValues(loaded = character(0))
    
    observe({
      current_selected <- input$samples
      previous_selected <- loaded_bams$loaded
      
      to_remove <- setdiff(previous_selected, current_selected)
      purrr::walk(to_remove, ~{
        track_name <- paste(.x, paste0("subsample_", input$subsample_percentage, "%"), sep = "_")
        try(igvShiny::removeTracksByName(session, ns("igvShiny"), track_name))
      })
      
      to_add <- setdiff(current_selected, previous_selected)
      if(length(to_add) > 0){
        showNotification("Loading samples...(It may take some time)", id = "loading", type = "message")
        
        filter_mapq <- input$filter_mapq
        subsample_percentage <- input$subsample_percentage / 100
        
        promises <- purrr::map(to_add, function(sample_name){
          future({
            set.seed(123)
            bam_path <- bam_files$Mouse[[sample_name]]
            
            param <- Rsamtools::ScanBamParam(
              what = c("mapq"),
              tag = c("pa", "gi"),
              mapqFilter = if (filter_mapq) 60 else NA_integer_
            )
            
            raw_data <- GenomicAlignments::readGAlignments(bam_path, param = param)
            if(length(raw_data) > 1000) raw_data <- raw_data[sample(length(raw_data), subsample_percentage*length(raw_data))]
            
            list(name = sample_name, data = raw_data)
          }, seed = TRUE)
        })
        
        promises_all <- promises::promise_all(.list = promises)
        promises_all %...>% {
          removeNotification("loading")
          purrr::walk(., ~{
            if(!is.null(.x$data)){
              igvShiny::loadBamTrackFromLocalData(
                session = session,
                id = ns("igvShiny"),
                trackName = paste(.x$name, paste0("subsample_", input$subsample_percentage, "%"), sep = "_"),
                data = .x$data,
                displayMode = "FULL",
              )
            }
          })
        } %...!% (function(err){
          showNotification(paste("Error:", err$message), type = "error")
        })
      }
      loaded_bams$loaded <- current_selected
    })
    
    observeEvent(input$searchButton, {
      req(input$roi)
      showGenomicRegion(session, ns("igvShiny"), input$roi)
    })
  }
)}