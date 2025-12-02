library(dplyr, help, pos = 2, lib.loc = NULL)

# Mapping to formal names
sample_name_mapping <- list(
  brain = "Brain",
  thyroid = "Thyroid",
  thymus = "Thymus",
  heart = "Heart",
  lung = "Lung",
  liver = "Liver",
  spleen = "Spleen",
  pancreas = "Pancreas",
  stomach = "Stomach",
  small = "Small Intestine",
  large = "Large Intestine",
  adrenal = "Adrenal Gland",
  kidney = "Kidney",
  muscle = "Muscle",
  adipose = "Adipose Tissue",
  bone = "Bone Marrow",
  testis = "Testis",
  sperm = "Sperm"
)

data_column_mapping <- list(
  polya_read_num = "Poly(A) Read Number",
  median_polya = "Median Poly(A) Length",
  mean_polya = "Mean Poly(A) Length",
  bam_count = "BAM Count",
  bam_cpm = "BAM CPM"
)

prefix_mapping <- list(
  polya_read_num = "fulllength_",
  median_polya = "fulllength_",
  mean_polya = "fulllength_",
  bam_count = "", 
  bam_cpm = ""
)

batchGeneQueryUI <- function(id) {
  ns <- NS(id)
  tagList(

    tags$head(
      tags$style(HTML("
        .dataTables_wrapper {
          overflow-x: auto !important;
          width: 98%;
          max-width: 98%; 
        }
        
        .dataTables_scrollBody {
          overflow-y: hidden !important;
        }
      "))
    ),

    sidebarLayout(
      sidebarPanel(
        class = "custom-sidebar",  # Add this line to apply the custom CSS class
        textAreaInput(ns("batch_gene"), 
        "Enter Batch Gene Names:", 
        value = "Actb,\nENSMUST00000092163\nENSMUSG00000063457", 
        placeholder = "Enter genes separated by comma or newline, e.g.\nActb,\nENSMUST00000092163\nENSMUSG00000063457", 
        width='100%', 
        height='300px'
        ),

        # tags$hr(),
        tags$div(style = "font-weight: bold;", "Samples:"),
        tags$div(
          style = "column-count: 2; max-height: 300px; overflow-y: auto;",
          checkboxGroupInput(
            ns("selected_samples"),
            label = NULL,
            choices = setNames(names(sample_name_mapping), unlist(sample_name_mapping)),
            selected = names(sample_name_mapping) 
          )
        ),

        tags$hr(),
        tags$div(style = "font-weight: bold;", "Data Columns:"),
        checkboxGroupInput(
          ns("selected_types"),
          label = NULL,
          choices = setNames(names(data_column_mapping), unlist(data_column_mapping)),
          selected = names(data_column_mapping) 
        ),
        actionButton(ns("submit_batch"), "Submit", class = "btn-primary"),

        width = 2
      ),
      mainPanel(
        fluidRow(
            column(
            width = 3,
            h3("Batch Gene Query Results")
            ),
            column(
            width = 7,
            div(style = "text-align: left; margin-top: 25px;",
                downloadButton(ns("download_data"), "DOWNLOAD CSV", class = "btn-default")
            )
            )
        ),
        DT::dataTableOutput(ns("batch_query_result")),
        width = 10
     )
    )
  )
}

batchGeneQueryServer <- function(id) {
  moduleServer(id, function(input, output, session) {
    # Load necessary packages
    library(DBI)
    library(pool)
    library(DT)
    library(writexl)
    
    pool <- dbPool(
      RSQLite::SQLite(),
      dbname = "data/gene_cpm_polya.sqlite"
    )
    
    onStop(function() {
      poolClose(pool)
    })
    
    batch_query <- eventReactive(input$submit_batch, {
      req(input$batch_gene, input$selected_samples, input$selected_types)
      
      gene_list <- unlist(strsplit(input$batch_gene, "[,\n]")) %>%
        trimws() %>%
        .[nchar(.) > 0]
      
      if (length(gene_list) == 0) return(data.frame())

      dynamic_cols <- unlist(lapply(input$selected_samples, function(sample) {
        lapply(input$selected_types, function(type) {
          prefix <- prefix_mapping[[type]]
          paste0(
            sample, 
            "_", 
            prefix, 
            type
          )
        })
      }))
      
      base_cols <- c("isoform_id", "gene_id", "gene_symbol")
      all_cols <- c(base_cols, dynamic_cols)

      if (length(all_cols) == 3) return(data.frame())

      query <- sprintf(
        "SELECT %s FROM gene_table 
        WHERE 
          isoform_id_no_version IN ('%s') OR
          gene_id_no_version IN ('%s') OR
          LOWER(gene_symbol) IN ('%s')",
        paste(all_cols, collapse = ", "),
        paste(sub("\\..*", "", gene_list), collapse = "','"),
        paste(sub("\\..*", "", gene_list), collapse = "','"),
        paste(tolower(gene_list), collapse = "','")
      )
      
      result <- dbGetQuery(pool, query)
      result
    })

    output$download_data <- downloadHandler(
        filename = function() {
            paste0("gene_query_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
        },
        content = function(file) {
            req(batch_query())
            data <- batch_query()
            
            if (nrow(data) == 0) {
            showNotification("No data for download", type = "warning")
            return()
            }
            
            readr::write_csv(data, file)
        }
    )
    
    output$batch_query_result <- DT::renderDT({
      req(batch_query())
      datatable(
        batch_query(),
        # extensions = 'Scroller',
        extensions = 'FixedColumns',
        options = list(
          deferRender = TRUE,
          scrollX = TRUE,
          scrollY = FALSE,
          paging = TRUE,
          pageLength = 10, 
          lengthMenu = c(10, 25, 50),
          fixedColumns = list(leftColumns = 2),
          columnDefs = list(
            list(width = '200px', targets = 0:2)
          )
        )
      ) %>% 
        formatStyle(columns = 1:3, `white-space` = 'nowrap')
    })
  })
}
