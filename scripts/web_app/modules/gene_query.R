# Gene Query Module
source("modules/plot_functions.R")  # Source the new plotting functions

geneQueryUI <- function(id) {
  ns <- NS(id)
  tagList(
    sidebarLayout(
      sidebarPanel(
        textInput(ns("gene_name"), "Enter Gene:", value = "Rps15", placeholder = "e.g. Rps15"),
        actionButton(ns("submit_gene"), "Submit"),
        
        uiOutput(ns("transcript_selector")), 
        tags$hr(),
        tags$div(style = "font-weight: bold;", "Params:"),
        sliderInput(ns("bin_length"), "Bin Length (nt):",
                    min = 1, max = 20, value = 5, step = 1),
        sliderInput(ns("vmax"), "Maximum Value:",
                    min = 0.02, max = 0.2, value = 0.08, step = 0.01),
        selectInput(ns("cmap"), "Color Map:",
                    choices = c("Reds", "Blues", "Greens", 
                              "Purples", "Oranges", "Greys",
                              "YlOrRd", "YlGnBu", "RdPu", "RdBu", 
                              'PuOr', "RdYlBu"),
                    selected = "Reds"),
        checkboxInput(ns("revert_cmap"), "Revert Colormap", value=FALSE),
        tags$hr(),
        sliderInput(
          ns("min_sample_num"), 
          "min sample_num: ", 
          min = 1, max = 18, value = 16, step = 1
        ),
        width = 2,
      ),
      mainPanel(
        h3("Gene Query Results"),
        textOutput(ns("result")),
        uiOutput(ns("plots_container")),
        uiOutput(ns("distance_tables")),
        width = 10
      )

    )
)}

geneQueryServer <- function(id, data_info, data_array, distance_matrix, gene_list, isoform2symbol, isoform2gene, gene_merged_info) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    query_dis_per_gene <- function(
      transcript_id, distance_matrix, gene_list,
      isoform2symbol, isoform2gene, gene_merged_info
    ){
      idx <- match(transcript_id, gene_list)
      if (is.na(idx)) return(NULL)

      tmp <- data.frame(
        distance  = as.numeric(distance_matrix[idx, ]),
        short_gene_info = gene_list,
        stringsAsFactors = FALSE
      )
      tmp$gene_symbol <- isoform2symbol[tmp$short_gene_info]
      tmp$gene_id     <- isoform2gene[tmp$short_gene_info]

      tmp <- dplyr::left_join(
              tmp,
              gene_merged_info[, c("short_gene_info", "sample_info")],
              by = "short_gene_info"
            ) |>
            dplyr::mutate(sample_num = lengths(sample_info))

      tmp[order(tmp$distance), ]
    }

    submitted_gene <- reactiveVal(NULL)
    submitted_params <- reactiveVal(NULL)
    all_transcripts <- reactiveVal(NULL)

    # Input type
    input_type <- reactive({
      gene <- input$gene_name
      if (is.null(gene) || gene == "") {
        return("gene_symbol") 
      }

      ifelse(
        startsWith(gene, "ENSMUST"), "transcript_id",
        ifelse(
          startsWith(gene, "ENSMUSG"), "gene_id",
          "gene_symbol"
        )
      )
    })

    observeEvent(input$submit_gene, {
      submitted_gene(input$gene_name)  
      submitted_params(list(           
        bin_length = input$bin_length,
        vmax = input$vmax,
        cmap = input$cmap,
        revert_cmap = input$revert_cmap,
        selected_transcripts = input$selected_transcripts
      ))
    })

    observeEvent(c(input$bin_length, input$vmax, input$cmap, input$selected_transcripts, input$revert_cmap), {
      submitted_params(list( 
        bin_length = input$bin_length,
        vmax = input$vmax,
        cmap = input$cmap,
        revert_cmap = input$revert_cmap,
        selected_transcripts = input$selected_transcripts
      ))
    })

    library(pool)

    pool <- dbPool(
      drv = RSQLite::SQLite(),
      dbname = "data/polyann_info.sqlite"
    )
    
    onStop(function() {
      poolClose(pool)
    })

    gene_data <- reactive({
      req(submitted_gene())
      
      query <- sprintf(
        "SELECT * FROM polyann_info 
        WHERE 
          isoform_id_no_version IN ('%s') OR
          gene_id_no_version IN ('%s') OR
          LOWER(gene_symbol) IN ('%s')",
        sub("\\..*", "", submitted_gene()),
        sub("\\..*", "", submitted_gene()),
        tolower(submitted_gene())
      )
      
      dbGetQuery(pool, query)
    })

    observeEvent(input$submit_gene, {
      req(gene_data())
      all_transcripts(unique(sub("\\..*", "", gene_data()$isoform_id)))
    })

    output$transcript_selector <- renderUI({
      req(all_transcripts())

      if(input_type() %in% c("gene_id", "gene_symbol")) {
        tagList(
          tags$hr(),
          checkboxGroupInput(
            ns("selected_transcripts"), "Select Transcripts:",
            choices = all_transcripts(),
            selected = ""
          )
        )
        
      }
    })

    colorbar_file <- reactive({
      req(submitted_params())
      generate_colorbar(
        vmax = submitted_params()$vmax,
        cmap = submitted_params()$cmap,
        revert = submitted_params()$revert_cmap
      )
    }) %>% debounce(500)

    plot_ids <- reactive({
      req(submitted_gene())


      main_id <- submitted_gene()
      if (input_type() %in% c("gene_id", "transcript_id")) {
        main_id <- sub("\\..*", "", submitted_gene())
      }

      valid_transcripts <- character(0) 
      if (input_type() %in% c("gene_id", "gene_symbol")) {
        valid_transcripts <- submitted_params()$selected_transcripts %>% 
                            na.omit() %>% 
                            as.character()
      }

      if (input_type() == "gene_symbol") {
        valid_transcripts <- unique(c(main_id, valid_transcripts))
      } else if (input_type() %in% c("gene_id", "transcript_id")) {
  
        valid_transcripts <- unique(c(main_id, valid_transcripts))
        valid_transcripts <- sub("\\..*", "", valid_transcripts)
      }

      ids <- unique(na.omit(c(main_id, valid_transcripts))) 

      if (length(ids) == 0)
        stop("Error: No valid plot IDs found!")

      ids
    })

    transcripts_for_table <- reactive({
      req(submitted_gene())

      if (input_type() == "transcript_id") {    
        return(sub("\\..*", "", submitted_gene()))
      }

      tr <- submitted_params()$selected_transcripts
      tr <- tr[!is.na(tr)]
      tr <- as.character(tr)

      if (length(tr) == 0) return(NULL)
      tr
    })


    plot_files <- reactive({
      req(plot_ids(), submitted_params())
      
      map_chr(plot_ids(), ~{
        paste0(
          "user/", .x, "_",
          submitted_params()$bin_length, "nt_",
          submitted_params()$vmax, "vmax_",
          ".png"
        )
      })
    })

    output$result <- renderText({
      req(input$submit_gene)
      
      if (nrow(gene_data()) == 0) {
        "No records found for this gene."
      } else {
        paste("Found", nrow(gene_data()), "records for gene:", input$gene_name)
      }
    })

    observeEvent(plot_files(), {
      req(gene_data(), plot_ids(), plot_files())
      
      walk2(plot_ids(), plot_files(), function(id, file) {
        if (is.na(id) || id == "") {
          stop("Error: plot_id is NA or empty!")
        }

        colname <- ifelse(
          id %in% all_transcripts(),
          "isoform_id_no_version",
          ifelse(
            input_type() == "transcript_id", 
            "isoform_id_no_version",
            ifelse(
               input_type() == "gene_id", 
               "gene_id_no_version",
                input_type()
            )

          )
        )

        if (is.na(colname) || colname == "") {
          stop("Error: colname is NA or empty!")
        }

        colname <- as.character(colname)

        # print(paste("Plotting for:", id, "with colname:", colname))

        p <- plot_mRNA_tissues(
          mRNA = id,
          data_info = gene_data(),
          data_array = data_array,
          colname = colname,
          bin_length = submitted_params()$bin_length,
          vmax = submitted_params()$vmax,
          cmap = submitted_params()$cmap,
          color_revert = submitted_params()$revert_cmap
        )

        save_plot_as_png(p, file, width = 700, height = 800)
      })
      
      existing_files <- list.files("user", full.names = TRUE)
      if(length(existing_files) > 200) {
        file.remove(head(sort(existing_files), -50))
      }
    })
    
    output$plots_container <- renderUI({
      req(plot_files())
      
      plot_output_list <- map(plot_files(), ~{
        image_id <- ns(digest::digest(.x))
        div(
          style = "margin-right: 10px;",
          imageOutput(
            image_id,
            height = "400px",
            width = "350px"
          )
        )
      })

      colorbar_output <- div(
        style = "margin-left: 5px; display: flex; flex-direction: column; justify-content: center;",
        imageOutput(
          ns("colorbar"),
        ), 
        div(
          style = "text-align: center; margin-top: 10px;",
        )
      )

      div(
        style = "display: flex; flex-wrap: wrap; gap: 10px;", 
        plot_output_list,
        colorbar_output
      )
    })

    output$colorbar <- renderImage({
      list(
        src = colorbar_file(),
        contentType = "image/png",
        width = "80px",
        height = "320px"
      )
    }, deleteFile = TRUE)
    
    observeEvent(plot_files(), {
      req(plot_files())
      
      walk(plot_files(), function(file) {
        image_id <- digest::digest(file)
        output[[image_id]] <- renderImage({
          list(
            src = file,
            contentType = "image/png",
            width = 350,
            height = 400
          )
        }, deleteFile = FALSE)
      })
    })
    
    output$result <- renderText({
      req(submitted_gene())
      
      if(nrow(gene_data()) == 0) {
        "No matching records found in database."
      } else {
        sprintf("Showing %d plots for: %s", 
               length(plot_ids()), 
               submitted_gene())
      }
    })

    distance_tables <- reactive({
      req(input$min_sample_num)
      trs <- transcripts_for_table()
      if (is.null(trs)) return(list())

      lapply(trs, function(tid){
        df <- query_dis_per_gene(
          tid, distance_matrix, gene_list,
          isoform2symbol, isoform2gene, gene_merged_info
        )
        if (is.null(df)) return(NULL)

        df <- df[df$sample_num >= input$min_sample_num, ]
        head(df[order(df$distance), ], 100)
      })
    })



    output$distance_tables <- renderUI({
      tbls <- distance_tables() 
      if (length(tbls) == 0) return(NULL) 

      tagList(
        tags$hr(),
        h3("Nearest-distance genes"),
        Map(function(tbl, i) {
          tid    <- transcripts_for_table()[i]
          tbl_id <- paste0("dt_", i)

          local({
            this_tbl <- tbl
            output[[tbl_id]] <- DT::renderDT({
              DT::datatable(
                this_tbl,
                options = list(pageLength = 20,
                              lengthMenu = c(20,50,100),
                              scrollX = TRUE),
                rownames = FALSE
              )
            })
          })

          tagList(
            h4(sprintf("Transcript: %s", tid)),
            DTOutput(ns(tbl_id)),
            tags$hr(),
          )
        }, tbls, seq_along(tbls))
      )
    })


  })
}
