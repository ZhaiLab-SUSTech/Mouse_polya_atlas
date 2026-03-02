library(shiny)
library(igvShiny)
library(GenomicAlignments)
library(Rsamtools)
library(future)
library(promises)
library(purrr)
library(shinyjs)
library(shinycssloaders)
library(RSQLite)
library(reshape2)
library(pheatmap)
library(grid)
library(gridExtra)
library(RColorBrewer)
library(ggplot2)
library(dplyr)
library(tidyr)
library(HDF5Array)  
library(DT)


source("config/genomes.R")
source("config/bam_files.R")
source("modules/gene_query.R")
source("modules/igv_module.R")
source("modules/help_module.R")
source("modules/batch_gene_query.R")
source("modules/homepage.R")

print('Loading polya data...')
data_array <- readRDS("data/list_data.rds")

message("Loading distance matrix ...")
h5_path <- "data/merge_dis_complete_bin5.h5"


distance_matrix <- HDF5Array(h5_path, "distance_matrix") 

info_df <- read.csv("data/info_df.csv", stringsAsFactors = FALSE)
gene_list <- unique(info_df$short_gene_info)


rev_df <- info_df[nrow(info_df):1, ]
rev_df <- rev_df[!duplicated(rev_df$short_gene_info), ]
isoform2symbol <- setNames(rev_df$gene_symbol, rev_df$short_gene_info)
isoform2gene <- setNames(rev_df$gene_id, rev_df$short_gene_info)
rev_df2 <- rev_df[!duplicated(rev_df$gene_id), ]
gene2symbol <- setNames(rev_df2$gene_symbol, rev_df2$gene_id)


gene_merged_info <- info_df |>
  dplyr::group_by(short_gene_info) |>
  dplyr::summarise(
      gene_symbol = first(gene_symbol),
      gene_id     = first(gene_id),
      sample_info = list(sample_info) 
  ) |>
  dplyr::ungroup()
head(gene_merged_info)

print('Data loaded!')


ui <- navbarPage(
  title = "Mouse Poly(A) Atlas",
  id = 'tabs',
  theme = "bootstrap",
  tags$head(
    tags$link(rel = "stylesheet", type = "text/css", href = "www/styles.css"),
    tags$style(HTML("
      .custom-sidebar {
        width: 105% !important;
      }
    "))
  ),
  

  tabPanel(
    "Home",
    homepageUI('home')
  ),
  
  # Gene Query
  tabPanel(
    "Gene Query",
    geneQueryUI("gene")
  ),
  
  # Batch Gene Query
  tabPanel(
    "Batch Gene Query",
    batchGeneQueryUI("batch_gene")
  ),
  
  # IGV
  # add feature: params: tags to be included
  tabPanel(
    "IGV",
    igvModuleUI("igv")
  ),

  # Help
  tabPanel(
    "Tutorial",
    helpModuleUI("help")
  ),

)

# Server
server <- function(input, output, session) {
  options(future.globals.maxSize = 1 * 1024^3)
  plan(multisession, workers = 3)

  observeEvent(input$tabs, {
    # shinyjs::runjs("Shiny.setInputValue('tab', $('#tabs .active').text())")
    if (input$tabs == "IGV") {
      showNotification("If IGV browser not loaded successfully, please refresh the page with `Ctrl/Command + Shift + R`!", type = "message")
    }
  })


  geneQueryServer(
    "gene",
    data_array       = data_array,
    distance_matrix  = distance_matrix,
    gene_list        = gene_list,
    isoform2symbol   = isoform2symbol,
    isoform2gene     = isoform2gene,
    gene_merged_info = gene_merged_info
  )
  igvModuleServer("igv")
  batchGeneQueryServer("batch_gene")


}

shinyApp(ui = ui, server = server)