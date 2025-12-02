# modules/help_module.R
library(shiny)
library(igvShiny)
library(shinycssloaders)

helpModuleUI <- function(id) {
  ns <- NS(id)
  fluidPage(
    # h1("Help"),
    tags$style(HTML("
        .custom-text {
          margin-left: 12%;
          margin-right: 12%;
          text-align: justify;
        }
        p {
          font-size: 16px; 
        }
        body {
          margin-bottom: 100px; 
        }
        figure {
          text-align: center;
          margin: 20px;
        }
        figcaption {
          font-style: bold;
          margin-top: 5px;
          text-align: center;
          color: #999999; 
        }
      ")),
    p("Welcome to the Mouse Poly(A) Tail Atlas Database, an online resource for exploring poly(A) tail length profiles across 18 different mouse organs under normal physiological conditions. Below is a guide on how to use this website effectively."),
    tags$hr(),

    # =================================
    
    h2("Gene Query"),
    
    h3("Input"),
    p("You can initiate a search by entering one of the following identifiers in the 'Enter Gene' field:"),
    tags$ul(
      tags$li("Gene Symbol (e.g., Rps15)"),
      tags$li("Gene ID (e.g., ENSMUSG00000064341)"),
      tags$li("Transcript ID (e.g., ENSMUST00000082908)")
    ),
    tags$img(), # enter genes
    p("Version numbers in Gene IDs and Transcript IDs are ignored during the search."),
    p("If the entered identifier matches multiple transcripts, a transcript selection panel will appear to allow detailed selection."),
    tags$img(), # transcript selection panel
    
    h3("Parameters"),
    p("Several parameters are available to customize the visualization and analysis:"),
    tags$ul(
      tags$li(tags$b("Bin Length (nt): "), "Size of bins for smoothing poly(A) tail length distribution. Default: 5 nt."),
      tags$li(tags$b("Maximum Value (vmax): "), "Upper limit of the color scale in the heatmap."),
      tags$li(tags$b("Color Map: "), "Color scheme options for heatmap visualization."),
      tags$li(tags$b("Revert Colormap: "), "Inverts the selected color map."),
      tags$li(tags$b("Min Sample Num: "), "Minimum number of tissues required for a gene to be included in nearest-distance analysis. Default: 16.")
    ),
    
    h3("Outputs"),
    p("Upon submitting a query, the following results will be generated:"),
    tags$ul(
      tags$li(tags$b("Heatmaps: "), "Poly(A) tail length profiles across tissues."),
      tags$li(tags$b("Nearest-Distance Gene Tables: "), "Top 100 genes with most similar poly(A) tail profiles to selected transcripts. The calculation of the distance can refer to the Method section of our article", tags$a(href="", "Mouse Pan-organ Full-length RNA Atlas Reveals Regulatory Gene Modules by Poly(A) Tail Patterns"), "."),
    ),
    div(style = "display: flex; justify-content: center; gap: 50px;",
        tags$figure(
            tags$img(
                src="helppage/heatmap.png", 
                style = "display: block; margin: 20px auto; width: 550px; height: 300px;"
            ),
            tags$figcaption("Heatmap")
        ),
        tags$figure(
            tags$img(
                src="helppage/distance_table.png", 
                style = "display: block; margin: 20px auto; width: 500px; height: 300px;"
            ),
            tags$figcaption("Nearest-Distance Gene Table")
        ),
    ),

    h3("Usage Tips"),
    tags$ul(
    tags$li("Please ensure that there are no spaces at the beginning or end of the input content."),
    tags$li("The input Gene symbol should have its first letter capitalized.")
    ),
    
    # h3("Example Walkthrough"),
    # tags$ol(
    #   tags$li("Enter 'Rps15' into the search box and click 'Submit'."),
    #   tags$li("If prompted, select one or more transcripts."),
    #   tags$li("Adjust parameters such as bin length or color map if needed."),
    #   tags$li("View generated heatmaps, and nearest-distance gene information.")
    # ),
    tags$hr(),

    # =================================

    h2("Batch Gene Query"),

    h3("Input"),
    p("The Batch Gene Query module allows users to retrieve poly(A) tail length and expression data for multiple genes at once."),
    p("You can enter multiple gene identifiers into the 'Enter Batch Gene Names' field. Accepted identifier types include:"),
    tags$ul(
    tags$li("Gene Symbol (e.g., Actb)"),
    tags$li("Gene ID (e.g., ENSMUSG00000063457)"),
    tags$li("Transcript ID (e.g., ENSMUST00000092163)")
    ),
    p("Entries can be separated by commas or newlines. Version numbers in Gene IDs and Transcript IDs are ignored during the search."),
    p("After inputting the gene list, select the tissue samples and data types you wish to retrieve, then click the 'Submit' button."),

    h3("Parameters"),
    p("Several customizable parameters are available:"),
    tags$ul(
    tags$li(tags$b("Samples: "), "Select one or more mouse organs/tissues for which poly(A) and expression data will be retrieved."),
    tags$li(tags$b("Data Columns: "), "Choose the type of information you want to retrieve, including median poly(A) length, mean poly(A) length, poly(A) read number, BAM count, and BAM CPM.")
    ),

    h3("Outputs"),
    p("After submission, a downloadable data table will be generated. The results can be downloaded as a csv file for offline analysis."),

    # h3("Example Walkthrough"),
    # tags$ol(
    # tags$li("Input multiple gene names (e.g., Actb, ENSMUST00000092163) into the batch input field."),
    # tags$li("Select the tissues and data types you want to retrieve."),
    # tags$li("Click 'Submit' to execute the batch query."),
    # tags$li("View and download the resulting data table.")
    # ),
    tags$hr(),

    # =================================

    h2("IGV Browser"),

    p("The IGV Browser module allows users to interactively visualize sequencing reads and poly(A) tail information across different mouse tissues, based on customized ", tags$a(href='https://github.com/igvteam/igv.js/',  "igv.js"), "project."), # hyper link

    h3("Input"),
    p("To explore a genomic region, enter a region of interest in the format 'chr:start-end' (e.g., chr5:142,902,115-142,907,754) into the search bar and click the 'Search' button. The genome is preloaded with the mm10 (GRCm38) mouse assembly."),
    p("Due to server load limitations, there may be occasions when the IGV browser does not load properly, especially during periods of high usage. If the expected browser window (see example below) does not appear, we recommend refreshing the page using:"),
    tags$ul(
        tags$li("For Windows/Linux users: ", tags$code("Ctrl + Shift + R")),
        tags$li("For macOS users: ", tags$code(HTML("Command(&#8984;) + Shift + R")))
    ),
    tags$img(src='helppage/image-5.png', style = "display: block; margin: 20px auto; width: 600px; height: 300px;"),
    p("Refreshing will typically resolve transient loading issues by forcing the browser to reload all resources."),

    h3("Parameters"),
    p("Several options are available to adjust the display:"),
    tags$ul(
    tags$li(tags$b("Samples: "), "Select one or more tissues to load corresponding sequencing reads as tracks. This may take some time. Your patience will be highly appreciated. In order to accelerate the loading process, the",  tags$code("seq"), "information has not been included, and all samples are, by default, pre-downsampled to approximately 180000 reads."),
    tags$li(tags$b("Subsample Percentage: "), "Specify the proportion of reads to load (from 10% to 100%) to optimize browser performance."),
    tags$li(tags$b("Filter out MAPQ<60: "), "Optionally exclude reads with mapping quality scores below 60 to enhance data reliability."),
    ),

    h3("Outputs"),
    p("The browser will display selected tissue tracks aligned to the specified genomic region. Each track represents reads mapped from a particular tissue, and subsampling settings apply to manage large data volumes."),

    h3("Usage Tips"),
    tags$ul(
    tags$li("If a large number of reads slows loading, consider reducing the 'Subsample Percentage' to improve responsiveness."),
    tags$li("Every sample-subsampling pair is considered as a unique object. If you want to remove a track, we always recomment using the 'Remove track' button in the IGV interface by clicking the gear icon on the right side of the track."),
    ),
    tags$img(src="helppage/remove-track.png", style = "display: block; margin: 20px auto; width: 550px; height: 300px;"),

    tags$hr(),

    # =========================================

    h2("FAQ"),
    
    # h3("1. What should I do if no records are found for my query?"),
    # p("Please double-check the gene symbol, gene ID, or transcript ID you entered. Ensure that it is a mouse gene and properly formatted. You may also try searching by a different identifier type."),
    
    h3("1. Why is the IGV browser not loading?"),
    p("If the IGV browser does not load properly, try refreshing the page."),
    
    tags$hr(),

    # ====================================

    h2("Acknowledgements"),

    p("We gratefully acknowledge the open-source projects and their developers whose contributions have made this database possible."),
    # p(" In particular, we recognize the essential roles of the following resources:")
    # tags$ul(
    # tags$li(tags$b("Shiny:"), " Providing the core framework for web application development in R."),
    # tags$li(tags$b("igvShiny:"), " Integrating the IGV genome browser within Shiny applications."),
    # tags$li(tags$b("DBI and RSQLite:"), " Enabling efficient management of local relational databases."),
    # tags$li(tags$b("DT:"), " Powering interactive and customizable data tables."),
    # tags$li(tags$b("Ggplot2, Plotly, and related visualization packages:"), " Supporting high-quality data visualizations throughout the platform.")
    # ),
    p("We extend our sincere thanks to the broader R and bioinformatics communities for their ongoing efforts in developing and maintaining these invaluable tools."),
    tags$hr(),

    # =====================================
    # h2("How to cite"),
    h2("Contact us"),
    p("Thank you for using the Mouse Poly(A) Tail Atlas Database! If you encounter any problem or have any question, please don't hesitate to contact us at ", tags$code("zhaijx@sustech.edu.cn"), "."),
  

  )
}
