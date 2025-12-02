homepageUI <- function(id) {
    ns <- NS(id)
    fluidPage(
      tags$style(HTML("
        .custom-text {
          margin-left: 10%;
          margin-right: 10%;
          text-align: justify;
        }
        p {
            font-size: 18px;
        }
      ")),
      h2("Welcome to Mouse Poly(A) Tail Atlas v1!"),
      tags$hr(),
      div(class = "custom-text",
        p("This website provides a comprehensive resource for exploring poly(A) tail length profiles across 18 different mouse organs under normal physiological conditions."),
        tags$img(
          src='homepage/mouse-atlas-home2.png',
          style = "display: block; margin: 20px auto; width: 550px; height: 400px;"
        ),
        p("Using FLEP-seq2, we have systematically mapped poly(A) tail dynamics at the gene and transcript level. Users can easily search by gene ID, transcript ID, or gene symbol to access poly(A) tail length and expression data across all organs. Interactive visualizations, including heatmaps, detailed tables, and an integrated IGV browser for poly(A) tail inspection, are available to support further exploration and analysis."),
        p("Our goal is to make poly(A) tail information readily accessible to the community, facilitating discoveries in gene regulation and mRNA biology. Detailed introduction of the website usage can be referred to the Help page.")
      )
    )
}