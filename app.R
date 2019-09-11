# TREND-DB app ------------------------------------------------------------

# Loading required packages -----------------------------------------------

library(shiny)
library(DT)
library(shinydashboard)
library(rentrez)
library(limma)
library(GO.db)
library(png)
library(ggmap)
library(Gviz)
library(org.Hs.eg.db)
library(shinyBS)
library(rtracklayer)
library(markdown)
library(igraph)
library(visNetwork)
library(magrittr)
library(rintrojs)


# Loading required input data ---------------------------------------------

message("Loading R data...")
load("trendseq.RData")

# these objects needed to be updated in newer Bioc versions
# txBygene <- updateObject(txBygene,verbose = FALSE)
utrs <- updateObject(utrs,verbose = FALSE)
# just in case, saving workspace
# save.image(file = "trendseq_bioc3_10.RData")
message("Done! \n")

# load annotation db
message("Loading TxDb...")
txdb <- loadDb("./data/hg38db_refGene.sqlite")
message("Done! \n")


# Loading the TRENDnetwork ------------------------------------------------

trendnet <-readRDS("TRENDnet.rds")

# processing that in brief
data <- toVisNetworkData(trendnet)
# visNetwork(nodes = data$nodes, edges = data$edges, height = "500px")
vis.nodes <- data$nodes
vis.links <- data$edges

vis.nodes$shape  <- "dot"  
vis.nodes$shadow <- TRUE # Nodes will drop shadow
vis.nodes$title  <- vis.nodes$label # Text on click
vis.nodes$size <- 4 * vis.nodes$size


# UI definition -----------------------------------------------------------
ui <- shinydashboard::dashboardPage(
  # header definition -------------------------------------------------------
  dashboardHeader(
    title = paste0("TREND-DB"),
    titleWidth = 900,
    # TODO: logo in the title?
    dropdownMenu(
      type="tasks",
      icon=icon("question-circle fa-1g"),
      badgeStatus=NULL,
      headerText="Help",
      notificationItem(
        text=actionButton(
          "dd_help", "First help",
          icon("hand-o-right")
        ),
        icon=icon(""), # tricking it to not have additional icon
        status="primary"
      ),
      notificationItem(
        text=actionButton(
          'dd_glossary', label="Open the glossary",
          icon=icon("book")
        ),
        icon=icon(""), status="primary"
      )
    ),
    dropdownMenu(
      type="tasks",
      icon=icon("info fa-1g"),
      badgeStatus=NULL,
      headerText="Additional information",
      notificationItem(
        text=actionButton(
          'session_info', label="About this session",
          icon=icon("window-maximize")
        ),
        icon=icon(""), status="primary"
      ),
      notificationItem(
        text=actionButton(
          'trenddb_info', label="About TREND-DB",
          icon=icon("heart")
        ),
        icon=icon(""), status="primary"
      )
    )
  ), # end of dashboardHeader
  # sidebar definition ------------------------------------------------------
  dashboardSidebar(
    disable = TRUE,
    width = 350), # end of dashboardSidebar
  # body definition ---------------------------------------------------------
  dashboardBody(
    introjsUI(),
    
    # adding some styling in the html elements
    shiny::tags$head(
      shiny::tags$style(
        HTML(
          '
          .skin-blue .wrapper, .content-wrapper {
          background-color: #e9e9e9;
          }
          
          #shiny-notification-panel {
          margin-bottom: 50px;
          margin-right: 7%;
          height: 50px;
          width: 400px;
          }
          
          .modal-lg {
          height: auto;
          width: 600px;
          }
          '
        )
      )
    ),
    
    ## main structure of the body for the dashboard
    tabBox(
      id = "tabs",
      width = 11,
      selected = "Welcome",
      
      #h4("Session Info"),
      #verbatimTextOutput("sessioninfo")),

      # ui panel welcome --------------------------------------------------------
      tabPanel(
        "Welcome",
        icon = icon("home"),
        actionButton(
          "tour_firststeps", "Click me for a quick tour",
          icon("hand-o-right")
        ),
        includeMarkdown("trenddb_welcomepage.md")
      ), # end of Welcome panel

      # ui panel data preview ---------------------------------------------------
      tabPanel(
        "Data Preview",
        icon = icon("eye"),
        actionButton(
          "tour_datapreview", "Click me for a quick tour",
          icon("hand-o-right")
        ),
        h2("Inspect Matrix"),
        fluidRow(
          column(
            width = 6,
            uiOutput(
              "inspectMatrix_desc"
            )
          ),
          column(
            width = 3,
            actionButton("continueMV", "Continue to Main View", icon = icon("arrow-right"))
          )
        ),
        shiny::tags$style(type = 'text/css', "#continueMV {position: absolute; right: 15px;}"),
        br(),
        fluidRow(
          column(
            width = 9,
            DT::dataTableOutput("inspectMatrix"),
            br(),
            div(
              align = "right",
              style = "margin-right:15px; margin-bottom:10px",
              sliderInput(
                "dirFilter",
                "Filter by shortening index:",
                min = 0,
                max = 15,
                value = 0,
                step = 0.5
              )
            ),
            htmlOutput("selection"),
            DT::dataTableOutput("overviewTable")
          )
        )
      ), # end of Data Preview panel

      # ui panel main view ------------------------------------------------------
      tabPanel(
        "Main View",
        icon = icon("table"),
        fluidRow(
          actionButton(
            "tour_mainview", "Click me for a quick tour",
            icon("hand-o-right")
          )
        ),
        fluidRow(
          column(
            width = 10,
            offset = 1,
            visNetworkOutput("trendnetwork",width = "800px",height = "800px")
          )
        ),
        fluidRow(
          column(
            width = 6,
            p(h2('Gene View')),
            fluidRow(column(
              width = 7,
              uiOutput("mainViewGene_desc"),
              br()
            )),
            fluidRow(
              column(
                width = 6,
                selectInput(
                  "geneInput",
                  "Gene",
                  choices = c("", genenames),
                  selected = NULL
                )
              ),
              column(
                width = 4,
                conditionalPanel(
                  condition = "output.geneInfoMain",
                  actionButton("continueGP", "Continue to Gene Plot", icon =
                                 icon("arrow-right"))
                )
              )
            ),
            shiny::tags$style(
              type = 'text/css',
              "#continueGP {position: absolute; right: 15px; margin-top: 25px;}"
            ),
            br(),
            conditionalPanel(
              condition = "output.geneInfoMain",
              h4("Gene Summary:")
            ),
            DT::dataTableOutput("geneInfoMain"),
            br(),
            textOutput("geneSummary"),
            br(),
            htmlOutput("NCBImain"),
            br(),
            conditionalPanel(
              condition = "output.geneTableOutput",
              br(),
              uiOutput("involvedConds_desc")
            ),
            DT::dataTableOutput("geneTableOutput"),
            br()
          ),
          column(
            width = 6,
            p(h2('Condition View')),
            fluidRow(
              column(
                width = 7,
                uiOutput("mainViewCond_desc"),
                br()
              )
            ),
            fluidRow(
              column(
                width = 6,
                selectInput(
                  "kdInput",
                  "Condition",
                  choices = c("", paste0(conditions, "_kd")),
                  selected = NULL
                )
              ),
              br(),
              conditionalPanel(
                condition = "output.kdTableOutput",
                column(
                  width = 2,
                  actionButton("showScatterplot", "Scatterplot", icon = icon("external-link"))
                ),
                column(
                  width = 2,
                  actionButton(
                    "goToGoana",
                    "GO Enrichment",
                    icon = icon("arrow-right"),
                    onclick = "location.href='#goana';"
                  )
                )
              )
            ),
            conditionalPanel(
              condition = "output.kdTableOutput",
              br(),
              uiOutput("affectedGenes_desc")
            ),
            DT::dataTableOutput("kdTableOutput"),
            bsModal(
              "scatterModal",
              "Plot",
              "showScatterplot",
              imageOutput("scatterplot"),
              size = "large"
            ),
            conditionalPanel(
              condition = "output.kdTableOutput",
              br(),br(),
              uiOutput("goana_desc"),
              br(),
              actionButton("goanaSubmit", "Goana"),
              br(),br(),
              DT::dataTableOutput("goanaTable")
            )
          )
        )
      ), # end of panel Main view

      # ui panel gene plot ------------------------------------------------------
      tabPanel(
        "Gene Plot",
        icon = icon("bar-chart"),
        fluidRow(
          column(
            width = 11,
            h2("Gene Plot"),
            actionButton(
              "tour_geneplot", "Click me for a quick tour",
              icon("hand-o-right")
            ),
            fluidRow(
              column(
                width = 7,
                uiOutput("genePlot_desc")
              ),
              column(
                width = 3,
                fluidRow(
                  column(
                    width = 8,
                    conditionalPanel(
                      condition = "!output.genePlotInfo",
                      actionButton("selectGene", "Select gene", icon =
                                     icon("search"))
                    ),
                    conditionalPanel(
                      condition = "output.genePlotInfo",
                      actionButton("switchGene", "Change gene", icon =
                                     icon("refresh"))
                    )
                  )
                ),
                br(),br(),
                fluidRow(
                  column(
                    width = 8,
                    conditionalPanel(
                      condition = "output.genePlotInfo",
                      actionButton("viewBrowser",
                                   "View in Genome Browser",
                                   icon("search"))
                    )
                  )
                )
              )
            ),
            # shiny::tags$style(type = 'text/css', "#selectGene {position: absolute; right: 15px;}"),
            #  shiny::tags$style(type = 'text/css', "#switchGene {position: absolute; right: 15px;}"),
            br(),
            fluidRow(
              column(
                width = 10,
                DT::dataTableOutput("genePlotInfo")
              )
            ),
            br(),
            htmlOutput("NCBIplot"),
            br(),
            conditionalPanel(
              condition = "output.genePlotInfo",
              h4("Plot Options:"),
              br(),
              fluidRow(
                column(
                  width = 9,
                  fluidRow(
                    column(
                      width = 4,
                      selectInput("genePlotCond",
                                  "Condition (default = max SI)",
                                  choices = "max SI")
                    ),
                    column(
                      width = 4,
                      radioButtons(
                        "selectView",
                        "Plot View (default = 3' UTR)",
                        choices = c("3' UTR", "Gene Body"),
                        selected = "3' UTR"
                      )
                    )
                  ),
                  br(),
                  actionButton(
                    "plotSubmit",
                    "Create/Update Plot",
                    width = 250,
                    icon("paint-brush"),
                    style = "color: #333; background-color: #88b5dd; border-color: #88b5dd; font-weight: bold; font-size: medium"
                  )
                )
              )
            ),
            br(),
            div(
              align = "left",
              style = "width:800px;",
              plotOutput("gviz"),
              br(),
              conditionalPanel(
                condition = "output.gviz",
                div(
                  align = "right",
                  style = "margin-right:15px; margin-bottom:10px",
                  downloadButton("download_genePlot", "Download Plot")
                  # textInput("filename_genePlot", label = "Save as...", value = "genePlot.pdf")
                )
              )
            ),
            conditionalPanel(
              condition = "output.genePlotInfo",
              sliderInput(
                "simAffectedSlider",
                "Number of similarily affected genes",
                min = 0,
                max = 30,
                value = 0,
                step = 1
              ),
              #br(),
              # actionButton("simInfoUpdate", "Update similar"),
              br(),
              htmlOutput("simAffected_desc"),
              br(),
              DT::dataTableOutput("simAffected_table")
            )
          )
        )
      ), # end of Gene Plot panel

      # ui panel genome browser -------------------------------------------------
      tabPanel(
        "Genome Browser",
        icon = icon("search"),
        h2("Genome Browser"),
        fluidRow(
          column(
            width = 8,
            actionButton(
              "tour_genomebrowser", "Click me for a quick tour",
              icon("hand-o-right")
            ),
            uiOutput(
              "genomeBrowser_desc"
            )
          )
        ),
        br(),br(),
        fluidRow(
          column(
            width = 2,
            uiOutput("geneSelection"),
            actionButton("changeGene", "Gene", icon = icon("search")),
            br(),
            br(),
            selectInput(
              "condsGB",
              "Condition",
              choices = c("", paste0(conditions, "_kd")),
              selected = NULL
            ),
            br(),
            #checkboxGroupInput("selectConditions",
            #                  "Conditions",
            #                   choices = c("", paste0(conditions, "_kd"))),
            checkboxGroupInput("selectTracks",
                               "Additional Tracks",
                               choices = list("miRNAs" = 1)),
            br(), br(),
            htmlOutput("browserNewWindow")
            # actionButton("tracksInput", "Update Tracks", icon = icon("refresh"))
          ),
          column(
            width = 10,
            htmlOutput("browserFrame")
          )
        )
      ), # end of Genome Browser panel

      # ui panel about ----------------------------------------------------------
      tabPanel(
        "About", 
        icon = icon("info-circle"),
        h2("About"),
        fluidRow(
          column(
            width = 11,
            uiOutput("about")
          )
        )
      ) # end of About panel
    ) # end of tabBox 
  ) # end of dashboardBody
)

# Server definition -------------------------------------------------------
server <- function(input, output, session) {
  v <- reactiveValues(clearGoana = TRUE, clearPlot = TRUE)
  
  # renders inspect matrix containing all genes + conditions with their respective dir indexes
  output$inspectMatrix <-
    DT::renderDataTable(server = TRUE, {
      cat(file = stderr(), "Loading inspect matrix...")
      table <- KDmat
      colnames(table) <-
        paste0(colnames(KDmat), "_kd")
      rownames(table) <- rownames(KDmat)
      table <- as.data.frame(table)
      cutoff <- input$dirFilter
      table <-
        subset(table, apply(table, MARGIN = 1, function(x)
          any(x > cutoff | x < -cutoff)))
      cat(file = stderr(), "Done! \n")
      DT::datatable(
        table,
        # selection = list(mode = 'single', target = 'cell'),
        colnames = c("Gene" = 1),
        rownames = TRUE,
        selection = "single",
        extensions = c('FixedColumns', 'Scroller'),
        options = list(
          scrollX = TRUE,
          deferRender = TRUE,
          scrollY = 407,
          scroller = TRUE,
          pageLength = 50,
          fixedColumns = FALSE
        )
      ) %>%
        formatRound(colnames(table), digits = 3) %>%
        formatStyle("Gene",
                    backgroundColor = "#ffeecc",
                    fontWeight = "bold")
    })
  
  # selection in inspect matrix (single)
  output$selection <- renderUI({
    req(input$inspectMatrix_rows_selected)
    cutoff <- input$dirFilter
    table <- as.data.frame(KDmat)
    if (input$dirFilter > 0) {
      table <-
        subset(table, apply(table, MARGIN = 1, function(x)
          any(x > cutoff | x < -cutoff)))
    }
    selectedRowIndex <-
      input$inspectMatrix_rows_selected
    selectedRow <- table[selectedRowIndex, ]
    gene <- rownames(table)[selectedRowIndex]
    h4("Selected Gene: ", gene)
  })
  
  output$trendnetwork <- renderVisNetwork({
    visnet <- visNetwork(vis.nodes, vis.links, height = "1000px", width = "1000px")
    visnet %>% 
      visOptions(highlightNearest = TRUE, selectedBy = "Group",nodesIdSelection = TRUE) %>% 
      visLegend(main="Legend", position="right", ncol=1) 
  })
  
  output$geneSelection <- renderUI({
    HTML("<div><strong>Gene</strong> <br>",
         input$geneInput,
         "<br><br></div>")
  })
  
  # genome browser iFrame
  output$browserFrame <- renderUI({
    iframe <-
      shiny::tags$iframe(
        src = browserURL(),
        height = 800,
        width = "100%",
        seamless = NA
      )
    iframe
  })
  
  # renders overview table for selected gene in inspect matrix
  output$overviewTable <- DT::renderDataTable({
    table <- overviewTable()
    DT::datatable(
      table,
      colnames = c("Condition" = 1),
      rownames = paste0(rownames(table), "_kd"),
      selection = "single",
      options = list(
        scrollY = 300,
        scrollCollapse = TRUE,
        paging = FALSE,
        searching = TRUE
      )
    ) %>%
      formatStyle("Condition",
                  backgroundColor = "#ffe0cc",
                  fontWeight = "bold")
  })
  
  # renders gene table
  output$geneTableOutput <- DT::renderDataTable({
    req(input$geneInput)
    table <- geneTable()
    DT::datatable(
      table,
      colnames = c("Condition" = 1),
      rownames = paste0(rownames(table), "_kd"),
      selection = "single",
      extensions = c('FixedColumns'),
      options = list(
        scrollX = TRUE,
        scrollY = 300,
        paging = FALSE,
        searching = TRUE,
        fixedColumns = FALSE
      )
    ) %>%
      formatStyle("Condition",
                  backgroundColor = "#ffe0cc",
                  fontWeight = "bold")
  })
  
  # renders gene info table in Main View
  output$geneInfoMain <- DT::renderDataTable({
    req(input$geneInput)
    DT::datatable(
      geneInfoTableMain(),
      extensions = c('FixedColumns'),
      options = list(
        paging = FALSE,
        searching = FALSE,
        bInfo = 0,
        scrollX = TRUE,
        fixedColumns = FALSE
      )
    )  %>%
      formatStyle("Gene",
                  backgroundColor = "#ffeecc",
                  fontWeight = "bold")
  })
  
  # description text for gene in Main View (input$geneInput, using rentrez)
  output$geneSummary <- renderPrint({
    req(input$geneInput)
    gene_id <- sym2eg(input$geneInput)
    sum <- ""
    try(x <- entrez_summary(db = "gene", id = gene_id))
    try(sum <- x$summary)
    cat(sum)
  })
  
  # renders NCBI link for Main View
  output$NCBImain <- renderUI({
    a("View Gene (NCBI)", href = NCBImain(), target = "_blank")
  })
  
  # renders GB link
  output$browserNewWindow <- renderUI({
    a("Open in new window!", href=browserURL(), target="_blank") 
  })
  
  # renders similarily affected genes table
  output$simAffected_table <- DT::renderDataTable({
    req(input$geneInput)
    req(input$simAffectedSlider != 0)
    table <- getSimAffected()
    colnames(table) <-
      paste0(colnames(table), "_kd")
    DT::datatable(
      table,
      colnames = c("Gene" = 1),
      options = list(
        scrollX = TRUE,
        # scrollY = 300,
        searching = TRUE,
        paging = FALSE
      )
    )  %>%
      formatStyle("Gene",
                  backgroundColor = "#ffeecc",
                  fontWeight = "bold")
  })
  
  # renders condition table
  output$kdTableOutput <- DT::renderDataTable({
    req(input$kdInput)
    DT::datatable(
      kdTable(),
      colnames = c("Gene" = 1),
      selection = "single",
      extensions = c('FixedColumns'),
      options = list(
        scrollX = TRUE,
        scrollY = 480,
        paging = FALSE,
        searching = TRUE,
        fixedColumns = FALSE
      )
    )  %>%
      formatStyle("Gene",
                  backgroundColor = "#ffeecc",
                  fontWeight = "bold")
  })
  
  output$goanaTable <- DT::renderDataTable({
    req(input$kdInput)
    if (v$clearGoana) {
      return()
    } else {
      DT::datatable(
        withProgress(
          message = 'Calculating',
          detail = 'This may take a while...',
          value = 0.8,
          {
            createGoana()
          }
        ),
        colnames = c("GO" = 1),
        extensions = c('FixedColumns'),
        options = list(
          bInfo = 0,
          scrollX = TRUE,
          scrollY = 480,
          fixedColumns = FALSE,
          order = list(5, 'asc')
        )
      )  %>%
        formatStyle("GO",
                    backgroundColor = "#b3ccff",
                    fontWeight = "bold")
    }
  })
  
  # renders info table for gene plot
  output$genePlotInfo <- DT::renderDataTable({
    req(input$geneInput)
    DT::datatable(genePlotTable(),
                  options = list(
                    paging = FALSE,
                    searching = FALSE,
                    bInfo = 0
                  ))  %>%
      formatStyle("Gene",
                  backgroundColor = "#ffeecc",
                  fontWeight = "bold") %>%
      formatStyle("Condition",
                  backgroundColor = "#ffe0cc",
                  fontWeight = "bold")
  })
  
  # renders gene plot
  output$gviz <- renderPlot({
    if (v$clearPlot) {
      return()
    } else {
      withProgress(message = 'Calculating',
                   detail = 'This may take a while...',
                   value = 0.8,
                   {
                     drawPlot()
                   })
      # }, height = function() {
      #   session$clientData$output_gviz_width
    }
  })
  
  # creates a download handler for gene plot
  output$download_genePlot <- downloadHandler(
    filename = "genePlot.png",
    content = function(file) {
      png(file, width = 800)
      gene <- input$geneInput
      if (input$genePlotCond == "max SI") {
        cond <- maxDI(gene)
      } else {
        cond <- strsplit(input$genePlotCond, "_")[[1]][1]
      }
      if (input$selectView == "3' UTR") {
        view = "utr"
      } else if (input$selectView == "Gene Body") {
        view = "gb"
      }
      try(genePlot(gene, cond, view))
      dev.off()
    }
  )
  
  # renders NCBI link for Gene Plot
  output$NCBIplot <- renderUI({
    a("View Gene (NCBI)", href = NCBImain(), target = "_blank")
  })
  
  # creates download link for gene plot
  output$pdflink <-
    downloadHandler(filename <- "myplot.pdf",
                    content <-
                      function(file) {
                        file.copy("plot.pdf", file)
                      })
  
  # description of About tab
  output$about <- renderUI({
    HTML(
      "<div><br>This page contains descriptions for the different features of the TREND-DB. </div><br>
      <h4>Data Preview</h4>
      <div>The <strong>Data Preview</strong> tab shows an overview of all shortening index values of each gene across all conditions. With the
      '<strong>Filter by shortening index</strong>' slider, the table can be subsetted to only show genes with a shortening index value (absolute) above the threshold
      in at least one condition. Upon selecting a row (= gene) in the table, a summary table of the selected gene is printed, showing the conditions in which the gene is affected
      with the respective shortening indices and p-values.</div> <br> <br>
      <h4>Main View</h4>
      <div>The <strong>Main View</strong> tab is divided into <strong>Condition View</strong> and <strong>Gene View</strong>.
      In the <strong>Condition View</strong>, users can select a condition to view a table containing all genes that are affected in this condition.
      By clicking the <strong>Goana</strong> button, an over-representation analysis for Gene Ontology terms in a list of Entrez Gene IDs
      is performed and the results are shown in a table.  The list of genes affected in the selected condition is used as a gene set with all genes as a background.
      In the <strong>Gene View</strong>, users can select a gene to print a short summary and a table containing all conditions where the selected gene is affected.</div><br><br>
      <h4>Gene Plot</h4>
      <div>In the <strong>Gene Plot</strong> tab, users can plot the selected gene (using the <strong>Gviz</strong> package). Aside from the reference track of the selected gene and its chromosome, the control track as well as the selected condition are shown.
      The default condition is the condition with the maximum SI for this gene, it can however be changed using a drop-down menu.
      Adjusting the <strong> Number of similarily affected</strong> slider will print genes that are affected in a similar way as the selected gene (with the amount of them determined by the slider).
      Similarily affected genes are obtained by computing the distance matrix (using the shortening indices of all genes across conditions) and looking for genes with the smallest distances.</div><br><br>
      <h4>Genome Browser</h4>
      <div> The <strong>Genome Browser</strong> tab includes an instance of the UCSC Genome Browser, where users can take a more in-depth look at the data. </div><br><br>
      <h4>Contact</h4>
      <div>Denise Scherzinger <a href='mailto:denscher@uni-mainz.de' target='_top'>(denscher@uni-mainz.de)</a></div>
      "
    )
  })
  
  # description of Inspect Matrix
  output$inspectMatrix_desc <- renderUI({
    HTML(
      "This is an overview of all shortening index values of each gene across all conditions. <br> Adjust
      the '<strong>Filter by shortening index</strong>' slider to only show genes with a shortening index value (absolute) above the threshold
      in at least one condition. <br> Select a row in the table to print a summary table of the selected gene, showing the conditions in which it is affected
      with the respective shortening indices and p-values."
    )
  })
  
  # description of Gene View in Main View tab
  output$mainViewGene_desc <- renderUI({
    HTML(
      "<div>Select a <strong>gene</strong> to print a gene summary and a table containing all conditions where this gene is affected.
      </div>"
    )
  })
  
  # description of Genome Browser tab
  output$genomeBrowser_desc <- renderUI({
    HTML(
      "<div>Please wait for data to be fetched from UCSC (~10s). The Genome Browser shows the currently selected gene and contains data tracks for all conditions which are grouped under 'Trendseq UCSC Hub'.
      To show additional tracks for microRNA target sites, select the checkboxes on the left.</div>"
    )
  })
  
  # description of affected genes
  output$affectedGenes_desc <- renderUI({
    HTML(
      "<div><h4>Affected Genes:</h4><br>",
      paste0(
        "The following genes are affected in the <strong>",
        input$kdInput
      ),
      "</strong> condition.</div>"
    )
  })
  
  # description of Condition View in Main View tab
  output$mainViewCond_desc <- renderUI({
    HTML(
      "<div>Select a <strong>condition</strong> to print a table showing all genes affected in this condition.</div>"
    )
  })
  
  # description of involved conditions
  output$involvedConds_desc <- renderUI({
    HTML(
      "<div><h4>Involved Conditions:</h4><br><strong>",
      paste0(
        input$geneInput,
        "</strong> is involved in the following conditions. </div>"
      )
    )
  })
  
  # description of Goana method
  output$goana_desc <- renderUI({
    HTML(
      "<div><h4><a name='goana'></a>GO Enrichment:</h4><br> The <strong>Goana</strong> (Limma Package) method performs an over-representation analysis for Gene Ontology terms in a list of Entrez Gene IDs.
      The list of genes affected in the selected condition is used as a gene set with all genes as a background. <br>
      Columns in the resulting table show the following: <br>
      <strong>Term</strong> GO Term <br> <strong>Ont</strong> ontology that the GO term belongs to ('BP', 'CC' or 'MF')<br>
      <strong>N</strong> number of genes in the GO term <br> <strong>DE</strong> number of genes in the provided gene set<br>
      <strong>P.DE</strong> p-value for over-representation of the GO term in the set<br><br>Click the <strong>Goana</strong>
      button to create the table.</div>"
    )
  })
  
  # description of similarily affected genes selection
  output$simAffected_desc <- renderUI({
    req(input$geneInput)
    req(input$simAffectedSlider != 0)
    HTML(
      "<div><h4>Similarily affected genes</h4>",
      paste0(
        "Similarily affected genes are obtained by computing the distance matrix (using the shortening indices of all genes across conditions) and looking for genes with the smallest distances.<br>", 
        "The following genes are affected in a similar way as <strong>",
        input$geneInput,
        "</strong>:</div>"
      )
    )
  })
  
  # description of Gene Plot tab
  output$genePlot_desc <- renderUI({
    HTML(
      "<div>Plots in this tab are created using the <strong>Gviz</strong> package. Aside from the reference track of the selected gene and its chromosome, the control track as well as the selected condition are shown. <br>
      The default condition is the condition with the maximum SI for this gene, it can however be changed using the drop-down menu.<br>
      Adjusting the <strong> Number of similarily affected</strong> slider will print genes that are affected in a similar way as the selected gene (with the amount of them determined by the slider). </div>"
    )
  })
  
  # sessioninfo
  output$sessioninfo <- renderPrint({
    sessionInfo()
  })
  
  # proxy object for similarily affected table (used to clear selection)
  proxy = dataTableProxy("simAffectedTable")
  
  # overview table for selected gene in inspect matrix, only created when a row is selected
  overviewTable <-
    eventReactive(input$inspectMatrix_rows_selected, {
      cutoff <- input$dirFilter
      table <- as.data.frame(KDmat)
      if (input$dirFilter > 0) {
        table <-
          subset(table, apply(table, MARGIN = 1, function(x)
            any(x > cutoff | x < -cutoff)))
      }
      selectedRowIndex <-
        input$inspectMatrix_rows_selected
      selectedRow <- table[selectedRowIndex, ]
      gene <- rownames(table)[selectedRowIndex]
      overviewTable <- createGeneTable(gene)
      
      updateSelectInput(session,
                        "geneInput",
                        choices = c("", genenames),
                        selected = gene)
      return(overviewTable)
    })
  
  # reactive function for affected genes, updated on input change (input$kdInput)
  affectedGenes <- reactive({
    kd <- strsplit(input$kdInput, "_")[[1]][1]
    getGenes(kd)
  })
  
  # reactive function for involved conditions, updated on input change (input$geneInput)
  conditionsInvolved <- reactive({
    getConds(input$geneInput)
  })
  
  # creates gene table, updated on input change in Main View (input$geneInput)
  geneTable <- reactive({
    createGeneTable(input$geneInput)
  })
  
  # creates condition table containing all affected genes with dir index, p-value prox + p-value dist, updated on input change in Main View (input$kdInput)
  kdTable <- reactive({
    affectedGenes <- affectedGenes()
    kdTable <- matrix(, nrow = 0, ncol = 3)
    colnames(kdTable) <-
      c("shortening_index", "pvalue_short", "pvalue_long")
    for (gene in affectedGenes) {
      x <- rep(NA, 3)
      kdTable <- rbind(kdTable, x)
      rownames(kdTable)[rownames(kdTable) == "x"] <-
        gene
      kd <- strsplit(input$kdInput, "_")[[1]][1]
      kdTable[gene, "shortening_index"] <-
        KDmat[gene, kd]
      kdTable[gene, "pvalue_short"] <-
        format(pValPmat[gene, kd], scientific = TRUE)
      kdTable[gene, "pvalue_long"] <-
        format(pValDmat[gene, kd], scientific = TRUE)
    }
    return(kdTable)
  })
  
  # creates output table for goana GO term enrichment analysis
  createGoana <- eventReactive(input$goanaSubmit, {
    req(input$kdInput)
    cat(file = stderr(), "Calculating GO Enrichment table...")
    affectedGenes <- affectedGenes()
    geneSet <- c()
    for (gene in affectedGenes) {
      eg <- sym2eg(gene)
      geneSet <- append(geneSet, eg)
    }
    background <- c()
    for (gene in genenames) {
      eg <- sym2eg(gene)
      background <- append(background, eg)
    }
    
    table <-
      limma::goana(geneSet,
                   species = "Hs",
                   FDR = 0.05,
                   universe = background)
    goanaTable <- table[table$P.DE < 0.05,]
    cat(file = stderr(), "Done! \n")
    return(goanaTable)
  })
  
  getSimAffected <- reactive({
    req(input$simAffectedSlider != 0 && !is.null(input$geneInput))
    return(simAffected(input$geneInput, input$simAffectedSlider))
  })
  
  # window for scatterplots
  scatterModal <- eventReactive(input$showScatterplot, {
    cond <- strsplit(input$kdInput, "_")[[1]][1]
    file <- paste0("./data/scatterplots/", cond, ".png")
    return(list(src = file, contentType = "image/png", width = 400, height = 400))
  })
  
  output$scatterplot <- renderImage({
    scatterModal()
  }, deleteFile = FALSE)
  
  # calls genePlot for selected gene (input$geneInput), only called after action button is used (input$plotSubmit)
  drawPlot <- eventReactive(input$plotSubmit, {
    req(input$geneInput)
    gene <- input$geneInput
    if (input$genePlotCond == "max SI") {
      cond <- maxDI(gene)
    } else {
      cond <- strsplit(input$genePlotCond, "_")[[1]][1]
    }
    if (input$selectView == "3' UTR") {
      view = "utr"
    } else if (input$selectView == "Gene Body") {
      view = "gb"
    }
    id <- paste0(gene, "_", cond, "_", view)
    file <- paste0("./data/plotCache/", id, ".png")
    cat(file = stderr(), "Checking plot cache...")
    if (isCached(id)) {
      cat(file = stderr(), "Plot cached! Loading...\n")
      img <- readPNG(file)
      myplot <- ggimage(img)
      print(myplot)
      cat(file = stderr(), "Done!\n")
    } else{
      cat(file = stderr(), "Plot not yet cached! Calculating...\n")
      try(genePlot(gene, cond, view))
      #cachePlot(file, gene, cond, view)
      cat(file = stderr(), "Done!\n")
    }
  })
  
  # info table for gene plot, updated on input change (input$geneInput)
  genePlotTable <- reactive({
    if (input$genePlotCond == "max SI") {
      cond <- maxDI(input$geneInput)
    } else {
      cond <- strsplit(input$genePlotCond, "_")[[1]][1]
    }
    genePlotInfoTable(input$geneInput, cond)
  })
  
  # creates gene info table in Main View, updated on input change (input$geneInput)
  geneInfoTableMain <- reactive({
    return(geneInfoTable(input$geneInput))
  })
  
  # creates NCBI link for gene (input$geneInput)
  NCBImain <- reactive({
    req(input$geneInput)
    return(paste0(
      "http://www.ncbi.nlm.nih.gov/gene/",
      sym2eg(input$geneInput)
    ))
  })
  
  browserURL <- reactive({
    return(createBrowserURL(input$geneInput))
  })
  

  # all observers -----------------------------------------------------------
  
  # update knockdown condition if selected via trendnetwork
  observe({
    x <- input$trendnetwork_selected
    
    # Can use character(0) to remove all choices
    if (is.null(x))
      x <- character(0)
    
    # Can also set the label and select items
    updateSelectInput(session, "kdInput",
                      selected = paste0(x,"_kd")
    )
  })
    
  # clear Goana table if condition selection is changed
  observeEvent(input$kdInput, {
    v$clearGoana <- TRUE
  }, priority = 10)
  
  observeEvent(input$goanaSubmit, {
    v$clearGoana <- FALSE
  }, priority = 10)
  
  # clear plot if condition selection is changed
  observeEvent(input$genePlotCond, {
    v$clearPlot <- TRUE
  }, priority = 10)
  
  observeEvent(input$geneInput, {
    v$clearPlot <- TRUE
    updateSelectInput(session, "selectView",
                      selected = "3' UTR")
  }, priority = 10)
  
  observeEvent(input$plotSubmit, {
    v$clearPlot <- FALSE
  }, priority = 10)
  
  # observer to change selected tab to Main View on button press (when gene has been selected in Inspect Matrix)
  observe({
    if (input$continueMV > 0) {
      updateTabItems(session, "tabs", selected = "Main View")
    }
  })
  
  # observer to change selected tab to Gene Plot on button press (when gene has been selected in Main View)
  observe({
    if (input$continueGP > 0) {
      updateTabItems(session, "tabs", selected = "Gene Plot")
    }
  })
  
  # observers to change selected tab to Main View (in order to select/change a gene to plot)
  observe({
    if (input$selectGene > 0) {
      updateTabItems(session, "tabs", selected = "Main View")
      updateSelectInput(session, "selectView",
                        selected = "3' UTR")
    }
  })
  
  # observer to change selected tab to Main View
  observe({
    if (input$switchGene > 0) {
      updateTabItems(session, "tabs", selected = "Main View")
    }
  })
  
  # observer to change selected tab to Main View
  observe({
    if (input$changeGene > 0) {
      updateTabItems(session, "tabs", selected = "Main View")
    }
  })
  
  # observer to change selected tab to Genome Browser
  observe({
    if (input$viewBrowser > 0) {
      updateTabItems(session, "tabs", selected = "Genome Browser")
    }
  })
  
  # observer to update gene selection in Main View when gene is selected in condition table
  observe({
    if (length(input$kdTableOutput_rows_selected) > 0) {
      table <- kdTable()
      selectedRowIndex <-
        input$kdTableOutput_rows_selected[1]
      selectedRow <- table[selectedRowIndex, ]
      gene <- rownames(table)[selectedRowIndex]
      updateSelectInput(session,
                        "geneInput",
                        selected = gene)
    }
  })
  
  observe({
    if (length(input$geneTableOutput_rows_selected) > 0) {
      table <- geneTable()
      rownames(table) <-
        paste0(rownames(table), "_kd")
      selectedRowIndex <-
        input$geneTableOutput_rows_selected[1]
      selectedRow <- table[selectedRowIndex, ]
      cond <- rownames(table)[selectedRowIndex]
      updateSelectInput(session,
                        "genePlotCond",
                        selected = cond)
      updateSelectInput(session,
                        "kdInput",
                        selected = cond)
      updateSelectInput(session,
                        "condsGB",
                        selected = cond)
    }
    else {
      updateSelectInput(session,
                        "genePlotCond",
                        selected = "max SI")
    }
  })
  
  observe({
    if (length(input$overviewTable_rows_selected) > 0) {
      table <- overviewTable()
      rownames(table) <-
        paste0(rownames(table), "_kd")
      selectedRowIndex <-
        input$overviewTable_rows_selected[1]
      selectedRow <- table[selectedRowIndex, ]
      cond <- rownames(table)[selectedRowIndex]
      updateSelectInput(session,
                        "kdInput",
                        selected = cond)
      updateSelectInput(session,
                        "condsGB",
                        selected = cond)
    }
  })
  
  # updates condition selection choices in Gene Plot for selected gene (input$genePlotInput)
  updatePlotConds <- observe({
    req(input$geneInput)
    conds <- getConds(input$geneInput)
    updateSelectInput(
      session,
      "genePlotCond",
      choices = c("max SI", paste0(conds, "_kd")),
      selected = "max SI"
    )
    updateSelectInput(session,
                      "condsGB",
                      choices = c("max SI", paste0(conds, "_kd")),
                      selected = "max SI")
  })
  
  observeEvent(input$tour_firststeps, {
    tour <- read.delim("tours/intro_firststeps.txt",
                       sep=";", stringsAsFactors=FALSE, row.names=NULL, quote="")
    introjs(session, options=list(steps=tour))
  })
  
  observeEvent(input$tour_datapreview, {
    tour <- read.delim("tours/intro_datapreview.txt",
                       sep=";", stringsAsFactors=FALSE, row.names=NULL, quote="")
    introjs(session, options=list(steps=tour))
  })
  
  observeEvent(input$tour_mainview, {
    tour <- read.delim("tours/intro_mainview.txt",
                       sep=";", stringsAsFactors=FALSE, row.names=NULL, quote="")
    introjs(session, options=list(steps=tour))
  })
  
  observeEvent(input$tour_geneplot, {
    tour <- read.delim("tours/intro_geneplot.txt",
                       sep=";", stringsAsFactors=FALSE, row.names=NULL, quote="")
    introjs(session, options=list(steps=tour))
  })
  
  observeEvent(input$tour_genomebrowser, {
    tour <- read.delim("tours/intro_genomebrowser.txt",
                       sep=";", stringsAsFactors=FALSE, row.names=NULL, quote="")
    introjs(session, options=list(steps=tour))
  })
  
  observeEvent(input$dd_help, {
    showModal(modalDialog(
      title="First help", size="l",fade=TRUE,
      footer=NULL, easyClose=TRUE,
      tagList(
        HTML("Some first help"),
        renderPrint({
          sessionInfo()
        })
      )
    ))
  })
  
  observeEvent(input$dd_glossary, {
    showModal(modalDialog(
      title="The TREND-DB glossary", size="l",fade=TRUE,
      footer=NULL, easyClose=TRUE,
      tagList(
        includeMarkdown("trenddb_glossary.md")
      )
    ))
  })
  
  
  observeEvent(input$session_info, {
    showModal(modalDialog(
      title="Session information", size="l",fade=TRUE,
      footer=NULL, easyClose=TRUE,
      tagList(renderPrint({
        sessionInfo()
      }))
    ))
  })
  
  observeEvent(input$trenddb_info, {
    showModal(modalDialog(
      title="About TREND-DB", size="m", fade=TRUE,
      footer=NULL, easyClose=TRUE,
      tagList(
        br(), br(),
        HTML("If you use this application for your projects, please use the following citation information:"),
        renderPrint({
          citation("iSEE")
        })
      )
    ))
  })

  # function definitions ----------------------------------------------------

  # function to get all affected genes in this condition (cond)
  getGenes <- function(cond) {
    allgenes_selectedKD <-
      KDmat[, colnames(KDmat) %in% cond, drop = FALSE] # force staying as a matrix
    affectedGenes <-
      rownames(allgenes_selectedKD)[!is.na(allgenes_selectedKD)]
    return(affectedGenes)
  }
  
  # function to determine specified number (n) of similarily affected genes for this gene
  simAffected <- function(gene, n) {
    if (!(gene %in% genenames))
      return()
    List <- list()
    for (i in 1:n) {
      x <-
        which(dist == sort(dist[, gene], partial = i)[i], arr.ind = TRUE)
      List[[i]] <- row.names(x)
    }
    top <- KDmat0[unlist(List),]
    top <- unique(top[, ])
    return(top[!rownames(top) %in% gene, , drop = FALSE])
  }
  
  # function to get all involved conditions for this gene
  getConds <- function(gene) {
    allconds_selectedGene <-
      KDmat[(rownames(KDmat) %in% gene), , drop = FALSE] # force staying as a matrix
    involvedConds <-
      colnames(allconds_selectedGene)[!is.na(allconds_selectedGene)]
    return(involvedConds)
  }
  
  # function to convert gene symbol (gene) into entrez ID
  sym2eg <- function(gene) {
    if (length(sym2eg_list) > 0) {
      gene_id <- sym2eg_list[gene]
      gene_id <- as.character(gene_id[1])
    }
    return(gene_id)
  }
  
  # get chromosome of this gene (using txdb)
  chromosome <- function(gene) {
    gene_id <- sym2eg(gene)
    # x <- entrez_summary(db = "gene", id = gene_id)
    # chrom <- paste0("chr", x$chromosome)
    chrom <-
      select(txdb,
             gene_id,
             columns = c("TXCHROM"),
             keytype = "GENEID")
    chrom <- chrom$TXCHROM[1]
    return(chrom)
  }
  
  # get strand of this gene (using txdb)
  strand <- function(gene) {
    gene_id <- sym2eg(gene)
    strand <-
      select(txdb,
             gene_id,
             columns = c("TXSTRAND"),
             keytype = "GENEID")
    strand <- strand$TXSTRAND[1]
    return(strand)
  }
  
  # get start coordinates of this gene (using txdb if possible, else using rentrez)
  getStart <- function(gene, view, strand) {
    cat(file=stderr(), "Fetching start coordinates...\n")
    gene_id <- sym2eg(gene)
    if (view == "utr") {
      txIDs <- IDs[IDs$geneID == gene_id, ]
      txID <- as.character(txIDs[1, 2])
      utrRanges <- ranges(utrs)
      gr <- utrRanges[[txID]]
      if (length(gr) > 1) {
        maxlength <- 0
        selectedUtr <- NULL
        for (i in 1:length(gr)) {
          length <- abs(end(gr[i]) - start(gr[i]))
          if (length > maxlength) {
            maxlength <- length
            selectedUtr <- gr[i]
          }
        }
        
        gr <- selectedUtr
      }
      start <- start(gr)
    } else {
      start <- start(ranges(allgenes[allgenes$gene_id == gene_id]))
      if (length(start) == 0) {
        x <- entrez_summary(db = "gene", id = gene_id)
        if (strand == "-") {
          start <- x$genomicinfo$chrstop
        } else {
          start <- x$genomicinfo$chrstart
        }
      }
    }
    cat(file=stderr(), "Start: ", start, "\n Done! \n")
    return(start)
  }
  
  # get end coordinates of this gene (using txdb if possible, else using rentrez)
  getEnd <- function(gene, view, strand) {
    cat(file=stderr(), "Fetching end coordinates...\n")
    gene_id <- sym2eg(gene)
    if (view == "utr") {
      txIDs <- IDs[IDs$geneID == gene_id, ]
      txID <- as.character(txIDs[1, 2])
      utrRanges <- ranges(utrs)
      gr <- utrRanges[[txID]]
      if (length(gr) > 1) {
        maxlength <- 0
        selectedUtr <- NULL
        for (i in 1:length(gr)) {
          length <- abs(end(gr[i]) - start(gr[i]))
          if (length > maxlength) {
            maxlength <- length
            selectedUtr <- gr[i]
          }
        }
        gr <- selectedUtr
      }
      end <- end(gr)
    } else {
      end <- end(ranges(allgenes[allgenes$gene_id == gene_id]))
      if (length(end) == 0) {
        x <- entrez_summary(db = "gene", id = gene_id)
        if (strand == "-") {
          end <- x$genomicinfo$chrstart
        } else {
          end <- x$genomicinfo$chrstop
        }
      }
    }
    cat(file=stderr(), "End: ", end, "\n Done! \n")
    return(end)
  }
  
  # creates table for this gene containing all affected conditions with dir index, p-value prox + p-value dist
  createGeneTable <- function(gene) {
    involvedConds <- getConds(gene)
    geneTable <- matrix(, nrow = 0, ncol = 3)
    colnames(geneTable) <-
      c("shortening_index", "pvalue_short", "pvalue_long")
    for (cond in involvedConds) {
      x <- rep(NA, 3)
      geneTable <- rbind(geneTable, x)
      rownames(geneTable)[rownames(geneTable) == "x"] <-
        cond
      geneTable[cond, "shortening_index"] <-
        KDmat[gene, cond]
      geneTable[cond, "pvalue_short"] <-
        format(pValPmat[gene, cond], scientific = TRUE)
      geneTable[cond, "pvalue_long"] <-
        format(pValDmat[gene, cond], scientific = TRUE)
    }
    return(geneTable)
  }
  
  # function to create ideogram track of this gene for Gviz
  ideoTrack <- function(gene, chrom) {
    itrack <-
      IdeogramTrack(genome = "hg38", chromosome = chrom)
    return(itrack)
  }
  
  # function to create transcript track of this gene for Gviz
  txTrack <- function(gene, chrom) {
    cat(file=stderr(), "Generating GeneRegionTrack...")
    txTr <-
      GeneRegionTrack(
        txdb,
        genome = "hg38",
        chromosome = chrom,
        #      symbol = gene,
        #     showId = TRUE,
        #    geneSymbol = TRUE,
        name = gene,
        # transcriptAnnotation = "symbol",
        background.title = "salmon"
      )
    displayPars(txTr) <-
      list(
        background.panel = "#ffeecc",
        col = NULL,
        size = 1.5,
        fontsize = 12
      )
    # symbols <-
    #   unlist(mapIds(org.Hs.eg.db, gene(txTr), "SYMBOL", "ENTREZID", multiVals = "first"))
    # symbol(txTr) <- symbols[gene(txTr)]
    cat(file=stderr(), "Done! \n")
    return(txTr)
  }
  
  # function to create data track of this gene (range = BigWig file, background = background color)
  createDataTrack <-
    function(gene, range, name, background, max, chrom) {
      cat(file=stderr(), "Generating data track...")
      if (name != "ctrl") {
        name = paste0(name, "_kd")
      }
      dtrack <-
        DataTrack(
          range = range,
          genome = "hg38",
          chromosome = chrom,
          type = "h",
          ylim = c(0, max),
          name = name,
          background.title = "salmon"
        )
      displayPars(dtrack) <-
        list(
          background.panel = background,
          col = NULL,
          size = 3,
          fontsize = 12
        )
      cat(file=stderr(), "Done! \n")
      return(dtrack)
    }
  
  # function to determine condition with max dir index (absolute) for this gene
  maxDI <- function(gene) {
    index <-
      which.max(abs(as.numeric(KDmat0[gene, ])))
    maxDI <- colnames(KDmat0)[index]
    return(maxDI)
  }
  
  # function to define data track for this gene, selects condition based on selection + BigWig file
  dataTrack <- function(gene, cond, range, max, chrom) {
    dtrack <-
      createDataTrack(gene, range, cond, "#ffe0cc", max, chrom)
    return(dtrack)
  }
  
  # function to define control data track for this gene
  ctrlTrack <- function(gene, range, max, chrom) {
    name <- "ctrl"
    ctrlTrack <- createDataTrack(gene, range, name, "#ccffcc", max, chrom)
    return(ctrlTrack)
  }
  
  # create poly a sites track
  polyAtrack <- function(gene, chrom, strand) {
    cat(file=stderr(), "Generating polyA track...")
    if (strand == "-") {
      range <- paste0("./data/", polyAneg)
    } else if (strand == "+") {
      range <- paste0("./data/", polyApos)
    }
    atrack <-
      AnnotationTrack(
        range = range,
        genome = "hg38",
        chromosome = chrom,
        type = "h",
        name = "APA",
        background.title = "salmon"
      )
    displayPars(atrack) <-
      list(
        background.panel = "#e6f2ff",
        col = NULL,
        size = 1.5,
        fontsize = 12
      )
    cat(file=stderr(), "Done! \n")
    return(atrack)
  }
  
  # function to create info table for gene plot (gene info + selected condition info)
  genePlotInfoTable <- function(gene, cond) {
    gi <- matrix(, nrow = 1, ncol = 10)
    strand <- strand(gene)
    colnames(gi) <-
      c(
        "Gene",
        "ID",
        "Chr",
        "Strand",
        "Start",
        "End",
        "Condition",
        "Shortening Index",
        "P-Value Short",
        "P-Value Long"
      )
    gi[1, "Gene"] <- gene
    gi[1, "ID"] <- sym2eg(gene)
    gi[1, "Chr"] <- chromosome(gene)
    gi[1, "Strand"] <- strand
    gi[1, "Start"] <- getStart(gene, "gb", strand)
    gi[1, "End"] <- getEnd(gene, "gb", strand)
    gi[1, "Condition"] <- paste0(cond, "_kd")
    gi[1, "Shortening Index"] <- KDmat0[gene, cond]
    gi[1, "P-Value Short"] <-
      format(pValPmat[gene, cond], scientific = TRUE)
    gi[1, "P-Value Long"] <-
      format(pValDmat[gene, cond], scientific = TRUE)
    return(gi)
  }
  
  # calculate the max values
  maxCovBw <- function(bw, gr) {
    ovlp <- subsetByOverlaps(bw, gr)
    if (length(ovlp) > 0) {
      max_cov <- max(ovlp$score)
    } else {
      print('WARNING: The selected genomic region has no coverage value in the BigWig')
      print('WARNING: Coverage value is arbitrary set to Zero.')
      max_cov <- 0
    }
    return(max_cov)
  }
  
  maxCovFiles <- function(bws, gr) {
    #bws <- lapply(bws, rtracklayer:::import)
    max_cov <- c()
    for (i in 1:length(gr)) {
      my_feat = gr[i,]
      max_cov[i] <- round(max(sapply(bws, maxCovBw, gr = my_feat))
                          , 2)
    }
    values(gr) <- max_cov
    return(gr)
  }
  
  # function to create list with condition and control BigWigs
  createBwList <- function(gene, cond, strand) {
    if (strand == "-") {
      if (!is.na(bw_filemat["ctrl", "neg"])) {
        ctrldata <-
          paste0("./data/modifiedBigWigs/", bw_filemat["ctrl", "neg"])
      }
      if (!is.na(bw_filemat[cond, "neg"])) {
        conddata <- paste0("./data/modifiedBigWigs/", bw_filemat[cond, "neg"])
      }
    } else if (strand == "+") {
      if (!is.na(bw_filemat["ctrl", "pos"])) {
        ctrldata <-
          paste0("./data/modifiedBigWigs/", bw_filemat["ctrl", "pos"])
      }
      if (!is.na(bw_filemat[cond, "pos"])) {
        conddata <- paste0("./data/modifiedBigWigs/", bw_filemat[cond, "pos"])
      }
    }
    bws <- list(ctrl = ctrldata, cond = conddata)
    return(bws)
  }
  
  # function to plot tracks of this gene using Gviz
  genePlot <- function(gene, cond, view) {
    chrom <- chromosome(gene)
    strand <- strand(gene)
    bws <- createBwList(gene, cond, strand)
    cat(file=stderr(), "Importing control and condition BigWigs...: \n")
    ctrldata <- import(bws$ctrl)
    cat(file=stderr(), "Control: ", bws$ctrl, "\n")
    conddata <- import(bws$cond)
    cat(file=stderr(), "Condition: ", bws$cond, "\n")
    cat(file=stderr(), "Done! \n")
    gene_id <- sym2eg(gene)
    if (view == "gb") {
      gr <- allgenes[allgenes$gene_id == gene_id]
    } else if (view == "utr") {
      txIDs <- IDs[IDs$geneID == gene_id, ]
      txID <- as.character(txIDs[1, 2])
      utrRanges <- ranges(utrs)
      gr <- utrs[[txID]]
      if (length(gr) > 1) {
        maxlength <- 0
        selectedUtr <- NULL
        for (i in 1:length(gr)) {
          length <- abs(end(gr[i]) - start(gr[i]))
          if (length > maxlength) {
            maxlength <- length
            selectedUtr <- gr[i]
          }
        }
        gr <- selectedUtr
      }
    }
    cat(file=stderr(), "Calculating max coverage... \n")
    maxCoverage <- as.data.frame(maxCovFiles(list(ctrldata, conddata), gr))[1,6]
    cat(file=stderr(), "Max coverage: ", maxCoverage, "\n")
    cat(file=stderr(), "Done! \n")
    cat(file=stderr(), "Plotting Gviz tracks... \n")
    plotTracks(
      list(
        #     ideoTrack(gene, chrom)
        gtrack,
        txTrack(gene, chrom),
        polyAtrack(gene, chrom, strand),
        dataTrack(gene, cond, bws$cond, maxCoverage, chrom),
        ctrlTrack(gene, bws$ctrl, maxCoverage, chrom)
      ),
      from = getStart(gene, view, strand),
      to = getEnd(gene, view, strand),
      extend.left = 0.5,
      extend.right = 0.5,
      showBandId = TRUE,
      add53 = TRUE,
      add35 = TRUE,
      title.width = 1.8
    )
    # }, height = function() {
    # session$clientData$output_gviz_width
  }
  
  # creates info table for this gene (name, ID, description, strand, chromosome, coords)
  geneInfoTable <- function(gene) {
    gene_id <- sym2eg(gene)
    gi <- matrix(, nrow = 1, ncol = 7)
    strand <- strand(gene)
    colnames(gi) <-
      c("Gene",
        "ID",
        "Description",
        "Chr",
        "Strand",
        "Start",
        "End")
    desc <- ""
    try(x <- entrez_summary(db = "gene", id = gene_id))
    try(desc <- x$description)
    gi[1, "Gene"] <- gene
    gi[1, "ID"] <- gene_id
    gi[1, "Description"] <- desc
    gi[1, "Chr"] <- chromosome(gene)
    gi[1, "Strand"] <- strand
    gi[1, "Start"] <- getStart(gene, "gb", strand)
    gi[1, "End"] <- getEnd(gene, "gb", strand)
    return(gi)
  }
  
  # function to create URL for genome browser
  createBrowserURL <- function(gene) {
    polyA <- ""
    miRNA <- ""
    ctrl <- ""
    if ('1' %in% input$selectTracks) {
      miRNA <- "&hub_114483_miRNA=full"
    } else {
      miRNA <- "&hub_114483_miRNA=hide"
    }
    
    if (nchar(gene) > 1) {
      strand <- strand(gene)
      if (input$condsGB == "max SI") {
        cond <- maxDI(input$geneInput)
      } else {
        cond <- strsplit(input$condsGB, "_")[[1]][1]
      }
      if (strand == "-") {
        if (!is.na(bw_filemat[cond, "neg"])) {
          file <- removeExt(bw_filemat[cond, "neg"], sep = ".")
        }
        polyA <-
          "&hub_114483_adj_contigs_neg=dense"
        ctrl <-
          "&hub_114483_UTR2_sum_ctrl_neg=full"
        
      } else if (strand == "+") {
        if (!is.na(bw_filemat[cond, "pos"])) {
          file <- removeExt(bw_filemat[cond, "pos"], sep = ".")
        }
        polyA <-
          "&hub_114483_adj_contigs_pos=dense"
        ctrl <-
          "&hub_114483_UTR2_sum_ctrl_pos=full"
      }
      pos <-
        paste0(chromosome(gene),
               ":",
               getStart(gene, "utr", strand),
               "-",
               getEnd(gene, "utr", strand))
      
      browserurl <<-
        paste0(
          "http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&hubUrl=ftp://ftp.imbei.uni-mainz.de/trendseqHub/hub.txt&position=",
          pos,
          "&hideTracks=1&knownGene=pack&refGene=pack&omimAvSnp=pack&all_mrna=dense&wgEncodeReg=show&snp147Common=dense&rmsk=full&hub_114483_",
          file,
          "=full",
          ctrl,
          polyA,
          miRNA
        )
    } else {
      browserurl <<-
        paste0(
          "http://genome.ucsc.edu/cgi-bin/hgTracks?db=hg38&hubUrl=ftp://ftp.imbei.uni-mainz.de/trendseqHub/hub.txt",
          polyA,
          miRNA
        )
    }
    return(browserurl)
  }
  
  # function to check if plot is already cached
  isCached <- function(id) {
    file <- paste0("./data/plotCache/", id, ".png")
    if (file.exists(file)) {
      return(TRUE)
    } else {
      return(FALSE)
    }
  }
  
  # function to create and cache plot
  cachePlot <- function(file, gene, cond, view) {
    try(genePlot(gene, cond, view))
    dev.print(png, file, width = 800)
    # dev.off()
  }
}


# Launching the app -------------------------------------------------------

shinyApp(ui = ui, server = server)



