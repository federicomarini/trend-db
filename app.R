## ------------------------------------------------------------------ ##
##                          KD Overview                               ##
## ------------------------------------------------------------------ ##

library(tools)

setwd("/Users/Denise/Documents/Master/IMBEI/data/")

# read KD data
KDdata <-
  read.csv(
    "KDmat.csv",
    dec = ".",
    sep = ",",
    header = FALSE,
    stringsAsFactors = FALSE,
    na.strings = "NA"
  )
pValP <-
  read.csv(
    "pProx.csv",
    dec = ".",
    sep = ",",
    header = FALSE,
    stringsAsFactors = FALSE,
    na.strings = "NA"
  )
pValD <-
  read.csv(
    "pDist.csv",
    dec = ".",
    sep = ",",
    header = FALSE,
    stringsAsFactors = FALSE,
    na.strings = "NA"
  )

KDmat <- as.matrix(KDdata)
pValPmat <- as.matrix(pValP)
pValDmat <- as.matrix(pValD)

# save conditions as colnames and genes as rownames
colnames(KDmat) <- KDmat[1,]
rownames(KDmat) <- KDmat[, 1]
colnames(pValPmat) <- pValPmat[1,]
rownames(pValPmat) <- pValPmat[, 1]
colnames(pValDmat) <- pValDmat[1,]
rownames(pValDmat) <- pValDmat[, 1]

KDmat <- KDmat[-1,-1]
pValPmat <- pValPmat[-1,-1]
pValDmat <- pValDmat[-1,-1]

class(KDmat) <- "numeric"
class(pValPmat) <- "numeric"
class(pValDmat) <- "numeric"

# create overview matrix
KDOverview <- KDmat
KDOverview[KDOverview > 0] <- 1
KDOverview[KDOverview < 0] <- 1
KDOverview[is.na(KDOverview)] <- 0

# create matrix with 0s instead of NAs
KDmat0 <- KDmat
KDmat0[is.na(KDmat0)] <- 0
KDmat0 <- data.matrix(KDmat0)

# list with all conditions/KDs
conditions <- sort(colnames(KDmat))

# list with genes
genenames <- rownames(KDmat)

cat(file = stderr(), str(KDmat0), "\n")

## ------------------------------------------------------------------ ##
##                          BigWigs                                   ##
## ------------------------------------------------------------------ ##

setwd("/Users/Denise/Documents/Master/IMBEI/data/IGV_vis_drops")

# get all bigwig files
bigwigs = list.files(pattern = ".bw")

# create matrix
bw_filemat <- matrix(, nrow = 0, ncol = 2)

# set colnames
colnames(bw_filemat) <- c("pos", "neg")

for (con in conditions) {
  x <- character(2)
  x <- NA
  
  bw_filemat <- rbind(bw_filemat, x)
  rownames(bw_filemat)[rownames(bw_filemat) == "x"] <- con
  
}

for (file in bigwigs) {
  x <- character(2)
  x <- NA
  
  filename = file_path_sans_ext(file)
  str = strsplit(filename, '_')
  # UTR2 / UTR3
  utr = str[[1]][1]
  # gene name
  gene = str[[1]][3]
  # strand (pos / neg)
  strand = str[[1]][4]
  
  # create row for each gene
  if (!(gene %in% rownames(bw_filemat))) {
    bw_filemat <- rbind(bw_filemat, x)
    rownames(bw_filemat)[rownames(bw_filemat) == "x"] <- gene
  }
  
  # add file for pos strand
  if (strand == "pos") {
    # if UTR3 exists, do not use UTR2 file
    if (utr == "UTR2" && !is.na(bw_filemat[gene, "pos"])) {
      
    }
    else
      bw_filemat[gene, "pos"] <- file
  }
  
  # add file for neg strand
  if (strand == "neg") {
    # if UTR3 exists, do not use UTR2 file
    if (utr == "UTR2" && !is.na(bw_filemat[gene, "neg"])) {
      
    }
    else
      bw_filemat[gene, "neg"] <- file
  }
}

# sort matrix by gene names
bw_filemat <- bw_filemat[order(rownames(bw_filemat)), ]


## ------------------------------------------------------------------ ##
##                          Gviz                                      ##
## ------------------------------------------------------------------ ##

library(Gviz)
library(org.Hs.eg.db)

setwd("/Users/Denise/Documents/Master/IMBEI/data/")

#hg38db <- makeTxDbFromUCSC(genome = "hg38", tablename ="knownGene")
#saveDb(hg38db, file="hg38db_knownGene.sqlite")

# load annotation db
txdb <- loadDb("hg38db_refGene.sqlite")

# extract all genes as GRanges object
allgenes <- genes(txdb)

# get the entrez gene identifiers that are mapped to a gene symbol & save as list
sym2eg <- org.Hs.egSYMBOL2EG
mapped_genes <- mappedkeys(sym2eg)
sym2eg_list <- as.list(sym2eg[mapped_genes])

# create genome axis track
gtrack <- GenomeAxisTrack()

## ------------------------------------------------------------------ ##
##                          Shiny                                     ##
## ------------------------------------------------------------------ ##

library(shiny)
library(DT)
library(shinydashboard)

ui <-
  shinydashboard::dashboardPage(
    dashboardHeader(
      title = paste0("Prototype of trendseq exploration app"),
      titleWidth = 900
    ),
    dashboardSidebar(
      disable = TRUE,
      width = 350,
      menuItem(
        "App settings",
        icon = icon("cogs"),
        textInput("inText1", "Type something")
      ),
      menuItem(
        "Plot settings",
        icon = icon("paint-brush"),
        numericInput(
          "export_width",
          label = "Width of exported figures (cm)",
          value = 30,
          min = 2
        ),
        numericInput(
          "export_height",
          label = "Height of exported figures (cm)",
          value = 30,
          min = 2
        )
      )
    ),
    dashboardBody(
      tabBox(
        width = 11,
        tabPanel("About",
                 h1("About"),
                 textOutput("about")),
        #h4("Session Info"),
        #verbatimTextOutput("sessioninfo")),
        tabPanel("Data Preview",
                 h1("Inspect Matrix"),
                 fluidRow(
                   column(
                     width = 9,
                     sliderInput(
                       "dirFilter",
                       "Filter by dir index:",
                       min = 0,
                       max = 15,
                       value = 0,
                       step = 0.5
                     ),
                     DT::dataTableOutput("inspectMatrix"),
                     br(),
                     textOutput("selection"),
                     DT::dataTableOutput("overviewTable")
                   )
                 )),
        tabPanel(
          "Condition Overview",
          h1("Overview matrix of the conditions, genes with 1 are affected"),
          fluidRow(column(
            width = 9,
            DT::dataTableOutput("overviewMatrix")
          ))
        ),
        tabPanel(
          "Main View",
          p(h1('Main View')),
          br(),
          fluidRow(column(
            width = 9,
            selectInput(
              "geneInput",
              "Select gene",
              choices = c("", genenames),
              selected = NULL
            ),
            h4("Involved Conditions:"),
            #textOutput("printConditionsInvolved"),
            DT::dataTableOutput("geneTableOutput")
          )),
          br(),
          fluidRow(column(
            width = 9,
            selectInput(
              "kdInput",
              "Select knockdown",
              choices = c("", conditions),
              selected = NULL
            ),
            h4("Affected genes:"),
            #textOutput("printAffectedGenes")
            DT::dataTableOutput("kdTableOutput")
          ))
        ),
        tabPanel(
          "Gviz Gene Plot",
          h1("Gene Plot"),
          fluidRow(
            column(
              width = 5,
              selectInput(
                "genePlotInput",
                "Select gene",
                choices = c("", genenames),
                selected = NULL
              ),
              br(),
              actionButton("plotSubmit", "Create Plot")
            ),
            column(
              width = 5,
              sliderInput(
                "extendUpstream",
                "Extend view upstream:",
                min = 0,
                max = 5000,
                value = 0,
                step = 100
              ),
              sliderInput(
                "extendDownstream",
                "Extend view downstream:",
                min = 0,
                max = 5000,
                value = 0,
                step = 100
              )
              #numericInput("extendUpstream", "Extend view upstream", value = 0),
              #numericInput("extendDownstream", "Extend view downstream", value = 0),
            )
          ),
          br(),
          br(),
          DT::dataTableOutput("geneInfo"),
          br(),
          htmlOutput("NCBIref"),
          br(),
          br(),
          checkboxInput('returnpdf', 'Save as PDF?', FALSE),
          conditionalPanel(condition = "input.returnpdf == true",
                           downloadLink('pdflink')),
          br(),
          plotOutput("gviz")
          # , height ="auto")
          
        )
      )
    ),
    skin = "blue"
  )

## ------------------------------------------------------------------ ##
##                          Define server                             ##
## ------------------------------------------------------------------ ##

server <- function(input, output, session) {
  output$inspectMatrix <- DT::renderDataTable(server = TRUE, {
    table <- KDmat
    colnames(table) <- paste0(colnames(KDmat), "_si")
    rownames(table) <- rownames(KDmat)
    table <- as.data.frame(table)
    cutoff <- input$dirFilter
    table <- subset(table, apply(table, MARGIN = 1, function(x) any(x > cutoff | x < -cutoff)))
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
        fixedColumns = TRUE
      )
    ) %>%
      formatRound(colnames(table), digits = 3) %>%
      formatStyle("Gene",
                  backgroundColor = "#ffeecc",
                  fontWeight = "bold")
  })
  
  proxy = dataTableProxy("inspectMatrix")
  
  output$selection <- renderPrint({
    req(input$inspectMatrix_rows_selected)
    cutoff <- input$dirFilter
    table <- as.data.frame(KDmat)
    if (input$dirFilter > 0) {
      table <- subset(table, apply(table, MARGIN = 1, function(x) any(x > cutoff | x < -cutoff)))
    }
    selectedRowIndex <- input$inspectMatrix_rows_selected
    selectedRow <- table[selectedRowIndex,]
    gene <- rownames(table)[selectedRowIndex]
    cat("Selected Gene: ", gene)
  })
  
  output$overviewTable <- DT::renderDataTable({
    table <- overviewTable()
    DT::datatable(
      table,
      colnames = c("Condition" = 1),
      rownames = paste0(rownames(table), "_si"),
      options = list(
        scrollY = 300,
        scrollCollapse = TRUE,
        paging = FALSE,
        searching = FALSE
      )
    ) %>%
      formatStyle("Condition",
                  backgroundColor = "#ffe0cc",
                  fontWeight = "bold")
  })
  
  overviewTable <-
    eventReactive(input$inspectMatrix_rows_selected, {
      cutoff <- input$dirFilter
      table <- as.data.frame(KDmat)
      if (input$dirFilter > 0) {
        table <- subset(table, apply(table, MARGIN = 1, function(x) any(x > cutoff | x < -cutoff)))
      }
      selectedRowIndex <- input$inspectMatrix_rows_selected
      selectedRow <- table[selectedRowIndex,]
      gene <- rownames(table)[selectedRowIndex]
      overviewTable <- createGeneTable(gene)
      return(overviewTable)
    })
  
  output$overviewMatrix <- DT::renderDataTable(server = TRUE, {
    table <- KDOverview
    colnames(table) <- paste0(colnames(KDOverview), "_si")
    DT::datatable(
      table,
      colnames = c("Gene" = 1),
      extensions = c('FixedColumns', 'Scroller'),
      options = list(
        scrollX = TRUE,
        deferRender = TRUE,
        scrollY = 407,
        scrollCollapse = TRUE,
        scroller = TRUE,
        pageLength = 50,
        fixedColumns = TRUE
      )
    )  %>%
      formatStyle("Gene",
                  backgroundColor = "#ffeecc",
                  fontWeight = "bold")  %>%
      formatStyle(
        colnames(table),
        fontWeight = styleInterval(0, c('normal', 'bold')),
        backgroundColor = styleInterval(0, c('normal', '#b3ffb3')))
  })
  
  getGenes <- function(cond) {
    allgenes_selectedKD <-
      KDmat[, colnames(KDmat) %in% cond, drop = FALSE] # force staying as a matrix
    affectedGenes <-
      rownames(allgenes_selectedKD)[!is.na(allgenes_selectedKD)]
    return(affectedGenes)
  }
  
  affectedGenes <- reactive({
    getGenes(input$kdInput)
  })
  
  output$printAffectedGenes <- renderPrint({
    cat(paste(affectedGenes(), collapse = ", "))
  })
  
  getConds <- function(gene) {
    allconds_selectedGene <-
      KDmat[(rownames(KDmat) %in% gene), , drop = FALSE] # force staying as a matrix
    involvedConds <-
      colnames(allconds_selectedGene)[!is.na(allconds_selectedGene)]
    return(involvedConds)
  }
  
  conditionsInvolved <- reactive({
    getConds(input$geneInput)
  })
  
  output$printConditionsInvolved <- renderPrint({
    cat(paste(conditionsInvolved(), collapse = ", ", "\n"))
  })
  
  sym2eg <- reactive({
    if (length(sym2eg_list) > 0) {
      gene_id <- sym2eg_list[input$genePlotInput]
      gene_id <- as.character(gene_id[1])
    }
    return(gene_id)
  })
  
  chromosome <- reactive({
    gene_id <- sym2eg()
    chrom <-
      select(txdb,
             gene_id,
             columns = c("TXCHROM"),
             keytype = "GENEID")
    chrom <- chrom$TXCHROM[1]
    cat(file = stderr(), "Chromosome: ", chrom, "\n")
    return(chrom)
  })
  
  strand <- reactive({
    gene_id <- sym2eg()
    strand <-
      select(txdb,
             gene_id,
             columns = c("TXSTRAND"),
             keytype = "GENEID")
    strand <- strand$TXSTRAND[1]
    cat(
      file = stderr(),
      "Gene: ",
      input$genePlotInput,
      " Gene ID: ",
      gene_id,
      " Strand: ",
      strand,
      "\n"
    )
    return(strand)
  })
  
  getStart <- reactive({
    gene_id <- sym2eg()
    start <-
      start(ranges(allgenes[allgenes$gene_id == gene_id]))
    cat(file = stderr(), "Gene start: ", start, "\n")
    if (is.na(input$extendUpstream)) {
      return(start)
    } else
      return(start - input$extendUpstream)
  })
  
  getEnd <- reactive({
    gene_id <- sym2eg()
    end <-
      end(ranges(allgenes[allgenes$gene_id == gene_id]))
    cat(file = stderr(), "Gene end: ", end, "\n")
    if (is.na(input$extendDownstream)) {
      return(end)
    } else
      return(end + input$extendDownstream)
  })
  
  createGeneTable <- function(gene) {
    involvedConds <- getConds(gene)
    geneTable <- matrix(, nrow = 0, ncol = 3)
    colnames(geneTable) <-
      c("dir_index", "pvalue_prox", "pvalue_dist")
    for (cond in involvedConds) {
      cat(file = stderr(), cond, "\n")
      x <- rep(NA, 3)
      geneTable <- rbind(geneTable, x)
      rownames(geneTable)[rownames(geneTable) == "x"] <- cond
      geneTable[cond, "dir_index"] <-
        KDmat[gene, cond]
      geneTable[cond, "pvalue_prox"] <-
        format(pValPmat[gene, cond], scientific = TRUE)
      geneTable[cond, "pvalue_dist"] <-
        format(pValDmat[gene, cond], scientific = TRUE)
    }
    return(geneTable)
  }
  
  geneTable <- reactive({
    createGeneTable(input$geneInput)
  })
  
  output$geneTableOutput <- DT::renderDataTable({
    req(input$geneInput)
    table <- geneTable()
    DT::datatable(
      table,
      colnames = c("Condition" = 1),
      rownames = paste0(rownames(table), "_si"),
      options = list(
        scrollY = 300,
        scrollCollapse = TRUE,
        paging = FALSE,
        searching = FALSE
      )
    ) %>%
      formatStyle("Condition",
                  backgroundColor = "#ffe0cc",
                  fontWeight = "bold")
  })
  
  kdTable <- reactive({
    affectedGenes <- affectedGenes()
    kdTable <- matrix(, nrow = 0, ncol = 3)
    colnames(kdTable) <-
      c("dir_index", "pvalue_prox", "pvalue_dist")
    for (gene in affectedGenes) {
      cat(file = stderr(), gene, "\n")
      x <- rep(NA, 3)
      kdTable <- rbind(kdTable, x)
      rownames(kdTable)[rownames(kdTable) == "x"] <- gene
      kdTable[gene, "dir_index"] <-
        KDmat[gene, input$kdInput]
      kdTable[gene, "pvalue_prox"] <-
        format(pValPmat[gene, input$kdInput], scientific = TRUE)
      kdTable[gene, "pvalue_dist"] <-
        format(pValDmat[gene, input$kdInput], scientific = TRUE)
    }
    return(kdTable)
  })
  
  output$kdTableOutput <- DT::renderDataTable({
    req(input$kdInput)
    DT::datatable(
      kdTable(),
      colnames = c("Gene" = 1),
      options = list(
        scrollY = 300,
        scrollCollapse = TRUE,
        paging = FALSE,
        searching = FALSE
      )
    )  %>%
      formatStyle("Gene",
                  backgroundColor = "#ffeecc",
                  fontWeight = "bold")
  })
  
  ideoTrack <- reactive({
    chr <- chromosome()
    itrack <-
      IdeogramTrack(genome = "hg38", chromosome = chr)
    return(itrack)
  })
  
  txTrack <- reactive({
    txTr <-
      GeneRegionTrack(
        txdb,
        genome = "hg38",
        chromosome = chromosome(),
        symbol = input$genePlotInput,
        showId = TRUE,
        geneSymbol = TRUE,
        name = "UCSC",
        background.title = "salmon"
      )
    displayPars(txTr) <-
      list(background.panel = "#ffeecc", col = NULL)
    symbols <-
      unlist(mapIds(org.Hs.eg.db, gene(txTr), "SYMBOL", "ENTREZID", multiVals = "first"))
    symbol(txTr) <- symbols[gene(txTr)]
    return(txTr)
  })
  
  createDataTrack <- function(range, name, background) {
    dtrack <-
      DataTrack(
        range = range,
        genome = "hg38",
        chromosome = chromosome(),
        type = "h",
        name = name,
        background.title = "salmon"
      )
    displayPars(dtrack) <-
      list(background.panel = background, col = NULL)
    return(dtrack)
  }
  
  maxDI <- reactive({
    index <-
      which.max(abs(as.numeric(KDmat0[input$genePlotInput,])))
    # cat(file = stderr(), "Max DI Index: ", index, "\n")
    maxDI <- colnames(KDmat0)[index]
    cat(file = stderr(), "Max DI: ", maxDI, "\n")
    return(maxDI)
  })
  
  dataTrackMaxDI <- reactive({
    setwd("/Users/Denise/Documents/Master/IMBEI/data/IGV_vis_drops")
    maxDI <- maxDI()
    if (strand() == "-") {
      if (!is.na(bw_filemat[maxDI, "neg"])) {
        range <- bw_filemat[maxDI, "neg"]
        dtrack <- createDataTrack(range, maxDI, "#ffe0cc")
      }
    } else if (strand() == "+") {
      if (!is.na(bw_filemat[maxDI, "pos"])) {
        range <- bw_filemat[maxDI, "pos"]
        dtrack <- createDataTrack(range, maxDI, "#ffe0cc")
      }
    }
    return(dtrack)
  })
  
  ctrlTrack <- reactive({
    setwd("/Users/Denise/Documents/Master/IMBEI/data/IGV_vis_drops")
    name <- "ctrl"
    if (strand() == "-") {
      if (!is.na(bw_filemat["ctrl", "neg"])) {
        range <- bw_filemat["ctrl", "neg"]
        ctrlTrack <-
          createDataTrack(range, name, "#ccffcc")
      }
    } else if (strand() == "+") {
      if (!is.na(bw_filemat["ctrl", "pos"])) {
        range <- bw_filemat["ctrl", "pos"]
        ctrlTrack <-
          createDataTrack(range, name, "#ccffcc")
      }
    }
    return(ctrlTrack)
  })
  
  genePlot <- reactive({
    req(input$genePlotInput)
    if (input$returnpdf) {
      pdf("plot.pdf")
      plotTracks(
        list(
          ideoTrack(),
          gtrack,
          txTrack(),
          dataTrackMaxDI(),
          ctrlTrack()
        ),
        from = getStart(),
        to = getEnd(),
        extend.left = 50,
        showBandId = TRUE,
        add53 = TRUE,
        add35 = TRUE,
        sizes = NULL
      )
      dev.off()
    }
    plotTracks(
      list(
        ideoTrack(),
        gtrack,
        txTrack(),
        dataTrackMaxDI(),
        ctrlTrack()
      ),
      from = getStart(),
      to = getEnd(),
      extend.left = 50,
      showBandId = TRUE,
      add53 = TRUE,
      add35 = TRUE,
      sizes = NULL
    )
    # }, height = function() {
    # session$clientData$output_gviz_width
  })
  
  drawPlot <- eventReactive(input$plotSubmit, {
    genePlot()
  })
  
  output$gviz <- renderPlot({
    drawPlot()
  })
  
  geneInfoTable <- reactive({
    gi <- matrix(, nrow = 1, ncol = 10)
    colnames(gi) <-
      c(
        "Gene",
        "ID",
        "Chr",
        "Strand",
        "Start",
        "End",
        "Condition",
        "Dir Index",
        "P-Value Prox",
        "P-Value Dist"
      )
    cond <- maxDI()
    gi[1, "Gene"] <- input$genePlotInput
    gi[1, "ID"] <- sym2eg()
    gi[1, "Chr"] <- chromosome()
    gi[1, "Strand"] <- strand()
    gi[1, "Start"] <- getStart()
    gi[1, "End"] <- getEnd()
    gi[1, "Condition"] <- paste0(cond, "_si")
    gi[1, "Dir Index"] <- KDmat0[input$genePlotInput, cond]
    gi[1, "P-Value Prox"] <-
      format(pValPmat[input$genePlotInput, cond], scientific = TRUE)
    gi[1, "P-Value Dist"] <-
      format(pValDmat[input$genePlotInput, cond], scientific = TRUE)
    return(gi)
  })
  
  output$geneInfo <- DT::renderDataTable({
    req(input$genePlotInput)
    DT::datatable(geneInfoTable(), options = list(paging = FALSE, searching = FALSE))  %>%
      formatStyle("Gene",
                  backgroundColor = "#ffeecc",
                  fontWeight = "bold") %>%
      formatStyle("Condition",
                  backgroundColor = "#ffe0cc",
                  fontWeight = "bold")
  })
  
  NCBIref <- reactive({
    req(input$genePlotInput)
    return(paste0("http://www.ncbi.nlm.nih.gov/gene/", sym2eg()))
  })
  
  output$NCBIref <- renderUI({
    a("View Gene (NCBI)", href = NCBIref(), target = "_blank")
  })
  
  output$pdflink <- downloadHandler(filename <- "myplot.pdf",
                                    content <- function(file) {
                                      file.copy("plot.pdf", file)
                                    })
  
  output$about <- renderPrint({
    cat("Prototype of the trendseq exploration app")
  })
  
  output$sessioninfo <- renderPrint({
    sessionInfo()
  })
}

shinyApp(ui = ui, server = server)
