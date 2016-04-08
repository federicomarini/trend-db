## ------------------------------------------------------------------ ##
##                          KD Overview                               ##
## ------------------------------------------------------------------ ##

library(tools)

setwd("/Users/Denise/Documents/Master/IMBEI/data/")

# read KD data
KDdata <- read.csv("KDmat.csv", header = FALSE, sep = ";")

KDmat <- as.matrix(KDdata)

head(KDmat)
dim(KDmat)

# save conditions as colnames and genes as rownames
colnames(KDmat) <- KDmat[1,]
rownames(KDmat) <- KDmat[, 1]

KDmat <- KDmat[-1,-1]

# create overview matrix
KDOverview <- KDmat
KDOverview[KDOverview > 0] <- 1
KDOverview[KDOverview < 0] <- 1
KDOverview[is.na(KDOverview)] <- 0

# list with all conditions/KDs
conditions <- sort(colnames(KDmat))

# list with genes
genenames <- rownames(KDmat)

## ------------------------------------------------------------------ ##
##                          BigWigs                                   ##
## ------------------------------------------------------------------ ##

setwd("/Users/Denise/Documents/Master/IMBEI/data/IGV_vis_drops")

# get all bigwig files
bigwigs = list.files(pattern=".bw")

# create matrix 
bw_filemat <- matrix(,nrow = 0, ncol = 2)

# set colnames
colnames(bw_filemat) <- c("pos", "neg")

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
  if(!(gene %in% rownames(bw_filemat))) {
    bw_filemat <- rbind(bw_filemat,x)
    rownames(bw_filemat)[rownames(bw_filemat)=="x"] <- gene
  }
  
  # add file for pos strand
  if(strand == "pos") {
    # if UTR3 exists, do not use UTR2 file
    if (utr == "UTR2" && !is.na(bw_filemat[gene,"pos"])) {}
    else bw_filemat[gene,"pos"] <- file
  }
  
  # add file for neg strand
  if(strand == "neg") {
    # if UTR3 exists, do not use UTR2 file
    if (utr == "UTR2" && !is.na(bw_filemat[gene,"neg"])) {}
    else bw_filemat[gene,"neg"] <- file
  }
}

# sort matrix by gene names
bw_filemat <- bw_filemat[order(rownames(bw_filemat)), ] 


## ------------------------------------------------------------------ ##
##                          Gviz                                      ##
## ------------------------------------------------------------------ ##

library(Gviz)
library(org.Hs.eg.db)

#hg38db <- makeTxDbFromUCSC(genome = "hg38", tablename ="knownGene")
#saveDb(hg38db, file="hg38db_knownGene.sqlite")

# load annotation db
txdb <-loadDb("hg38db_knownGene.sqlite")

#txdb <- makeTxDbFromGFF("Homo_sapiens.GRCh38.84.gtf", format="gtf")

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

library("shiny")
library("DT")
library(shinydashboard)

ui <-
  shinydashboard::dashboardPage(
    dashboardHeader(
      title = paste0("Prototype of trendseq exploration app"),
      titleWidth = 900
    ),
    
    dashboardSidebar(
      width = 350,
      menuItem(
        "App settings",
        icon = icon("cogs"),
        textInput("inText1", "Type something")
        
      ),
      
      menuItem(
        "Plot settings",
        icon = icon("paint-brush"),
        # bstooltip: Use the widgets below to setup general parameters for exporting produced plots
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
        width = 12,
        
        tabPanel(
          "Data Preview",
          h1(
            "Here will go some head of the count dataset, the samples design/covariates and so"
          ),
          DT::dataTableOutput("inspectMatrix")
        ),
        
        tabPanel(
          "Main View",
          p(h1('MAIN VIEW')),
          selectInput(
            "geneInput",
            "Select gene",
            choices = c("", genenames),
            selected = NULL
          ),
          selectInput(
            "kdInput",
            "Select knockdown",
            choices = c("", conditions),
            selected = NULL
          ),
          
          h2("AFFECTED GENES"),
          textOutput("printAffectedGenes"),
          
          h2("INVOLVED CONDITIONS"),
          textOutput("printConditionsInvolved")
          
        ),
        
        tabPanel(
          "Condition Overview",
          h1("Overview matrix of the conditions, genes with 1 are affected"),
          DT::dataTableOutput("overviewMatrix")
        ),
   
      tabPanel("Gviz Gene View",
               h1("Gene Plot"),
               selectInput(
                 "genePlotInput",
                 "Select gene",
                 choices = c("", genenames),
                 selected = NULL
               ), 
               selectInput(
                 "strandInput", 
                 "Select strand", 
                 choices = c(""," pos", "neg"), 
                 selected = NULL
              ), 
              submitButton(
                "Update View", 
                icon("refresh")
              ),
              plotOutput("gviz"))
      
    ),
    
    skin = "blue"
  ))

## ------------------------------------------------------------------ ##
##                          Define server                             ##
## ------------------------------------------------------------------ ##

server <- function(input, output) {
  
  output$inspectMatrix <- renderDataTable({
    table <- head(KDmat, 20)
    DT::datatable(table, options = list(scrollX = TRUE))
  })
  
  output$overviewMatrix <- renderDataTable({
    table <- KDOverview
    DT::datatable(table, options = list(scrollX = TRUE, pageLength = 15))
  })
  
  affectedGenes <- reactive({
    allgenes_selectedKD <-
      KDmat[, colnames(KDmat) %in% input$kdInput, drop = FALSE] # force staying as a matrix
    affectedGenes <-
      rownames(allgenes_selectedKD)[!is.na(allgenes_selectedKD)]
    return(affectedGenes)
  })
  
  output$printAffectedGenes <- renderPrint({
    cat(
      paste0(
        "The following genes are affected in the ",
        input$kdInput,
        " condition:\n",
        paste(affectedGenes(), collapse = ", ")
      )
    )
    
  })
  
  conditionsInvolved <- reactive({
    allconds_selectedGene <-
      KDmat[(rownames(KDmat) %in% input$geneInput), , drop = FALSE] # force staying as a matrix
    involvedConds <-
      colnames(allconds_selectedGene)[!is.na(allconds_selectedGene)]
    return(involvedConds)
  })
  
  
  output$printConditionsInvolved <- renderPrint({
    cat(
      paste0(
        "The following conditions are involved in the regulation of the ",
        input$geneInput,
        " gene:\n",
        paste(conditionsInvolved(), collapse = ", ")
      )
    )
    
  })
  
  chromosome <- reactive({
    if(length(sym2eg_list) > 0) {
      gene_id <- sym2eg_list[input$genePlotInput]
      gene_id <- as.character(gene_id[1])
    }
      chrom <- select(txdb, gene_id, columns=c("TXCHROM"), keytype = "GENEID")
      chrom <- chrom$TXCHROM[1]
      return(chrom)
  })
  
  ideoTrack <- reactive({
    chr <- chromosome() 
    itrack <- IdeogramTrack(genome = "hg38", chromosome = chr)
    return(itrack)
  })
  
  getStart <- reactive({
    if(length(sym2eg_list) > 0) {
      gene_id <- sym2eg_list[input$genePlotInput]
      gene_id <- as.character(gene_id[1])
    }
    start <- start(ranges(allgenes[allgenes$gene_id == gene_id]))
    return(start)
  })
  
  getEnd <- reactive({
    if(length(sym2eg_list) > 0) {
      gene_id <- sym2eg_list[input$genePlotInput]
      gene_id <- as.character(gene_id[1])
    }
    end <- end(ranges(allgenes[allgenes$gene_id == gene_id]))
    return(end)
  })
  
  txTrack <- reactive({
    txTr <- GeneRegionTrack(txdb, genome = "hg38", chromosome = chromosome(), symbol = input$genePlotInput, showId = TRUE, geneSymbol = TRUE, name = "UCSC", background.title = "salmon")
    displayPars(txTr) <- list(background.panel = "#FFFEDB",col = NULL)
    symbols <- unlist(mapIds(org.Hs.eg.db, gene(txTr), "SYMBOL", "ENTREZID", multiVals = "first"))
    symbol(txTr) <- symbols[gene(txTr)]
    return(txTr)
  })
  
  dataTracks <- reactive({
    allconds <- KDmat[(rownames(KDmat) %in% input$genePlotInput), , drop = FALSE] # force staying as a matrix
    involved <- colnames(allconds)[!is.na(allconds)]
    
    dTracks <- list()
    
    for (con in involved) {
      
        if(input$strandInput == "neg") {
        range <- bw_filemat[con,"neg"]
        cat(file=stderr(), "Using range: ", range, "\n")
        if (!(is.na(range))) {
          dtrack_neg <- DataTrack(range = range, genome = "hg38", chromosome = chromosome(), type = "h", name = con, background.title = "salmon")
          displayPars(dtrack_neg) <- list(background.panel = "#FFFEDB",col = NULL)
          dTracks[[length(dTracks)+1]] <- dtrack_neg
        }
        }
      
      if(input$strandInput == "pos") {
        range <- bw_filemat[con,"pos"]
        cat(file=stderr(), "Using range: ", range, "\n")
        if (!(is.na(range))) {
          dtrack_pos <- DataTrack(range = range, genome = "hg38", chromosome = chromosome(), type = "h", name = con, background.title = "salmon")
          displayPars(dtrack_pos) <- list(background.panel = "#FFFEDB",col = NULL)
          dTracks[[length(dTracks)+1]] <- dtrack_pos
        }
      }
     }
   return(dTracks)
  })
  
  output$gviz <- renderPlot({
    plotTracks(append(list(ideoTrack(), gtrack, txTrack()), dataTracks()), from = getStart(), to = getEnd(), extend.left = -0.2, extend.right = -0.2, showBandId = TRUE, add53 = TRUE, add35 = TRUE)
  })
  
  
}

shinyApp(ui = ui, server = server)


