library("shiny")
load("../F11-juliane/annoHuman_rel82.RData")


set.seed(42)
allcoding <- annoHuman[annoHuman$gene_biotype=="protein_coding",]
genenames <- allcoding$gene_name %>% sort %>% unique()

alltheknockdowns <- paste0("KD_",1:20,"_",sample(genenames,20))

proto_kdmat <- matrix(0,length(genenames),length(alltheknockdowns))
dim(proto_kdmat)

for (i in 1:length(alltheknockdowns)){
  proto_kdmat[sample((1:nrow(proto_kdmat)),200),i] <- 1 # set 200 random genes to 1
}

head(proto_kdmat)
rownames(proto_kdmat) <- genenames
colnames(proto_kdmat) <- alltheknockdowns

library("DT")
library(shinydashboard)

protoTrend <- function(){

  newuiui <-
    shinydashboard::dashboardPage(
      dashboardHeader(
        title = paste0("Quick prototype of trendseq exploration app"),
        titleWidth = 900),

      dashboardSidebar(
        width = 350,
        menuItem("App settings",icon = icon("cogs"),
                 textInput("inText1","Type something")



        ),
        menuItem("Plot settings", icon = icon("paint-brush"),
                 # bstooltip: Use the widgets below to setup general parameters for exporting produced plots
                 numericInput("export_width",label = "Width of exported figures (cm)",value = 30,min = 2),
                 numericInput("export_height",label = "Height of exported figures (cm)",value = 30,min = 2)

        )
      ),

      dashboardBody(
        tabBox(
          width=12,

          tabPanel("About",
                   includeMarkdown(system.file("extdata", "about.md",package = "pcaExplorer"))),

          tabPanel("Instructions",
                   includeMarkdown(system.file("extdata", "instructions.md",package = "pcaExplorer"))),

          tabPanel("Data Preview",
                   h1("Here will go some head of the count dataset, the samples design/covariates and so"),
                   DT::dataTableOutput("inspectMatrix")),

          tabPanel("Main View",
                   p(h1('MAIN VIEW')),
                   selectInput("geneInput","type the name of a gene",choices = c("",genenames),selected = NULL),
                   selectInput("kdInput","type the name of a knockdown condition",choices = c("",alltheknockdowns),selected = NULL),

                   h2("AFFECTED GENES"),
                   textOutput("printAffectedGenes"),

                   h2("INVOLVED CONDITIONS"),
                   textOutput("printConditionsInvolved")



                   )

        )

      ),

      skin="blue"
    )

  ## ------------------------------------------------------------------ ##
  ##                          Define server                             ##
  ## ------------------------------------------------------------------ ##

  newserver <- function(input, output) {



    output$inspectMatrix <- renderDataTable({
      head(proto_kdmat,20)
    })

    affectedGenes <- reactive({
      allgenes_selectedKD <- proto_kdmat[,colnames(proto_kdmat) %in% input$kdInput,drop=FALSE] # force staying as a matrix
      affectedGenes <- rownames(allgenes_selectedKD)[allgenes_selectedKD==1]
      return(affectedGenes)
    })

    output$printAffectedGenes <- renderPrint({
            cat(paste0("The following genes are affected in the",input$kdInput," condition:\n",
                 paste(affectedGenes(), collapse=", ")))

    })

    conditionsInvolved <- reactive({
      allconds_selectedGene <- proto_kdmat[(rownames(proto_kdmat) %in% input$geneInput),,drop=FALSE] # force staying as a matrix
      involvedConds <- colnames(allconds_selectedGene)[allconds_selectedGene==1]
      return(involvedConds)
    })


    output$printConditionsInvolved <- renderPrint({

      cat(paste0("The following conditions are involved in the regulation of the ",input$geneInput," gene:\n",
                 paste(conditionsInvolved(), collapse=", ")))

    })



  }
  shinyApp(ui = newuiui, server = newserver)
}
