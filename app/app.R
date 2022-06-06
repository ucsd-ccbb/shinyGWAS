# RNA Editing Shiny App for plotting and subsetting

source("data_prep.R")
source("PlottingFunctions.R")

# Define UI for data upload app ----
ui <- fluidPage(
  navbarPage("GWAS",
             tabPanel("Summary",
                      fluidRow(
                        column(width = 12,
                          plotOutput("plot"))),
                      fluidRow(
                        column(width = 4,
                               selectInput("plotType", "Plot Type",
                                           c("minChrState" = "minChrState",
                                             "CADD Score" = "CADD",
                                             "negLogP" = "negLogP")),
                               # actionButton(
                               #   inputId = "plot_now",
                               #   label = "Plot"
                               # ),
                               downloadButton('downFile',"Save Plot")),
                        column(width = 8,
                          plotOutput("plot2")))
             ),
             tabPanel("Manhattan",
                      h4("Zoom into a locus by selecting a window and double clicking."),
                             fluidRow(
                               column(width = 10,
                                      plotOutput("manhattan", height = 600,
                                                 dblclick = "manhattan_dblclick",
                                                 brush = brushOpts(
                                                   id = "manhattan_brush",
                                                   resetOnNew = TRUE
                                                 )
                                      )
                               ),
                             )
             ),
             tabPanel("Circos",
                      verbatimTextOutput("circos")
             )
             ))
  
# Define server logic to read selected file ----
server <- function(input, output) {
  p1 <- ggplot(fuma_snps_df, aes(factor(GenomicLocus), fill=func)) + geom_bar(position="stack") + theme_classic() + xlab("Genomic Loci") + theme(axis.text.x = element_text(size=10, angle=90))
  p2 <- ggplot(fuma_snps_df, aes(func, CADD, fill=func)) + geom_boxplot() + theme_classic() + theme(axis.text.x = element_text(angle=90)) + xlab("")
  p3 <- ggplot(fuma_snps_df, aes(func, minChrState, fill=func)) + geom_boxplot() + theme_classic() + theme(axis.text.x = element_text(angle=90)) + xlab("")
  p4 <- ggplot(fuma_snps_df, aes(func, negLogP, fill=func)) + geom_boxplot() + theme_classic() + theme(axis.text.x = element_text(angle=90)) + xlab("")
  
  # plotnow <- eventReactive(input$plot_now, {
  # 
  # })
  
  myPlot <- function() {

    if(input$plotType == "CADD"){
      p2
    }
    else if(input$plotType == "minChrState") {
      p3
    }
    
    else if(input$plotType == "negLogP") {
      p4
    }
  }
  
  output$plot <- renderPlot({
    p1
  })
  output$plot2 <- renderPlot({
    myPlot()
  })
  
  output$downFile <- downloadHandler(
    
    filename = function() {
      paste0(input$plotType, "_", gsub("-", "", Sys.Date()), sep=".pdf")
    },
    content = function(file){
      pdf(file)
      print(myPlot())
      dev.off()
    }
  )
  output$manhattan <- renderPlot({
    gg.manhattan(gwas, threshold=1e-6, hlight=NULL, col=mypalette, xlims=ranges$x, ylims=c(0,9), title="")
  })
  
  # When a double-click happens, check if there's a brush on the plot.
  # If so, zoom to the brush bounds; if not, reset the zoom.
  ranges <- reactiveValues(x = NULL, y = NULL)
  observeEvent(input$manhattan_dblclick, {
    brush <- input$manhattan_brush
    if (!is.null(brush)) {
      ranges$x <- c(brush$xmin, brush$xmax)
      ranges$y <- c(brush$ymin, brush$ymax)
      
    } else {
      ranges$x <- NULL
      ranges$y <- NULL
    }
  })
  
  output$table1 <- DT::renderDataTable({
    DT::datatable(fuma_snps_df)
  })
  # 
  # output$table2 <- DT::renderDataTable({
  #   DT::datatable(fuma_eqtl)
  # })
  # 
  # output$table3 <- DT::renderDataTable({
  #   DT::datatable(fuma_genes)
  # })
  
}

# Create Shiny app ----
shinyApp(ui, server)