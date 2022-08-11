#libraries
library(shiny)
library(igvShiny)
library(htmlwidgets)
library(shinyWidgets)
library(data.table)
# library(dplyr)
# setwd("/Users/adammark/projects/shiny/shinyGWAS/app")
source("scripts/data_prep.R")
source("scripts/CircosFunctions.R")
source("scripts/CircosPlots.R")
source("scripts/GWASSummaryBarPlot.R")
# we need a local directory to write files - for instance, a vcf file representing a genomic
# region of interest.  we then tell shiny about that directory, so that shiny's built-in http server
# can serve up files we write there, ultimately consumed by igv.js
if(!dir.exists("tracks"))
  dir.create("tracks")
addResourcePath("tracks", "tracks")

# Define UI for data upload app ----
ui <- fluidPage(
  navbarPage("GWAS",
             tabPanel("Summary",
                      fluidRow(
                        column(width = 12,
                               plotOutput("plot"))),
                      #GWAS Summary bar plot
                      hr(),
                      fluidRow(
                        column(width = 12,
                               height="700px",
                               plotOutput("GwasCatSummBarPlot", height =700))),
                      hr(),
                      fluidRow(
                        column(width = 4,
                               selectInput("plotType", "Plot Type",
                                           c("minChrState" = "minChrState",
                                             "CADD Score" = "CADD",
                                             "negLogP" = "negLogP"))
                        ),
                        column(width = 4,
                               selectInput("gl", "Genomic Loci",
                                           c("All", sort(unique(fuma_snps_df$GenomicLocus))))),
                        column(width = 4,
                               style = "margin-top: 25px;",
                               downloadButton('downFilePdf',"Save .pdf"),
                               downloadButton('downFilePng',"Save .png")),
                        column(width = 12,
                               plotOutput("plot2")))
             ),
             tabPanel("Manhattan",
                      fluidRow(
                        #add Tracks tracks
                        actionButton("addGWASTrackButton", "Add GWAS Track"),
                        actionButton("addGWASCatalogTrackButton", "Add GWAS Catalog Track"),
                        dropdownButton(
                          inputId = "mydropdown",
                          label="eQTL Tissues",
                          icon=NULL,
                          inline = TRUE,
                          circle=FALSE,
                          
                          ##drop down menu
                          selectInput("selectEQTLtissueTrack", "Add eQTL Tissue track",
                                      lapply(1:length(eQTL_tissues), function(i) {
                                        choices= eQTL_tissues[i] =  eQTL_tissues[i]})
                          )),
                        actionButton("addCAADScoreTrackButton", "Add CADD Score Track"), 
                        actionButton("addRDBScoreTrackButton", "Add RDB Score Track"), 
                        # remove tracks
                        actionButton("removeUserTracksButton", "Remove All Tracks"),
                        br(),
                        
                        #lads IGV tracks 
                        igvShinyOutput('igvShiny_tracks'),
                        br(),
                        br()
                      )
             ),
             tabPanel("Circos",
                      fluidRow(
                        column(width = 12,
                               plotOutput("circos", width=600, height = 600))),
                      fluidRow(
                        column(width = 4,
                               selectInput("chromosome", "Chromosome",
                                           # c(seq(1,22), "X", "Y", "MT")),
                                           names(dataList$chrFileList)),
                               # seq(1,22)[seq(1,22) %in% as.numeric(gsub("^chr","", names(dataList[["chrFileList"]])))]
                               
                               downloadButton('downLoadCircosPdf',"Save .pdf"),
                               downloadButton('downLoadCircosPng',"Save .png"))
                      )
             ),
             tabPanel("Session Info",
                      verbatimTextOutput("sessionInfo")
                      )
  )
)

# Define server logic to read selected file ----
server <- function(input, output, session) {
  #### ----------------Summary  Tab -----------------------------
  
  p1 <- ggplot(fuma_snps_df, aes(factor(GenomicLocus), fill=func)) + geom_bar(position="stack") + theme_classic() + xlab("Genomic Loci") + theme(axis.text.x = element_text(size=10, angle=90))
  
  myPlot <- function() {
    if(input$gl == "All"){
      fuma_plot_data <- fuma_snps_df
    }else{
      fuma_plot_data <- subset(fuma_snps_df, GenomicLocus == input$gl)
    }
    p2 <- ggplot(fuma_plot_data, aes(func, CADD, fill=func)) + geom_boxplot() + theme_classic() + theme(axis.text.x = element_text(angle=90)) + xlab("")
    p3 <- ggplot(fuma_plot_data, aes(func, minChrState, fill=func)) + geom_boxplot() + theme_classic() + theme(axis.text.x = element_text(angle=90)) + xlab("")
    p4 <- ggplot(fuma_plot_data, aes(func, negLogP, fill=func)) + geom_boxplot() + theme_classic() + theme(axis.text.x = element_text(angle=90)) + xlab("")
    
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
  
  output$downFilePdf <- downloadHandler(
    filename = function() {
      paste0(input$plotType, "_", gsub("-", "", Sys.Date()), sep=".pdf")
    },
    content = function(file){
      pdf(file)
      print(myPlot())
      dev.off()
    }
  )
  
  output$downFilePng <- downloadHandler(
    filename = function() {
      paste0(input$plotType, "_", gsub("-", "", Sys.Date()), sep=".png")
    },
    content = function(file){
      png(file, width = 2200, height = 1500, units = "px", res = 300)
      print(myPlot())
      dev.off()
    }
  )
  
  output$manhattan <- renderPlot({
    gg.manhattan(gwas, threshold=1e-6, hlight=NULL, col=mypalette, xlims=ranges$x, ylims=c(0,9), title="")
  })
  
  
  
  output$circos <- renderPlot({
    # print(paste0("chr", input$chromosome))
    plotCircosByChr(input$chromosome, dataList)
 
  })
  
  output$downLoadCircosPdf <- downloadHandler(

    filename = function() {
      paste0(input$chromosome, "_circos_", gsub("-", "", Sys.Date()), sep=".pdf")
    },
    content = function(file){
      plotCircosByChr(input$chromosome, dataList)
      
    
      dev.off()
    }
  )
  
  output$downLoadCircosPng <- downloadHandler(
    
    filename = function() {
      paste0("chr", input$chromosome, "_circos_", gsub("-", "", Sys.Date()), sep=".png")
    },
    content = function(file){
      png(file, width = 3000, height = 3000, units = "px", res = 300)
      plotCircosByChr(paste0("chr", input$chromosome), dataList)
      dev.off()
    }
  )
  
  #GWAS Catalag Summary Bar plot
  output$GwasCatSummBarPlot <- renderPlot({
    gwasBarplot
  })
  
 
  
  
  #### ----------------Manhattan Tab -----------------------------
  #add tracks
  #manhattan plot
  observeEvent(input$addGWASTrackButton, {
    print("____Adding GWAS Track ____")
    showGenomicRegion(session, id="igvShiny_tracks", "all")
    loadGwasTrack(session, 
                  id="igvShiny_tracks", 
                  trackName="GWAS", 
                  tbl=new_gwas, 
                  ymin= -log10(max(new_gwas$P)),
                  ymax= -log10(min(new_gwas$P)),
                  deleteTracksOfSameName=FALSE)
  })
  
  #GWAS Catalog track
  observeEvent(input$addGWASCatalogTrackButton, {
    print("____Adding GWAS Catalog ____")
    loadBedTrack(session, id="igvShiny_tracks", trackName="GWAS Catalog", tbl=gwasCatalog, color="random", deleteTracksOfSameName=FALSE);
  })
  head(fuma_eqtl)
  #eQTL tissue track
  observeEvent(input$selectEQTLtissueTrack,{
    print("____Adding eQTL Track ____")
    tissue <- input$selectEQTLtissueTrack
    tissue<-tissue[length(tissue)]
    index <- which(fuma_eqtl$tissue==tissue)
    loadGwasTrack(session, 
                  id="igvShiny_tracks", 
                  trackName=paste0("eQTL Tissue: ", tissue), 
                  tbl=fuma_eqtl[index, ],
                  deleteTracksOfSameName=FALSE)
    # updateCheckboxInput(session, inputId = "selectEQTLtissueTrack", value = FALSE) #resets checkboxes
    
  })
  
  #CADD Score track
  observeEvent(input$addCAADScoreTrackButton, {
    print("____Adding CADD Score Track ____")
    maxCADDval <- max(CADD_scores$CADD, na.rm=TRUE)
    # loadBedTrack(session, id="igvShiny_tracks", trackName="CADD Score", tbl=CADD_scores_df, color="red", deleteTracksOfSameName=FALSE);
    loadBedGraphTrack(session, 
                      id="igvShiny_tracks",
                      trackName="CADD Score", 
                      tbl=CADD_scores, 
                      color="random",  
                      # max=maxCADDval,
                      autoscale=TRUE, 
                      trackHeight = 75,
                      deleteTracksOfSameName=FALSE);
  })
  
  
  #RDB Catalog track
  observeEvent(input$addRDBScoreTrackButton, {
    print("____Adding RDB Score Track ____")
    loadBedTrack(session, id="igvShiny_tracks",trackName="RDB Score", tbl=RDB_score, color="blue", deleteTracksOfSameName=FALSE);
    # loadBedGraphTrack(session, 
    #                   id="igvShiny_tracks",
    #                   trackName="RegulomeDB score", 
    #                   tbl=RDB_score, 
    #                   color="random",  
    #                   # max=maxCADDval,
    #                   autoscale=TRUE, 
    #                   trackHeight = 75,
    #                   deleteTracksOfSameName=FALSE);
  })
  
  #Remove Tracks 
  observeEvent(input$removeUserTracksButton, {
    removeUserAddedTracks(session, id="igvShiny_tracks")
    removeUserAddedTracks(session, id="selectEQTLtissueTrack")
    updateCheckboxInput(session, inputId = "selectEQTLtissueTrack", value = FALSE) #resets checkboxes
  }) 
  
  #manhattan plot & GWAS Catalog (this will load by default)
  observeEvent(input$igvReady, {
    containerID <- input$igvReady
    showGenomicRegion(session, id="igvShiny_tracks", "all")
    loadGwasTrack(session, 
                  id="igvShiny_tracks", 
                  trackName="GWAS", 
                  tbl=new_gwas, 
                  ymin= -log10(max(new_gwas$P)),
                  ymax= -log10(min(new_gwas$P)),
                  deleteTracksOfSameName=TRUE)
    loadBedTrack(session, 
                 id="igvShiny_tracks", 
                 trackName="GWAS Catalog",
                 tbl=gwasCatalog,
                 color="random", 
                 deleteTracksOfSameName=TRUE);  })
  
  # Pop up info box
  observeEvent(input$trackClick, { #add popup window when a SNP is clicked
    printf("--- igv-trackClick popup")
    x <- input$trackClick
    print(length(x))
    print(x)
    maxCols=length(x)/2

    if (length(x)>14 && x[15]=="GTEX"){ #add GTEX protal hyperlink
      x[16] <- paste0("<a href='",x[16],"'>", x[4]," URL</a>")
    }
    #
    if (length(x)>8 &&x[9]=="PUBMEDID"){ #add hyperlink to pubmed
      x[10]<- paste0("<a href='https://pubmed.ncbi.nlm.nih.gov/",x[10],"'>", x[10],"</a>")
      x[11] <- "Chr"
    }
    #RBD score track
    if ( x[1]=="RDB"){
      x[3]<-"RDB score Info"
    }

    attribute.name.positions <- grep("name", names(x))
    attribute.value.positions <- grep("value", names(x))
    attribute.names <- as.character(x)[attribute.name.positions][1:maxCols] #if different SNPs points overlap, one click will list several of them
    attribute.values <- as.character(x)[attribute.value.positions][1:maxCols]

    tbl <- data.frame(name=attribute.names,
                      value=attribute.values,
                      stringsAsFactors=FALSE)

    dialogContent <- renderTable(tbl,
                                 striped=TRUE,
                                 hover=TRUE,
                                 aliged="c",
                                 bordered = TRUE,
                                 width="100%",
                                 sanitize.text.function = function(x) x)

    html <- HTML(dialogContent())
    showModal(modalDialog(html,
                          size="m",
                          easyClose=TRUE))
  })

  
  # render IGV 
  output$igvShiny_tracks <- renderIgvShiny({
    cat("--- starting renderIgvShiny\n");
    genomeOptions <- parseAndValidateGenomeSpec(genomeName="hg38",  initialLocus="all")
    igvShiny(genomeOptions)
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
  
  output$sessionInfo <- renderPrint({
    #capture.output(sessionInfo())
    sessionInfo()
  })
  
}

# Create Shiny app ----
shinyApp(ui, server)