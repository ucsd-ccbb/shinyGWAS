#libraries
library(shiny)
library(igvShiny)
library(htmlwidgets)
library(shinyWidgets)
library(data.table)
library(dplyr)
library(qqman)
library(plotly)
# setwd("/Users/adammark/projects/shiny/shinyGWAS/app")
source("CircosFunctions.R")
source("data_prep.R")
getwd()
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
                               downloadButton('downFile',"Save Plot")),
                        column(width = 12,
                               plotOutput("plot2")))
             ),
             tabPanel("Manhattan",
                      fluidRow(
                        #default manhattan plot
                        # plotOutput("manhttanPlot_Default"), #TESTING
                        br(),
                        
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
                        
                        #loads IGV tracks 
                        
                        igvShinyOutput('igvShiny_tracks'),
                        br(),
                        br()
                      )
             ),
             tabPanel("Circos",
                      fluidRow(
                        column(width = 12,
                               plotOutput("circos", height = 800))),
                      fluidRow(
                        column(width = 4,
                               selectInput("chromosome", "Chromosome",
                                           c(seq(1,22), "X", "Y", "MT")),
                               downloadButton('downLoadCircos',"Save Plot")),
                      )
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
  
  output$circos <- renderPlot({
    plotCircosByChr(paste0("chr", input$chromosome), dataList)
  })
  
  output$downLoadCircos <- downloadHandler(
    
    filename = function() {
      paste0("chr", input$chromosome, "_circos_", gsub("-", "", Sys.Date()), sep=".pdf")
    },
    content = function(file){
      pdf(file)
      plotCircosByChr(paste0("chr", input$chromosome), dataList)
      dev.off()
    }
  )
  
  
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
    maxCADDval <- max(CADD_scores_df$CADD, na.rm=TRUE)
    # loadBedTrack(session, id="igvShiny_tracks", trackName="CADD Score", tbl=CADD_scores_df, color="red", deleteTracksOfSameName=FALSE);
    loadBedGraphTrack(session, 
                      id="igvShiny_tracks",
                      trackName="CADD Score", 
                      tbl=CADD_scores_df, 
                      color="random",  
                      # max=maxCADDval,
                      autoscale=TRUE, 
                      trackHeight = 75,
                      deleteTracksOfSameName=FALSE);
  })
  
  #RDB Catalog track
  observeEvent(input$addRDBScoreTrackButton, {
    print("____Adding RDB Score Track ____")
    loadBedTrack(session, 
                 id="igvShiny_tracks",
                 trackName="RDB Score", 
                 tbl=RDB_score, 
                 color="random", 
                 deleteTracksOfSameName=FALSE);
  })
  
  #Remove Tracks 
  observeEvent(input$removeUserTracksButton, {
    removeUserAddedTracks(session, id="igvShiny_tracks")
    removeUserAddedTracks(session, id="selectEQTLtissueTrack")
    updateCheckboxInput(session, inputId = "selectEQTLtissueTrack", value = FALSE) #resets checkboxes
  }) 
  
  #######################################################
  #manhattan plot (not IGV) (TESTING)
  # highlightSNPs <- subset(new_gwas, P< 1E-5)
  # output$manhttanPlot_Default <- renderPlot({
  #   manhttanPlot <-manhattan(new_gwas, chr="CHR", bp="BP", snp="SNP", p="P" ,
  #                            col=c("gray", "lightgray"),
  #                            highlight = highlightSNPs$SNP,
  #                            annotatePval = TRUE,
  #                            cex.axis = 1,
  #                            ylim=c(-log10(max(new_gwas$P)), -log10(min(new_gwas$P)))
  #                           )
  # })

  ######################################################################################
  # ## ggplot manhattan plot (TESTING)
  # #compute cumulative position of SNP
  output$manhttanPlot_Default <- renderPlot({
    don <- new_gwas %>%
      # Compute chromosome size
      group_by(CHR) %>%
      summarise(chr_len=max(POS)) %>%
      # Calculate cumulative position of each chromosome
      mutate(tot=cumsum(chr_len)-chr_len) %>%
      dplyr::select(-chr_len) %>%
    #   # Add this info to the initial dataset
      left_join(new_gwas, ., by=c("CHR"="CHR")) %>%
    #   # Add a cumulative position of each SNP
      arrange(CHR, POS) %>%
      mutate( Chromosome=POS+tot)
    # prepare x-axis
    axisdf = don %>% group_by(CHR) %>% summarize(center=( max(Chromosome) + min(Chromosome) ) / 2 )
    #plot
    ggplot(don, aes(x=Chromosome, y=-log10(P))) +
      # Show all points
      geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1.3) +
      scale_color_manual(values = rep(c("grey", "darkgray"), 22 )) +
      # custom X axis:
      scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
      scale_y_continuous(expand = c(0, 0) )     +# remove space between plot area and x axis
      # Custom the theme:
      theme_bw() +
      theme(
        legend.position="none",
        panel.border = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank()
      )
  })
  
  
  ############### MANHATTAN via PLOTLY -----T E S T I N G ------- #####################3
  # Prepare the dataset
  # don <- new_gwas %>% 
  #   # Compute chromosome size
  #   group_by(CHR) %>% 
  #   summarise(chr_len=max(BP)) %>% 
  #   # Calculate cumulative position of each chromosome
  #   mutate(tot=cumsum(chr_len)-chr_len) %>%
  #   select(-chr_len) %>%
  #   # Add this info to the initial dataset
  #   left_join(new_gwas, ., by=c("CHR"="CHR")) %>%
  #   # Add a cumulative position of each SNP
  #   arrange(CHR, BP) %>%
  #   mutate( BPcum=BP+tot) %>%
  #   # Add highlight and annotation information
  #   mutate( is_highlight=ifelse(SNP %in% snpsOfInterest, "yes", "no")) %>%
  #   # Filter SNP to make the plot lighter
  #   filter(-log10(P)>0.5)
  # # Prepare X axis
  # axisdf <- don %>% group_by(CHR) %>% summarize(center=( max(BPcum) + min(BPcum) ) / 2 )
  # # Prepare text description for each SNP:
  # don$text <- paste("SNP: ", don$SNP, "\nPosition: ", don$BP, "\nChromosome: ", don$CHR, "\nLOD score:", -log10(don$P) %>% round(2), "\nWhat else do you wanna know", sep="")
  # # Make the plot
  # p <- ggplot(don, aes(x=BPcum, y=-log10(P), text=text)) +
  #   # Show all points
  #   geom_point( aes(color=as.factor(CHR)), alpha=0.8, size=1.3) +
  #   scale_color_manual(values = rep(c("grey", "darkgray"), 22 )) +
  #   # custom X axis:
  #   scale_x_continuous( label = axisdf$CHR, breaks= axisdf$center ) +
  #   scale_y_continuous(expand = c(0, 0) ) +     # remove space between plot area and x axis
  #   ylim(0,9) +
  #   # Add highlighted points
  #   geom_point(data=subset(don, is_highlight=="yes"), color="orange", size=2) +
  #   # Custom the theme:
  #   theme_bw() +
  #   theme( 
  #     legend.position="none",
  #     panel.border = element_blank(),
  #     panel.grid.major.x = element_blank(),
  #     panel.grid.minor.x = element_blank()
  #   )
  # 
  # output$manhttanPlot_Default <- renderPlot({
  #   ggplotly(p, tooltip="text")
  # })
  ########### end of plot ################
  
  
  
  
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
    #CADD score track
    if (x[1]=="CADD" || x[1]=="RDB"){
      x[3]<-"rsID"
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
  
}

# Create Shiny app ----
shinyApp(ui, server)
